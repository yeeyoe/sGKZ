import Foundation
import Combine
import AppKit
import UniformTypeIdentifiers

@MainActor
final class AppState: ObservableObject {
    @Published var rows: [PointRow] = Array(repeating: PointRow(), count: 8)
    @Published var showLabels = true
    @Published var showLattice = true
    @Published var showAxes = true
    @Published var solveState: SolveState = .idle
    @Published var plotState = PlotState()
    @Published var status = "Enter at least three integer points"

    private var solving = false
    private var cancellationToken: GKZCancellationToken?
    private var solveGeneration = 0
    private var undoStack: [[PointRow]] = []
    private var redoStack: [[PointRow]] = []
    private var activeEditKey: String?
    private var lastPointSignature = ""

    var canUndo: Bool { !undoStack.isEmpty }
    var canRedo: Bool { !redoStack.isEmpty }

    var validPoints: [CGPoint] { rows.compactMap(\.coordinatePoint) }

    func normalizeRows() {
        let currentPointSignature = pointSignature
        let pointsChanged = currentPointSignature != lastPointSignature
        if pointsChanged {
            invalidateSolve()
            plotState = PlotState()
            for index in rows.indices { rows[index].sigma = "" }
            solveState = .idle
            status = "Not calculated"
        }
        lastPointSignature = currentPointSignature
        var next = 1
        for index in rows.indices {
            rows[index].number = rows[index].validCoordinates ? next : nil
            if rows[index].validCoordinates { next += 1 }
            if !rows[index].validCoordinates { rows[index].sigma = "" }
        }
        if rows.last?.validCoordinates == true { rows.append(PointRow()) }
    }

    func update(row: PointRow, at index: Int) {
        rows[index] = row
        normalizeRows()
    }

    func setCoordinate(at index: Int, axis: CoordinateAxis, value: String) {
        guard rows.indices.contains(index) else { return }
        let key = "\(index)-\(axis.rawValue)"
        if activeEditKey != key {
            recordUndo()
            activeEditKey = key
        }
        switch axis {
        case .x: rows[index].x = value
        case .y: rows[index].y = value
        }
        normalizeRows()
    }

    func commitTextEdit() {
        activeEditKey = nil
    }

    func deleteRow(at index: Int) {
        guard rows.indices.contains(index) else { return }
        recordUndo()
        activeEditKey = nil
        rows.remove(at: index)
        if rows.isEmpty { rows.append(PointRow()) }
        normalizeRows()
    }

    func undo() {
        guard let previous = undoStack.popLast() else { return }
        redoStack.append(rows)
        activeEditKey = nil
        rows = previous
        normalizeRows()
    }

    func redo() {
        guard let next = redoStack.popLast() else { return }
        undoStack.append(rows)
        activeEditKey = nil
        rows = next
        normalizeRows()
    }

    func solve() {
        guard !solving else { return }
        let points = validPoints
        guard points.count >= 3 else { status = "At least three valid points are required"; solveState = .failed(status); return }
        guard rows.filter(\.validCoordinates).allSatisfy({ $0.integerPoint != nil }) else {
            status = "Exact sGKZ requires integer coordinates"
            solveState = .failed(status)
            return
        }
        var seen = Set<String>()
        for point in points {
            let key = "\(Int(point.x)),\(Int(point.y))"
            if !seen.insert(key).inserted { status = "Duplicate point: \(key)"; solveState = .failed(status); return }
        }
        let origin = points[0]
        let direction = CGPoint(x: points[1].x - origin.x, y: points[1].y - origin.y)
        let hasArea = points.dropFirst(2).contains { point in
            direction.x * (point.y - origin.y) - direction.y * (point.x - origin.x) != 0
        }
        guard hasArea else { status = "All points are collinear; hull area is zero"; solveState = .failed(status); return }
        solving = true
        solveGeneration += 1
        let generation = solveGeneration
        let token = GKZCancellationToken()
        cancellationToken = token
        solveState = .running("Solving exactly")
        status = "Exact verification in progress"
        let boxed = points.map { NSValue(point: NSPoint(x: $0.x, y: $0.y)) }
        DispatchQueue.global(qos: .userInitiated).async { [weak self] in
            let result = GKZBridge.solvePoints(boxed, progress: { iterations, active, text in
                Task { @MainActor [weak self] in
                    guard let self, self.solveGeneration == generation, !token.isCancelled else { return }
                    self.solveState = .running("\(text) (iteration \(iterations), active set \(active))")
                }
            }, token: token)
            Task { @MainActor [weak self] in self?.apply(result, generation: generation) }
        }
    }

    func exportCSV() {
        let panel = NSSavePanel()
        panel.allowedContentTypes = [.commaSeparatedText]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = "GKZ-points.csv"
        guard panel.runModal() == .OK, let url = panel.url else { return }

        var lines = ["No.,x,y,sGKZ"]
        for row in rows where row.validCoordinates {
            let values = [row.number.map(String.init) ?? "", row.x, row.y, row.sigma]
            lines.append(values.map(csvField).joined(separator: ","))
        }
        do {
            try (lines.joined(separator: "\r\n") + "\r\n").write(to: url, atomically: true, encoding: .utf8)
            status = "CSV exported"
        } catch {
            status = "CSV export failed: \(error.localizedDescription)"
        }
    }

    func cancelSolve() {
        cancellationToken?.cancel()
        solving = false
        cancellationToken = nil
        solveGeneration += 1
        solveState = .idle
        status = "Cancelled"
    }

    private func invalidateSolve() {
        if solving { cancellationToken?.cancel() }
        solving = false
        cancellationToken = nil
        solveGeneration += 1
    }

    private func apply(_ result: GKZGuiResult, generation: Int) {
        guard generation == solveGeneration else { return }
        solving = false
        cancellationToken = nil
        if result.cancelled {
            status = "Cancelled"
            solveState = .idle
            return
        }
        guard result.certified else {
            status = result.message
            solveState = .failed(result.message)
            return
        }
        var validIndex = 0
        for index in rows.indices where rows[index].validCoordinates {
            if validIndex < result.sigmaExact.count { rows[index].sigma = result.sigmaExact[validIndex] }
            validIndex += 1
        }
        plotState = PlotState(sigmaVee: result.sigmaVeeExact,
            subdivision: result.subdivisionCells.map { cell in cell.map { value in
                let point = value.pointValue
                return CGPoint(x: point.x, y: point.y)
            } },
                              triangulation: result.triangulationFaces.map { $0.map(\.intValue) })
        status = "Exact verification passed"
        solveState = .success
    }

    /// Update a point during canvas dragging without rebuilding the row model.
    /// Keeping this operation small prevents SwiftUI from interrupting the mouse drag.
    func beginPointDrag() {
        recordUndo()
        activeEditKey = nil
        invalidateSolve()
        plotState = PlotState()
        for index in rows.indices { rows[index].sigma = "" }
        solveState = .idle
        status = "Not calculated"
    }

    func movePoint(at validIndex: Int, to point: CGPoint) {
        var current = 0
        for index in rows.indices where rows[index].validCoordinates {
            if current == validIndex {
                var row = rows[index]
                row.x = String(format: "%.8g", point.x)
                row.y = String(format: "%.8g", point.y)
                row.sigma = ""
                rows[index] = row
                lastPointSignature = pointSignature
                return
            }
            current += 1
        }
    }

    private func recordUndo() {
        guard undoStack.last != rows else { return }
        undoStack.append(rows)
        if undoStack.count > 10 { undoStack.removeFirst() }
        redoStack.removeAll()
    }

    private var pointSignature: String {
        rows.compactMap(\.coordinatePoint)
            .map { "\($0.x.bitPattern),\($0.y.bitPattern)" }
            .joined(separator: ";")
    }

    private func csvField(_ value: String) -> String {
        guard value.contains(",") || value.contains("\"") || value.contains("\n") || value.contains("\r") else {
            return value
        }
        return "\"\(value.replacingOccurrences(of: "\"", with: "\"\""))\""
    }
}

enum CoordinateAxis: String {
    case x, y
}
