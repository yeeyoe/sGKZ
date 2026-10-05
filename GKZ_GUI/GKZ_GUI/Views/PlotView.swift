import SwiftUI
import AppKit

struct PlotView: NSViewRepresentable {
    @ObservedObject var state: AppState

    func makeNSView(context: Context) -> GKZPlotNSView { GKZPlotNSView(state: state) }
    func updateNSView(_ view: GKZPlotNSView, context: Context) {
        view.state = state
        view.needsDisplay = true
    }
}

private struct PlotTransform {
    let scale: CGFloat
    let center: CGPoint
    let offset: CGSize
    let bounds: CGRect

    func screen(_ point: CGPoint) -> CGPoint {
        CGPoint(x: bounds.midX + (point.x - center.x) * scale + offset.width,
                y: bounds.midY + (point.y - center.y) * scale + offset.height)
    }

    func world(_ point: CGPoint) -> CGPoint {
        CGPoint(x: center.x + (point.x - bounds.midX - offset.width) / scale,
                y: center.y + (point.y - bounds.midY - offset.height) / scale)
    }
}

final class GKZPlotNSView: NSView {
    var state: AppState
    private var zoom: CGFloat = 1
    private var offset = CGSize.zero
    private var draggingIndex: Int?
    private var lastPoint: NSPoint?
    private var dragTransform: PlotTransform?
    private var dragStartPoint: CGPoint?
    private var dragStartScreen: CGPoint?
    private var dragPreviewPoint: CGPoint?
    private var pointDragStarted = false
    private var viewportScale: CGFloat?
    private var viewportCenter = CGPoint.zero
    private var viewportSize: CGSize = .zero
    private var lastPointSignature = ""
    private var didFitInitialData = false

    init(state: AppState) {
        self.state = state
        super.init(frame: .zero)
        wantsLayer = true
        layer?.backgroundColor = NSColor.windowBackgroundColor.cgColor
    }

    required init?(coder: NSCoder) { fatalError() }

    override var acceptsFirstResponder: Bool { true }

    override func acceptsFirstMouse(for event: NSEvent?) -> Bool { true }

    override func draw(_ dirtyRect: NSRect) {
        guard let context = NSGraphicsContext.current?.cgContext else { return }
        context.saveGState()
        context.clip(to: bounds)
        defer { context.restoreGState() }
        context.setFillColor(NSColor.windowBackgroundColor.cgColor)
        context.fill(bounds)

        var points = state.validPoints
        if let index = draggingIndex, points.indices.contains(index), let preview = dragPreviewPoint {
            points[index] = preview
        }
        let validRows = state.rows.filter(\.validCoordinates)
        let transform = makeTransform(points)
        if state.showAxes { drawAxes(context, transform) }
        if state.showLattice { drawLattice(context, transform) }
        if !state.plotState.subdivision.isEmpty { drawSubdivision(context, transform) }
        drawHull(context, transform, points)
        for (index, point) in points.enumerated() {
            let screen = transform.screen(point)
            context.setFillColor(NSColor.systemBlue.cgColor)
            context.fillEllipse(in: CGRect(x: screen.x - 4, y: screen.y - 4, width: 8, height: 8))
            if state.showLabels, index < validRows.count, let number = validRows[index].number {
                let numberLabel = "\(number)"
                let coordinateLabel = " (\(coordinateString(point.x)), \(coordinateString(point.y)))"
                let numberAttributes: [NSAttributedString.Key: Any] = [
                    .font: NSFont.boldSystemFont(ofSize: 11),
                    .foregroundColor: NSColor.systemRed
                ]
                let coordinateAttributes: [NSAttributedString.Key: Any] = [
                    .font: NSFont.systemFont(ofSize: 11),
                    .foregroundColor: NSColor.labelColor
                ]
                NSString(string: numberLabel).draw(
                    at: CGPoint(x: screen.x + 7, y: screen.y + 5),
                    withAttributes: numberAttributes)
                let numberWidth = NSString(string: numberLabel).size(withAttributes: numberAttributes).width
                NSString(string: coordinateLabel).draw(
                    at: CGPoint(x: screen.x + 7 + numberWidth, y: screen.y + 5),
                    withAttributes: coordinateAttributes)
            }
        }
    }

    private func makeTransform(_ points: [CGPoint]) -> PlotTransform {
        let signature = pointSignature(points)
        let needsFit = viewportScale == nil
            || viewportSize != bounds.size
            || (!didFitInitialData && !points.isEmpty)
        if needsFit {
            if points.isEmpty {
                // Keep the empty canvas at a modest, useful scale until the first points arrive.
                viewportScale = 48 * zoom
                viewportCenter = .zero
            } else {
                let xs = points.map(\.x), ys = points.map(\.y)
                let minX = (xs.min() ?? -1) - 1, maxX = (xs.max() ?? 1) + 1
                let minY = (ys.min() ?? -1) - 1, maxY = (ys.max() ?? 1) + 1
                let width = max(maxX - minX, 1)
                let height = max(maxY - minY, 1)
                // Use the smaller ratio so neither dimension can run outside the canvas.
                let fitScale = min(bounds.width / width, bounds.height / height) * 0.72
                viewportScale = max(0.0001, min(fitScale, 80) * zoom)
                viewportCenter = CGPoint(x: (minX + maxX) / 2, y: (minY + maxY) / 2)
                didFitInitialData = true
            }
            viewportSize = bounds.size
            lastPointSignature = signature
            offset = .zero
        }
        let scale = viewportScale ?? 1
        return PlotTransform(scale: scale, center: viewportCenter, offset: offset, bounds: bounds)
    }

    private func pointSignature(_ points: [CGPoint]) -> String {
        points.map { "\($0.x),\($0.y)" }.joined(separator: ";")
    }

    private func fitTransform(_ points: [CGPoint]) -> PlotTransform {
        let xs = points.map(\.x), ys = points.map(\.y)
        let minX = (xs.min() ?? -1) - 1, maxX = (xs.max() ?? 1) + 1
        let minY = (ys.min() ?? -1) - 1, maxY = (ys.max() ?? 1) + 1
        let scale = max(0.0001, min(bounds.width / max(maxX - minX, 1),
                                    bounds.height / max(maxY - minY, 1)) * 0.8 * zoom)
        return PlotTransform(scale: scale,
                             center: CGPoint(x: (minX + maxX) / 2, y: (minY + maxY) / 2),
                             offset: offset, bounds: bounds)
    }

    private func drawAxes(_ context: CGContext, _ transform: PlotTransform) {
        let origin = transform.screen(.zero)
        let xAxisY = min(max(origin.y, bounds.minY + 18), bounds.maxY - 18)
        let yAxisX = min(max(origin.x, bounds.minX + 34), bounds.maxX - 4)
        context.setStrokeColor(NSColor.separatorColor.cgColor)
        context.setLineWidth(1)
        context.move(to: CGPoint(x: bounds.minX, y: xAxisY))
        context.addLine(to: CGPoint(x: bounds.maxX, y: xAxisY))
        context.move(to: CGPoint(x: yAxisX, y: bounds.minY))
        context.addLine(to: CGPoint(x: yAxisX, y: bounds.maxY))
        context.strokePath()

        let step = tickStep(for: transform.scale)
        let xRange = visibleWorldRange(horizontal: true, transform: transform)
        let yRange = visibleWorldRange(horizontal: false, transform: transform)
        context.setStrokeColor(NSColor.tertiaryLabelColor.withAlphaComponent(0.55).cgColor)
        context.setLineWidth(1)
        if let xRange {
            for tick in integerTicks(xRange, step: step) {
                let p = transform.screen(CGPoint(x: tick, y: 0))
                context.move(to: CGPoint(x: p.x, y: xAxisY - 4))
                context.addLine(to: CGPoint(x: p.x, y: xAxisY + 4))
            }
        }
        if let yRange {
            for tick in integerTicks(yRange, step: step) {
                let p = transform.screen(CGPoint(x: 0, y: tick))
                context.move(to: CGPoint(x: yAxisX - 4, y: p.y))
                context.addLine(to: CGPoint(x: yAxisX + 4, y: p.y))
            }
        }
        context.strokePath()

        let attributes: [NSAttributedString.Key: Any] = [
            .font: NSFont.monospacedDigitSystemFont(ofSize: 10, weight: .regular),
            .foregroundColor: NSColor.secondaryLabelColor
        ]
        if let xRange {
            for tick in integerTicks(xRange, step: step) {
                guard tick != 0 else { continue }
                let p = transform.screen(CGPoint(x: tick, y: 0))
                let label = String(tick) as NSString
                label.draw(at: CGPoint(x: p.x - label.size(withAttributes: attributes).width / 2,
                                       y: xAxisY + 6), withAttributes: attributes)
            }
        }
        if let yRange {
            for tick in integerTicks(yRange, step: step) {
                guard tick != 0 else { continue }
                let p = transform.screen(CGPoint(x: 0, y: tick))
                let label = String(tick) as NSString
                label.draw(at: CGPoint(x: yAxisX + 6, y: p.y - label.size(withAttributes: attributes).height / 2),
                           withAttributes: attributes)
            }
        }
    }

    private func tickStep(for scale: CGFloat) -> Int {
        // Keep individual integer ticks readable at the normal zoom level.
        if scale >= 24 { return 1 }
        let desired = 58 / max(scale, 0.0001)
        let magnitude = pow(10, floor(log10(max(desired, 1))))
        for factor in [1.0, 2.0, 5.0, 10.0] where factor * magnitude >= desired {
            return max(1, Int(min(factor * magnitude, 1_000_000_000)))
        }
        return max(1, Int(min(10 * magnitude, 1_000_000_000)))
    }

    private func visibleWorldRange(horizontal: Bool, transform: PlotTransform) -> ClosedRange<Double>? {
        let a = horizontal ? transform.world(CGPoint(x: bounds.minX, y: bounds.midY)).x
                           : transform.world(CGPoint(x: bounds.midX, y: bounds.minY)).y
        let b = horizontal ? transform.world(CGPoint(x: bounds.maxX, y: bounds.midY)).x
                           : transform.world(CGPoint(x: bounds.midX, y: bounds.maxY)).y
        guard a.isFinite, b.isFinite else { return nil }
        return min(a, b)...max(a, b)
    }

    private func integerTicks(_ range: ClosedRange<Double>, step: Int) -> [Int] {
        let stride = Double(step)
        let first = ceil(range.lowerBound / stride) * stride
        let last = floor(range.upperBound / stride) * stride
        guard first.isFinite, last.isFinite, first >= Double(Int.min), last <= Double(Int.max),
              last >= first, (last - first) / stride <= 2_000 else { return [] }
        return Swift.stride(from: Int(first.rounded()), through: Int(last.rounded()), by: step).map { $0 }
    }

    private func drawHull(_ context: CGContext, _ transform: PlotTransform, _ points: [CGPoint]) {
        let hull = convexHull(points)
        guard hull.count > 1 else { return }
        context.setStrokeColor(NSColor.systemBlue.cgColor)
        context.setLineWidth(2)
        context.move(to: transform.screen(hull[0]))
        for point in hull.dropFirst() { context.addLine(to: transform.screen(point)) }
        context.closePath()
        context.strokePath()
    }

    private func drawLattice(_ context: CGContext, _ transform: PlotTransform) {
        guard let xRange = visibleWorldRange(horizontal: true, transform: transform),
              let yRange = visibleWorldRange(horizontal: false, transform: transform),
              xRange.lowerBound >= Double(Int.min), xRange.upperBound <= Double(Int.max),
              yRange.lowerBound >= Double(Int.min), yRange.upperBound <= Double(Int.max) else { return }
        let x0 = Int(ceil(xRange.lowerBound)), x1 = Int(floor(xRange.upperBound))
        let y0 = Int(ceil(yRange.lowerBound)), y1 = Int(floor(yRange.upperBound))
        guard x1 >= x0, y1 >= y0,
              Double(x1) - Double(x0) <= 100_000,
              Double(y1) - Double(y0) <= 100_000,
              (Double(x1) - Double(x0) + 1) * (Double(y1) - Double(y0) + 1) <= 100_000 else { return }
        context.setFillColor(NSColor.systemGray.withAlphaComponent(0.65).cgColor)
        for x in x0...x1 {
            for y in y0...y1 {
                let p = transform.screen(CGPoint(x: x, y: y))
                context.fillEllipse(in: CGRect(x: p.x - 2, y: p.y - 2, width: 4, height: 4))
            }
        }
    }

    private func drawSubdivision(_ context: CGContext, _ transform: PlotTransform) {
        context.setStrokeColor(NSColor.systemOrange.cgColor)
        context.setLineWidth(1.5)
        for cell in state.plotState.subdivision where cell.count > 1 {
            context.move(to: transform.screen(cell[0]))
            for point in cell.dropFirst() { context.addLine(to: transform.screen(point)) }
            context.closePath()
            context.strokePath()
        }
    }

    override func scrollWheel(with event: NSEvent) {
        let factor: CGFloat = event.scrollingDeltaY > 0 ? 1.1 : 0.9
        let oldZoom = zoom
        zoom = min(max(zoom * factor, 0.2), 20)
        if oldZoom != zoom {
            viewportScale = max(0.0001, (viewportScale ?? 1) * (zoom / oldZoom))
            needsDisplay = true
        }
    }

    override func mouseDown(with event: NSEvent) {
        window?.makeFirstResponder(self)
        let point = convert(event.locationInWindow, from: nil)
        lastPoint = point
        let points = state.validPoints
        let transform = makeTransform(points)
        dragTransform = transform
        draggingIndex = points.indices.min {
            distance(transform.screen(points[$0]), point) < distance(transform.screen(points[$1]), point)
        }
        if let index = draggingIndex {
            if distance(transform.screen(points[index]), point) > 14 { draggingIndex = nil }
            else {
                dragStartPoint = points[index]
                dragStartScreen = point
                dragPreviewPoint = points[index]
                pointDragStarted = false
            }
        }
    }

    override func mouseDragged(with event: NSEvent) {
        guard let previous = lastPoint else { return }
        let current = convert(event.locationInWindow, from: nil)
        let delta = CGSize(width: current.x - previous.x, height: current.y - previous.y)
        if let index = draggingIndex, let transform = dragTransform,
           let start = dragStartPoint, let startScreen = dragStartScreen {
            let totalDelta = CGSize(width: current.x - startScreen.x, height: current.y - startScreen.y)
            let worldDelta = CGPoint(x: totalDelta.width / transform.scale, y: totalDelta.height / transform.scale)
            let updated = state.showLattice
                ? CGPoint(x: (start.x + worldDelta.x).rounded(), y: (start.y + worldDelta.y).rounded())
                : CGPoint(x: start.x + worldDelta.x, y: start.y + worldDelta.y)
            if updated != start && !pointDragStarted {
                state.beginPointDrag()
                pointDragStarted = true
            }
            dragPreviewPoint = updated
        } else {
            offset.width += delta.width
            offset.height += delta.height
        }
        lastPoint = current
        needsDisplay = true
    }

    override func mouseUp(with event: NSEvent) {
        if pointDragStarted, let index = draggingIndex, let finalPoint = dragPreviewPoint {
            state.movePoint(at: index, to: finalPoint)
        }
        draggingIndex = nil
        dragStartPoint = nil
        dragStartScreen = nil
        dragPreviewPoint = nil
        pointDragStarted = false
        dragTransform = nil
        lastPoint = nil
        needsDisplay = true
    }

    private func coordinateString(_ value: CGFloat) -> String {
        value.rounded() == value ? String(Int(value)) : String(format: "%.6g", Double(value))
    }

    private func distance(_ a: CGPoint, _ b: CGPoint) -> CGFloat { hypot(a.x - b.x, a.y - b.y) }

    private func convexHull(_ points: [CGPoint]) -> [CGPoint] {
        let sorted = points.sorted { $0.x == $1.x ? $0.y < $1.y : $0.x < $1.x }
        func cross(_ o: CGPoint, _ a: CGPoint, _ b: CGPoint) -> CGFloat {
            (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
        }
        var lower = [CGPoint]()
        for point in sorted {
            while lower.count >= 2 && cross(lower[lower.count - 2], lower[lower.count - 1], point) <= 0 { lower.removeLast() }
            lower.append(point)
        }
        var upper = [CGPoint]()
        for point in sorted.reversed() {
            while upper.count >= 2 && cross(upper[upper.count - 2], upper[upper.count - 1], point) <= 0 { upper.removeLast() }
            upper.append(point)
        }
        return Array(lower.dropLast() + upper.dropLast())
    }

    private func insideOrOnBoundary(_ point: CGPoint, _ polygon: [CGPoint]) -> Bool {
        guard polygon.count >= 3 else { return false }
        var sign: CGFloat = 0
        for index in polygon.indices {
            let a = polygon[index], b = polygon[(index + 1) % polygon.count]
            let dx = b.x - a.x, dy = b.y - a.y
            let cross = dx * (point.y - a.y) - dy * (point.x - a.x)
            let tolerance = max(1, hypot(dx, dy)) * 1e-9
            if abs(cross) <= tolerance { continue }
            let current: CGFloat = cross < 0 ? -1 : 1
            if sign == 0 { sign = current }
            else if current != sign { return false }
        }
        return true
    }
}
