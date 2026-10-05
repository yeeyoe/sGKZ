import Foundation

struct PointRow: Identifiable, Equatable {
    let id = UUID()
    var x = ""
    var y = ""
    var sigma = ""
    var number: Int? = nil

    var validCoordinates: Bool {
        guard let x = Double(x.trimmingCharacters(in: .whitespacesAndNewlines)),
              let y = Double(y.trimmingCharacters(in: .whitespacesAndNewlines)) else { return false }
        return x.isFinite && y.isFinite
    }

    var coordinatePoint: CGPoint? {
        guard let x = Double(x.trimmingCharacters(in: .whitespacesAndNewlines)),
              let y = Double(y.trimmingCharacters(in: .whitespacesAndNewlines)),
              x.isFinite, y.isFinite else { return nil }
        return CGPoint(x: x, y: y)
    }

    var integerPoint: CGPoint? {
        guard let point = coordinatePoint,
              point.x.rounded() == point.x, point.y.rounded() == point.y,
              abs(point.x) <= 9_007_199_254_740_991,
              abs(point.y) <= 9_007_199_254_740_991 else { return nil }
        return point
    }
}
