import Foundation

enum SolveState: Equatable {
    case idle
    case running(String)
    case success
    case failed(String)
}
