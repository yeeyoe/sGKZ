import SwiftUI

struct ToolbarView: View {
    @ObservedObject var state: AppState

    var body: some View {
        HStack(spacing: 14) {
            Button("sGKZ") { state.solve() }
                .buttonStyle(.borderedProminent)
                .disabled({ if case .running = state.solveState { return true }; return false }())
            Button { state.undo() } label: {
                Image(systemName: "arrow.uturn.backward")
            }
            .accessibilityLabel("Undo")
            .keyboardShortcut("z", modifiers: .command)
            .disabled(!state.canUndo)
            Button { state.redo() } label: {
                Image(systemName: "arrow.uturn.forward")
            }
            .accessibilityLabel("Redo")
            .keyboardShortcut("z", modifiers: [.command, .shift])
            .disabled(!state.canRedo)
            Button { state.exportCSV() } label: {
                Image(systemName: "square.and.arrow.down")
            }
            .accessibilityLabel("Export CSV")
            .help("Export CSV")
            if case .running = state.solveState {
                Button("Cancel") { state.cancelSolve() }
            }
            Toggle("Labels", isOn: $state.showLabels)
            Toggle("Lattice points", isOn: $state.showLattice)
            Toggle("Axes", isOn: $state.showAxes)
            Spacer()
            Text(state.status).foregroundStyle(.secondary).lineLimit(1)
        }
        .padding(.horizontal, 12).padding(.vertical, 8)
    }
}
