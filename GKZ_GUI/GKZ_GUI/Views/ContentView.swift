import SwiftUI

struct ContentView: View {
    @StateObject private var state = AppState()
    var body: some View {
        VStack(spacing: 0) {
            ToolbarView(state: state)
            Divider()
            HStack(spacing: 0) {
                PointTableView(state: state).frame(width: 370)
                Divider()
                PlotView(state: state).frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
        .frame(minWidth: 900, minHeight: 600)
        .onAppear { state.normalizeRows() }
    }
}
