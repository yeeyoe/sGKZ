import SwiftUI

struct PointTableView: View {
    @ObservedObject var state: AppState

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            Text("Point Coordinates").font(.headline)
            ScrollView {
              Grid(alignment: .leading, horizontalSpacing: 5, verticalSpacing: 5) {
                GridRow {
                    Text("No.").frame(width: 38, alignment: .leading)
                    Text("x").frame(width: 76, alignment: .leading)
                    Text("y").frame(width: 76, alignment: .leading)
                    Text("sGKZ").frame(width: 86, alignment: .leading)
                    Color.clear.frame(width: 28)
                }.font(.caption).foregroundStyle(.secondary)
                ForEach(state.rows.indices, id: \.self) { index in
                    GridRow {
                        Text(state.rows[index].number.map(String.init) ?? "")
                            .frame(width: 38, alignment: .leading)
                        TextField("", text: Binding(get: { state.rows[index].x }, set: { state.setCoordinate(at: index, axis: .x, value: $0) }))
                            .textFieldStyle(.roundedBorder).frame(width: 76)
                            .onSubmit { state.commitTextEdit() }
                        TextField("", text: Binding(get: { state.rows[index].y }, set: { state.setCoordinate(at: index, axis: .y, value: $0) }))
                            .textFieldStyle(.roundedBorder).frame(width: 76)
                            .onSubmit { state.commitTextEdit() }
                        Text(state.rows[index].sigma)
                            .frame(width: 86, alignment: .leading)
                            .monospacedDigit()
                            .textSelection(.enabled)
                        Button {
                            state.deleteRow(at: index)
                        } label: {
                            Image(systemName: "trash")
                        }
                        .buttonStyle(.borderless)
                        .frame(width: 28)
                        .accessibilityLabel("Delete row")
                    }
                }
              }
            }
        }
        .padding(12)
    }
}
