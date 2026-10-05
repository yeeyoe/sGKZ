import SwiftUI

@main
struct GKZ_GUIApp: App {
    @Environment(\.openWindow) private var openWindow

    var body: some Scene {
        WindowGroup("sGKZ") { ContentView() }
            .commands {
                CommandGroup(replacing: .appInfo) {
                    Button("About sGKZ") { openWindow(id: "about") }
                }
            }

        Window("About sGKZ", id: "about") {
            AboutView()
                .frame(width: 340, height: 190)
        }
        .windowResizability(.contentSize)
    }
}

private struct AboutView: View {
    var body: some View {
        VStack(spacing: 10) {
            Text("sGKZ")
                .font(.system(size: 24, weight: .semibold))
            Text("Version: 1.0")
            Text("Vibe Coding by Yi Yao 姚懿")
            Link("github.com/yeeyoe", destination: URL(string: "https://github.com/yeeyoe")!)
            Text("Date: 2026-10-05")
        }
        .multilineTextAlignment(.center)
        .padding(24)
    }
}
