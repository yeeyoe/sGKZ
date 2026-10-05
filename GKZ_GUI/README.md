# GKZ_GUI

Native macOS SwiftUI/AppKit front end for the two-dimensional point-mode shortest GKZ solver.

The app accepts integer points in a table, renders the convex hull and lattice points, supports zoom/pan/drag with integer snapping, and displays exact sGKZ values only after rational certification. Certified subdivision data is drawn over the point set.

Build the core adapter with CMake, then open the Xcode project on a full macOS development installation. The intended release target is Universal 2 with a macOS 13 deployment target. GMP and MPFR are linked statically into the app; CGAL and Eigen are header-only at runtime.

`Scripts/build_dependencies.sh` builds the C++ libraries for the selected `GKZ_ARCHS` (defaulting to the host architecture). `Scripts/build_app.sh` requires a full Xcode installation and development headers/libraries for the selected architecture. `Scripts/make_dmg.sh` creates an architecture-labelled drag-to-Applications disk image after `Scripts/verify_bundle.sh` confirms there are no developer-machine paths or dynamically linked GMP/MPFR dependencies in the app.

The current source tree can compile and type-check the C++/Objective-C++ layers with Command Line Tools, but producing the `.app` requires selecting full Xcode with `sudo xcode-select --switch /Applications/Xcode.app/Contents/Developer` (or the installed Xcode path). A Homebrew arm64-only prefix builds an Apple Silicon app only; build or supply arm64+x86_64 dependencies to create a Universal 2 release.
