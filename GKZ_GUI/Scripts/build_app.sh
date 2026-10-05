#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
if ! command -v xcodebuild >/dev/null 2>&1 || ! xcodebuild -version >/dev/null 2>&1; then
  echo "Full Xcode is required to build the SwiftUI app; Command Line Tools alone are insufficient." >&2
  exit 1
fi
if [ -z "${GKZ_DEPENDENCY_PREFIX:-}" ] && ! command -v brew >/dev/null 2>&1; then
  echo "A developer dependency prefix is required (CGAL, Eigen, GMP, MPFR). Set GKZ_DEPENDENCY_PREFIX." >&2
  exit 1
fi
DEP_PREFIX=${GKZ_DEPENDENCY_PREFIX:-$(brew --prefix)}
ARCHS=${GKZ_ARCHS:-$(uname -m)}
export GKZ_DEPENDENCY_PREFIX="$DEP_PREFIX"
export GKZ_ARCHS="$ARCHS"
"$ROOT/Scripts/build_dependencies.sh"
xcodebuild -project "$ROOT/GKZ_GUI.xcodeproj" -scheme GKZ_GUI -configuration Release \
  -derivedDataPath "$ROOT/build/xcode" ARCHS="$ARCHS" \
  MACOSX_DEPLOYMENT_TARGET=13.0 GKZ_DEPENDENCY_PREFIX="$DEP_PREFIX" build

APP_PATH="$ROOT/build/xcode/Build/Products/Release/sGKZ.app"
rm -f "$APP_PATH/Contents/Frameworks/libgmpxx.dylib" \
  "$APP_PATH/Contents/Frameworks/libgmp.dylib" \
  "$APP_PATH/Contents/Frameworks/libmpfr.dylib"
codesign --force --deep --sign - "$APP_PATH"
"$ROOT/Scripts/verify_bundle.sh"
