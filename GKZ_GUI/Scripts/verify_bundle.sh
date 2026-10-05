#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
APP_PATH=${1:-"$ROOT/build/xcode/Build/Products/Release/sGKZ.app"}
EXECUTABLE="$APP_PATH/Contents/MacOS/sGKZ"
[ -d "$APP_PATH" ] || { echo "missing app: $APP_PATH" >&2; exit 1; }
[ -x "$EXECUTABLE" ] || { echo "missing executable: $EXECUTABLE" >&2; exit 1; }
ARCHS=$(lipo -archs "$EXECUTABLE")
EXPECTED_ARCHS=${GKZ_ARCHS:-arm64}
for arch in $EXPECTED_ARCHS; do
  case " $ARCHS " in *" $arch "*) ;; *) echo "app lacks requested architecture $arch: $ARCHS" >&2; exit 1;; esac
done

mkdir -p "$ROOT/build"
DEPENDENCIES="$ROOT/build/dependencies.txt"
otool -L "$EXECUTABLE" | sed -n 's/^[[:space:]]*\([^[:space:]]*\) (compatibility version.*$/\1/p' > "$DEPENDENCIES"
if rg -n '/opt/homebrew|/usr/local|/Users/' "$DEPENDENCIES"; then
  echo "app still refers to a developer-machine dependency path" >&2
  exit 1
fi
if rg -n 'lib(gmp|gmpxx|mpfr)' "$DEPENDENCIES"; then
  echo "math libraries must be statically linked into the app" >&2
  exit 1
fi

codesign --verify --deep --strict "$APP_PATH"
case "$ARCHS" in
  "arm64 x86_64"|"x86_64 arm64") echo "Verified Universal 2 app bundle: $APP_PATH" ;;
  *) echo "Verified $ARCHS app bundle: $APP_PATH" ;;
esac
