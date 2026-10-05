#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
APP_PATH="$ROOT/build/xcode/Build/Products/Release/sGKZ.app"
"$ROOT/Scripts/verify_bundle.sh" "$APP_PATH"
mkdir -p "$ROOT/dist"
STAGING="$ROOT/build/dmg-staging"
rm -rf "$STAGING"
mkdir -p "$STAGING"
ditto "$APP_PATH" "$STAGING/sGKZ.app"
ln -sfn /Applications "$STAGING/Applications"
ARCHS=${GKZ_ARCHS:-$(lipo -archs "$APP_PATH/Contents/MacOS/sGKZ")}
case "$ARCHS" in
  "arm64 x86_64"|"x86_64 arm64") IMAGE_NAME=sGKZ-macOS-universal.dmg ;;
  arm64) IMAGE_NAME=sGKZ-macOS-arm64.dmg ;;
  x86_64) IMAGE_NAME=sGKZ-macOS-x86_64.dmg ;;
  *) echo "Unsupported app architecture: $ARCHS" >&2; exit 1 ;;
esac
hdiutil create -volname sGKZ -srcfolder "$STAGING" -ov -format UDZO "$ROOT/dist/$IMAGE_NAME"
echo "Created $ROOT/dist/$IMAGE_NAME"
