#!/bin/sh
set -eu
GUI_ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
PREFIX=${GKZ_DEPENDENCY_PREFIX:-$(brew --prefix)}
ARCHS=${GKZ_ARCHS:-$(uname -m)}
cmake -S "$GUI_ROOT" -B "$GUI_ROOT/build-root" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_OSX_ARCHITECTURES="$ARCHS"
cmake --build "$GUI_ROOT/build-root" --target gkz_core gkz_gui_adapter --parallel
echo "Built GKZ_GUI C++ adapter for $ARCHS using $PREFIX"
