#!/bin/bash

########################################
############# CSCI 2951-O ##############
########################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
ORTOOLS_PREFIX="$(brew --prefix or-tools)"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$ORTOOLS_PREFIX;/opt/homebrew"

cmake --build "$BUILD_DIR" --config Release
