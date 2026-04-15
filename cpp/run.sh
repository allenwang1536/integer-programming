#!/bin/bash

########################################
############# CSCI 2951-O ##############
########################################
set -euo pipefail

E_BADARGS=65
if [ $# -ne 1 ]
then
	echo "Usage: `basename $0` <input>"
	exit $E_BADARGS
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT="$1"
BIN="$SCRIPT_DIR/build/solver"

if [ ! -x "$BIN" ]; then
    "$SCRIPT_DIR/compile.sh"
fi

"$BIN" "$INPUT"
