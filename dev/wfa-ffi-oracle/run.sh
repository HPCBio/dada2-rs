#!/usr/bin/env bash
# Build the C++ WFA2-lib at a pinned commit and run the #102 oracle against it.
# Nothing here ships; this is a throwaway test oracle (see README.md).
set -euo pipefail
cd "$(dirname "$0")"

# Pinned to the commit tested when results were recorded in the README. Bump
# deliberately (and re-record results) to re-check a newer WFA2-lib.
WFA_COMMIT="${WFA_COMMIT:-bcf473a}"
SRC=wfa2lib_src

if [ ! -d "$SRC" ]; then
  git clone https://github.com/smarco/WFA2-lib "$SRC"
fi
git -C "$SRC" fetch -q origin && git -C "$SRC" checkout -q "$WFA_COMMIT"

# The Makefile expects these output dirs to exist.
mkdir -p "$SRC/build/cpp" "$SRC/lib"
make -C "$SRC" lib_wfa >/dev/null

cc -Wall -I"$SRC" oracle.c "$SRC/lib/libwfa.a" -o oracle -lm -lpthread
echo "WFA2-lib commit: $(git -C "$SRC" rev-parse --short HEAD)"
./oracle
