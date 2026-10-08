#!/usr/bin/env bash
# install_r_trace.sh — build R DADA2 with the #277 per-read trace (r_trace.patch)
# into its own library, for trace_read.sh. The user's DADA2 install is not touched.
#
# Usage:
#   install_r_trace.sh <lib-dir> [dada2-source-dir]
# Without a source dir, DADA2 at REF (default ab84434, = 1.36.0 on GitHub) is
# downloaded. Build with the same R that produced the reference you compare
# against: compiler flags are part of what is being tested.
set -euo pipefail
LIB=${1:?lib-dir}; SRC=${2:-}
REF=${REF:-ab84434}
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

if [ -n "$SRC" ]; then
  cp -R "$SRC" "$WORK/dada2"
else
  curl -fsSL "https://github.com/benjjneb/dada2/archive/$REF.tar.gz" | tar -xz -C "$WORK"
  mv "$WORK"/dada2-* "$WORK/dada2"
fi
rm -f "$WORK"/dada2/src/*.o "$WORK"/dada2/src/*.so
patch -d "$WORK/dada2" -p1 < "$HERE/r_trace.patch"
mkdir -p "$LIB"
R CMD INSTALL --no-test-load -l "$LIB" "$WORK/dada2"
echo "installed traced DADA2 into $LIB with $(R --version | head -n 1)"
echo "compiler flags: $(R CMD config CXX) $(R CMD config CXXFLAGS)"
