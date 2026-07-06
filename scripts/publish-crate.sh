#!/usr/bin/env bash
# publish-crate.sh — publish dada2-rs to crates.io as a Needleman-Wunsch-only crate.
#
# The experimental WFA backend depends on `wfa2lib-rs` via a *git* dependency,
# which crates.io rejects outright (issue #63). This script temporarily strips
# that git dependency and its off-by-default `wfa` feature from Cargo.toml, runs
# `cargo publish`, then restores the original manifest. The published crate is
# therefore NW-only; WFA stays a source-checkout developer build
# (`cargo build --features wfa`).
#
# Usage:
#   scripts/publish-crate.sh --dry-run     # verify packaging without uploading
#   scripts/publish-crate.sh               # real publish (needs `cargo login`)
#
# Any extra arguments are passed straight through to `cargo publish`.
set -euo pipefail
cd "$(dirname "$0")/.."

MANIFEST="Cargo.toml"
BACKUP="$(mktemp)"

restore() {
    if [ -f "$BACKUP" ]; then
        cp "$BACKUP" "$MANIFEST"
        rm -f "$BACKUP"
        echo "==> Restored original $MANIFEST"
    fi
}
trap restore EXIT
cp "$MANIFEST" "$BACKUP"

echo "==> Stripping the experimental WFA git dependency + feature for publishing"
python3 - "$MANIFEST" <<'PY'
import re, sys
path = sys.argv[1]
s = open(path).read()
before = s
# Remove the git `wfa2lib-rs` dependency line...
s = re.sub(r'(?m)^wfa2lib-rs\s*=.*\n', '', s)
# ...and the off-by-default `wfa` feature that referenced it.
s = re.sub(r'(?m)^wfa\s*=\s*\[.*\]\s*\n', '', s)
if s == before:
    sys.exit("ERROR: expected WFA dependency/feature lines were not found; "
             "manifest layout may have changed — update publish-crate.sh.")
open(path, "w").write(s)
PY

echo "==> cargo publish --allow-dirty $*"
# --allow-dirty because we have just edited Cargo.toml in place (and callers may
# publish with other working-tree changes present).
cargo publish --allow-dirty "$@"
