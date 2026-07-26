#!/usr/bin/env bash
# Drive the whole LOESS oracle comparison from one learn-errors JSON.
#
# Usage:
#   ./run.sh <learn_errors.json> [outdir]
#
# Produces, in outdir (default ./out):
#   triples.tsv          the (q, rlogp, tot) inputs, full 40-column grid
#   r_pred.tsv           R stats::loess, both surfaces
#   rust_pred.tsv        dada2-rs + loess-rs arms
#   n{3..8}.tsv          anchors-only sparsity sweep (n_valid = 3..8)
#   and prints the comparison tables
#
# Requires: R with no extra packages (stats only), python3, cargo.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
json="${1:?usage: run.sh <learn_errors.json> [outdir]}"
out="${2:-$here/out}"
mkdir -p "$out"

if ! command -v Rscript >/dev/null; then
  echo "Rscript not found — the ground-truth arm needs R (stats only)." >&2
  exit 1
fi

echo "==> building harness"
cargo build --quiet --manifest-path "$here/Cargo.toml"
bin="$here/target/debug/loess-oracle"

echo "==> full grid"
python3 "$here/extract_triples.py" "$json" "$out/triples.tsv"
Rscript "$here/r_truth.R" "$out/triples.tsv" "$out/r_pred.tsv" >/dev/null
"$bin" "$out/triples.tsv" "$out/rust_pred.tsv"
python3 "$here/compare.py" "$out/r_pred.tsv" "$out/rust_pred.tsv" --floor

echo
echo "==> production-shaped binned case (counts only at q=12,24,38)"
python3 "$here/extract_triples.py" "$json" "$out/triples_i100.tsv" --bins 12,24,38
Rscript "$here/r_truth.R" "$out/triples_i100.tsv" "$out/r_i100.tsv" >/dev/null 2>&1 || true
"$bin" "$out/triples_i100.tsv" "$out/rust_i100.tsv"
python3 "$here/compare.py" "$out/r_i100.tsv" "$out/rust_i100.tsv" --floor

echo
echo "==> sparsity sweep (anchors only, n_valid = 3..8)"
for n in 3 4 5 6 7 8; do
  python3 "$here/extract_triples.py" "$json" "$out/n$n.tsv" \
    --nbins "$n" --anchors-only >/dev/null
  Rscript "$here/r_truth.R" "$out/n$n.tsv" "$out/rn$n.tsv" >/dev/null 2>&1 || true
  "$bin" "$out/n$n.tsv" "$out/sn$n.tsv" 2>/dev/null
  echo "--- n_valid = $n"
  python3 "$here/compare.py" "$out/rn$n.tsv" "$out/sn$n.tsv" --floor | tail -n +3
done

echo
echo "outputs in $out"
