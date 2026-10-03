#!/usr/bin/env bash
# run_pool_tiebreak.sh — run one pooled dataset under both DADA2RS_POOL_TIEBREAK
# arms (issue #260): `lexical` (current default) and `first-seen` (R's
# combineDereps2 rule).
#
# Under `first-seen`, ties among equal-abundance pooled uniques keep their order
# of first appearance across samples, in INPUT ORDER. To compare against an R
# `dada(pool=TRUE)` baseline, pass the inputs in the order R received them, and
# use the error model R used (or the pinned one both sides were fitted from;
# check that `sum(trans)` matches first).
#
# Usage:
#   run_pool_tiebreak.sh <binary> <err.json> <out-dir> <input> [<input> ...]
#
# Environment:
#   ARMS     default "lexical first-seen"
#   THREADS  default 8
#   EXTRA    extra dada-pooled flags, e.g. the platform's "--band 32 --kmer-size 7"
#
# Writes <out-dir>/<arm>/ with the per-sample JSON, _pooled.json
# (--pooled-record), log.txt (--verbose, which compare_member_order.py reads
# for birth order) and inputs.txt (the input order this arm used).
#
# Then:
#   # arm vs arm: ASV churn and birth order
#   python3 dev/compare_member_order.py <out-dir>/lexical <out-dir>/first-seen
#   # each arm vs R (R long CSV from dev/concordance/write_reference.R)
#   python3 dev/concordance/rtab_to_seqtab.py r_dada.csv r_dada.json
#   python3 dev/compare_asvs.py --baseline R=r_dada.json \
#       --compare lexical=<out-dir>/lexical --compare first-seen=<out-dir>/first-seen
set -euo pipefail

BIN="${1:?usage: run_pool_tiebreak.sh <binary> <err.json> <out-dir> <input>...}"
ERR="${2:?missing err.json}"
OUT="${3:?missing out-dir}"
shift 3
[ "$#" -gt 0 ] || { echo "no inputs" >&2; exit 2; }
ARMS="${ARMS:-lexical first-seen}"
THREADS="${THREADS:-8}"
read -r -a extra <<< "${EXTRA:-}"

for arm in $ARMS; do
  o="$OUT/$arm"
  mkdir -p "$o"
  printf '%s\n' "$@" > "$o/inputs.txt"
  echo "==> $arm ($# input(s))"
  DADA2RS_POOL_TIEBREAK="$arm" "$BIN" dada-pooled "$@" --error-model "$ERR" \
    --output-dir "$o" --pooled-record "$o/_pooled.json" \
    --threads "$THREADS" --verbose ${extra[@]+"${extra[@]}"} 2> "$o/log.txt"
  # Refuse a run whose log does not name the arm: a gate that silently did not
  # take effect makes the arms identical, which is also what a null looks like.
  if ! grep -q "pool tiebreak=$arm " "$o/log.txt"; then
    echo "ERROR: $o/log.txt does not confirm pool tiebreak=$arm" >&2; exit 1
  fi
done
echo "done: $OUT"
