#!/usr/bin/env bash
# run_member_order.sh — run one dataset under every DADA2RS_MEMBER_ORDER arm
# (issue #157), for dev/compare_member_order.py.
#
# The error model is PINNED: the gate also acts inside learn-errors (it runs
# through run_dada), so learning per arm would vary the model as well as the
# order. Learn it once, without the gate, and pass it here.
#
# Usage:
#   run_member_order.sh <binary> <err.json> <out-dir> <input> [<input> ...]
#
# Environment:
#   MODE     pooled (default) | pseudo | per-sample
#   ARMS     default "insertion sorted shuffle:1 shuffle:2 shuffle:3 shuffle:4 shuffle:5"
#   THREADS  default 8
#   EXTRA    extra dada flags, e.g. the platform's "--band 32 --kmer-size 7"
#
# Writes <out-dir>/<arm>/ (":" -> "_") with the per-sample JSON and log.txt
# (the --verbose log, which compare_member_order.py reads for birth order).
# Timing is not the point, so no NUMA pinning; arms can run as separate jobs.
#
# Then:
#   python3 dev/compare_member_order.py <out-dir>/insertion <out-dir>/sorted \
#       <out-dir>/shuffle_{1..5}
set -euo pipefail

BIN="${1:?usage: run_member_order.sh <binary> <err.json> <out-dir> <input>...}"
ERR="${2:?missing err.json}"
OUT="${3:?missing out-dir}"
shift 3
[ "$#" -gt 0 ] || { echo "no inputs" >&2; exit 2; }
MODE="${MODE:-pooled}"
ARMS="${ARMS:-insertion sorted shuffle:1 shuffle:2 shuffle:3 shuffle:4 shuffle:5}"
THREADS="${THREADS:-8}"
read -r -a extra <<< "${EXTRA:-}"

for arm in $ARMS; do
  o="$OUT/${arm/:/_}"
  mkdir -p "$o"
  echo "==> $arm ($MODE, $# input(s))"
  export DADA2RS_MEMBER_ORDER="$arm"
  case "$MODE" in
    pooled)
      "$BIN" dada-pooled "$@" --error-model "$ERR" --output-dir "$o" \
        --threads "$THREADS" --verbose ${extra[@]+"${extra[@]}"} 2> "$o/log.txt" ;;
    pseudo)
      "$BIN" dada-pseudo "$@" --error-model "$ERR" --output-dir "$o" \
        --threads "$THREADS" --verbose ${extra[@]+"${extra[@]}"} 2> "$o/log.txt" ;;
    per-sample)
      : > "$o/log.txt"
      for f in "$@"; do
        s=$(basename "$f"); s=${s%.fastq.gz}; s=${s%.fastq}
        "$BIN" dada "$f" --error-model "$ERR" -o "$o/$s.json" \
          --threads "$THREADS" --verbose ${extra[@]+"${extra[@]}"} 2>> "$o/log.txt"
      done ;;
    *) echo "unknown MODE=$MODE" >&2; exit 2 ;;
  esac
  # Refuse a run whose log does not name the arm: a gate that silently did not
  # take effect would make the arms identical, which is also what a null looks like.
  if [ "$arm" != insertion ] && ! grep -q "member order=$arm (OVERRIDDEN)" "$o/log.txt"; then
    echo "ERROR: $o/log.txt does not confirm member order=$arm" >&2; exit 1
  fi
done
unset DADA2RS_MEMBER_ORDER
echo "done: $OUT"
