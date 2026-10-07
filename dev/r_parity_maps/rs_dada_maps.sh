#!/usr/bin/env bash
# rs_dada_maps.sh — the dada2-rs side of r_dada_maps.R: per-sample dada under a
# fixed error model, keeping the unique -> cluster map and transitions (#277).
#
# Usage:
#   rs_dada_maps.sh <binary> <err.json> <out-dir> <fastq>...
# Env: BAND (32), KMER (5, R's compile-time KMER_SIZE), THREADS (1), EXTRA.
#
# Per sample <s>, in <out-dir>:
#   <s>.rs.derep.json  the uniques, in the order dada's `map` indexes
#   <s>.rs.dada.json   dada output with --aux-outputs (map + transitions)
# dada reads the derep JSON rather than the FASTQ so that `map` indexes a file
# we hold, instead of a derep made inside dada that is never written.
set -euo pipefail
BIN=${1:?binary}; ERR=${2:?err.json}; OUT=${3:?out-dir}; shift 3
[ $# -gt 0 ] || { echo "usage: rs_dada_maps.sh <binary> <err.json> <out-dir> <fastq>..." >&2; exit 1; }
BAND=${BAND:-32}; KMER=${KMER:-5}; THREADS=${THREADS:-1}
# shellcheck disable=SC2206
extra=(${EXTRA:-})
mkdir -p "$OUT"
"$BIN" --version

for fq in "$@"; do
  s=$(basename "$fq" .fastq.gz); s=${s%_filt}
  "$BIN" derep "$fq" -o "$OUT/$s.rs.derep.json" 2> "$OUT/$s.rs.derep.log"
  "$BIN" dada "$OUT/$s.rs.derep.json" --error-model "$ERR" --band "$BAND" --kmer-size "$KMER" \
      --threads "$THREADS" --aux-outputs ${extra[@]+"${extra[@]}"} \
      -o "$OUT/$s.rs.dada.json" 2> "$OUT/$s.rs.dada.log"
  echo "$s: done"
done
