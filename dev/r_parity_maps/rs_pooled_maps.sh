#!/usr/bin/env bash
# rs_pooled_maps.sh — the dada2-rs side of r_pooled_maps.R: lay out an existing
# dada-pooled run for compare_maps.py (#277).
#
# Usage:
#   rs_pooled_maps.sh <binary> <dada-dir> <out-dir> <fastq>...
# <dada-dir> is dada-pooled's --output-dir (<s>.json or <s>.json.gz per sample).
#
# Per sample <s>, in <out-dir>:
#   <s>.rs.derep.json  `derep` of the FASTQ: the uniques the sample's `map` indexes
#   <s>.rs.dada.json   the pooled per-sample JSON, decompressed
# compare_maps.py checks that map + derep reproduce every ASV's reads, which a
# derep in a different order than dada-pooled's own would fail.
set -euo pipefail
BIN=${1:?binary}; DADA=${2:?dada-dir}; OUT=${3:?out-dir}; shift 3
[ $# -gt 0 ] || { echo "usage: rs_pooled_maps.sh <binary> <dada-dir> <out-dir> <fastq>..." >&2; exit 1; }
mkdir -p "$OUT"

for fq in "$@"; do
  # dada-pooled names its output by the FASTQ stem, _filt included; the
  # comparison keys on the name without it, as r_pooled_maps.R does.
  stem=$(basename "$fq" .fastq.gz); s=${stem%_filt}
  if [ -f "$DADA/$stem.json.gz" ]; then
    gzip -dc "$DADA/$stem.json.gz" > "$OUT/$s.rs.dada.json"
  elif [ -f "$DADA/$stem.json" ]; then
    cp "$DADA/$stem.json" "$OUT/$s.rs.dada.json"
  else
    echo "rs_pooled_maps.sh: no $stem.json[.gz] in $DADA" >&2; exit 1
  fi
  "$BIN" derep "$fq" -o "$OUT/$s.rs.derep.json" 2> "$OUT/$s.rs.derep.log"
  echo "$s: done"
done
