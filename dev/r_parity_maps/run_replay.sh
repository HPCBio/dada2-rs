#!/usr/bin/env bash
# run_replay.sh — replay each learnErrors pass in R and dada2-rs from R's own
# per-pass input models, and report where the two first differ (#277).
#
# Usage:
#   run_replay.sh <binary> <err.rds> <out-dir> <training fastq>...
# Env: THREADS (1), PASSES (all, e.g. "1 2 6"), BAND (32), KMER (5).
#
# <err.rds> is R's learnErrors object for exactly these training FASTQs, all of
# them used (nbases above the set's total). Every pass runs per-sample dada with
# OMEGA_C = 0, as learnErrors does. For each pass k, in <out-dir>/pass_<k>/:
#   r/, rs/      both sides' per-sample maps and transitions
#   compare.txt  compare_maps.py output; diffs.tsv lists what differs
# The last pass also checks R's per-sample sums against the model's own trans;
# if that check fails, this run is not replaying learnErrors and nothing else
# it reports means anything.
set -euo pipefail
BIN=${1:?binary}; ERR=${2:?err.rds}; OUT=${3:?out-dir}; shift 3
[ $# -gt 0 ] || { echo "usage: run_replay.sh <binary> <err.rds> <out-dir> <fastq>..." >&2; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-1}; BAND=${BAND:-32}; KMER=${KMER:-5}
mkdir -p "$OUT"

Rscript "$HERE/export_err_in.R" "$ERR" "$OUT/models"
Rscript "$HERE/../../scripts/learnerrors_to_dada2rs.R" "$ERR" "$OUT/models/model.json"
n=$(ls "$OUT"/models/err_in_*.json | wc -l | tr -d ' ')
PASSES=${PASSES:-$(seq 1 "$n")}

for k in $PASSES; do
  d="$OUT/pass_$k"
  echo "==> pass $k of $n"
  Rscript "$HERE/r_dada_maps.R" "$ERR" "$d/r" "$@" --band="$BAND" --threads="$THREADS" \
      --use-err-in --err-in-iter="$k" --omega-c=0
  THREADS=$THREADS BAND=$BAND KMER=$KMER EXTRA="--omega-c 0" \
      bash "$HERE/rs_dada_maps.sh" "$BIN" "$OUT/models/err_in_$k.json" "$d/rs" "$@"
  model=()
  [ "$k" = "$n" ] && model=(--r-model "$OUT/models/model.json")
  python3 "$HERE/compare_maps.py" "$d/r" "$d/rs" -o "$d/diffs.tsv" ${model[@]+"${model[@]}"} \
      > "$d/compare.txt" 2>&1 || true
  grep -E 'model check|sample\(s\);|^==>' "$d/compare.txt"
done
