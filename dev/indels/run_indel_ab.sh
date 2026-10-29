#!/usr/bin/env bash
# run_indel_ab.sh — A/B of the indel-aware λ (issue #291) on a full run, under
# dada-pooled: `default` (substitution-only λ, R's) against `indel`
# (--indel-rates from learn-errors --indel-rates-out).
#
# Both arms use one error model: learn-errors with --indel-rates-out learns the
# substitution model exactly as without it (byte-identical JSON) and the indel
# rates after it. Only λ differs between the arms, through merge, table and
# chimera removal.
#
# Usage:
#   run_indel_ab.sh <binary> <out-dir> <manifest.tsv>
#
# The manifest is tab-separated, no header: `sample  fwd.fastq.gz  [rev.fastq.gz]`.
# Two columns is single-end (PacBio); three is paired (Illumina), denoised per
# direction and merged.
#
# Environment:
#   THREADS       default 8
#   EXTRA         flags for learn-errors and dada-pooled, e.g. PacBio's
#                 "--band 32 --homo-gap-p -1 --kmer-size 7"
#   LEARN_EXTRA   flags for learn-errors only, e.g. "--errfun pacbio"
#   INDEL_EXTRA   flags for the indel arm's dada-pooled only, e.g.
#                 "--indel-end-zone 0" (it requires --indel-rates)
#   MERGE_EXTRA   flags for merge-pairs
#   FIT_CMD       --indel-fit-cmd for learn-errors (default: the native fit)
#   ARMS          default "default indel"
#   SKIP_LEARN    1 reuses <out-dir>/learn/ from an earlier run
#
# Writes <out-dir>/learn/ (error models, indel rates and history) and
# <out-dir>/<arm>/ (per-sample JSON, seqtab.json, nochim.json, logs, and
# time_*.txt with wall time and peak memory). Then:
#   python3 dev/indels/indel_ab_compare.py --binary <binary> <out-dir>
set -euo pipefail

BIN="${1:?usage: run_indel_ab.sh <binary> <out-dir> <manifest.tsv>}"
OUT="${2:?missing out-dir}"
MAN="${3:?missing manifest.tsv}"
THREADS="${THREADS:-8}"
ARMS="${ARMS:-default indel}"
read -r -a extra <<< "${EXTRA:-}"
read -r -a learn_extra <<< "${LEARN_EXTRA:-}"
read -r -a indel_extra <<< "${INDEL_EXTRA:-}"
read -r -a merge_extra <<< "${MERGE_EXTRA:-}"

samples=(); fwd=(); rev=()
while IFS=$'\t' read -r s f r; do
  [ -z "$s" ] && continue
  samples+=("$s"); fwd+=("$f"); [ -n "${r:-}" ] && rev+=("$r")
done < "$MAN"
n=${#samples[@]}
[ "$n" -gt 0 ] || { echo "no samples in $MAN" >&2; exit 2; }
paired=0
if [ "${#rev[@]}" -gt 0 ]; then
  [ "${#rev[@]}" -eq "$n" ] || { echo "$MAN: some rows lack a reverse file" >&2; exit 2; }
  paired=1
fi
dirs=(fwd); [ "$paired" = 1 ] && dirs=(fwd rev)
names=$(IFS=,; echo "${samples[*]}")
# (Plain loops and copies throughout: macOS ships bash 3.2, without mapfile.)
set_inputs() { if [ "$1" = fwd ]; then in=("${fwd[@]}"); else in=("${rev[@]}"); fi; }

# Wall time and peak memory: GNU time on Linux, BSD time on macOS; untimed
# (with a note in the log) where /usr/bin/time is missing.
timed() {
  local log=$1; shift
  if /usr/bin/time -v true 2> /dev/null; then /usr/bin/time -v -o "$log" "$@"
  elif /usr/bin/time -l true 2> /dev/null; then /usr/bin/time -l "$@" 2> >(tee "$log" >&2)
  else echo "no /usr/bin/time: untimed" > "$log"; "$@"; fi
}

mkdir -p "$OUT/learn"
cp "$MAN" "$OUT/manifest.tsv"
if [ "${SKIP_LEARN:-0}" != 1 ]; then
  for d in "${dirs[@]}"; do
    fit=(); [ -n "${FIT_CMD:-}" ] && fit=(--indel-fit-cmd "$FIT_CMD")
    echo "==> learn-errors $d ($n input(s))"
    set_inputs "$d"
    timed "$OUT/learn/time_$d.txt" "$BIN" learn-errors "${in[@]}" --threads "$THREADS" \
      ${extra[@]+"${extra[@]}"} ${learn_extra[@]+"${learn_extra[@]}"} \
      -o "$OUT/learn/$d.err.json" --indel-rates-out "$OUT/learn/$d.indel_rates.tsv" \
      --indel-history "$OUT/learn/$d.indel_history.tsv" ${fit[@]+"${fit[@]}"} \
      --verbose 2> "$OUT/learn/$d.log"
    grep "indel iter" "$OUT/learn/$d.log" | tail -1
  done
fi

for arm in $ARMS; do
  o="$OUT/$arm"
  mkdir -p "$o"
  for d in "${dirs[@]}"; do
    lam=(); [ "$arm" = indel ] && lam=(--indel-rates "$OUT/learn/$d.indel_rates.tsv" ${indel_extra[@]+"${indel_extra[@]}"})
    mkdir -p "$o/$d"
    echo "==> $arm dada-pooled $d"
    set_inputs "$d"
    timed "$o/time_$d.txt" "$BIN" dada-pooled "${in[@]}" --sample-names "$names" \
      --error-model "$OUT/learn/$d.err.json" --output-dir "$o/$d" --threads "$THREADS" \
      ${extra[@]+"${extra[@]}"} ${lam[@]+"${lam[@]}"} --verbose 2> "$o/$d.log"
    # Refuse an arm whose output does not record the λ it was asked for: arms
    # that silently ran alike would look exactly like a null result.
    first=$(ls "$o/$d"/*.json | grep -v '/_' | head -1)
    has=$(python3 -I -c "import json,sys; print('indel_lambda' in json.dumps(json.load(open(sys.argv[1]))))" "$first")
    want=False; [ "$arm" = indel ] && want=True
    [ "$has" = "$want" ] || { echo "ERROR: $first: indel_lambda recorded=$has, expected $want" >&2; exit 1; }
  done
  fj=(); rj=()
  for s in "${samples[@]}"; do fj+=("$o/fwd/$s.json"); rj+=("$o/rev/$s.json"); done
  if [ "$paired" = 1 ]; then
    echo "==> $arm merge-pairs"
    "$BIN" merge-pairs --fwd-dada "${fj[@]}" --rev-dada "${rj[@]}" --fwd-fastq "${fwd[@]}" \
      --rev-fastq "${rev[@]}" --sample-names "$names" --threads "$THREADS" \
      ${merge_extra[@]+"${merge_extra[@]}"} -o "$o/merged.json" 2> "$o/merge.log"
    "$BIN" make-sequence-table "$o/merged.json" -o "$o/seqtab.json" 2> /dev/null
  else
    "$BIN" make-sequence-table "${fj[@]}" -o "$o/seqtab.json" 2> /dev/null
  fi
  echo "==> $arm remove-bimera-denovo"
  "$BIN" remove-bimera-denovo "$o/seqtab.json" --method consensus --threads "$THREADS" \
    -o "$o/nochim.json" 2> "$o/chimera.log"
done
echo "done: $OUT"
