#!/usr/bin/env bash
# trace_pairs.sh — re-run each sample in a find_pairs.py table under the arms
# its pairs need, with --cluster-trace, for summarize.py (issue #246).
#
# Each (sample, arm) runs once however many pairs it holds: the baseline arm
# for every sample, plus the first arm each pair flipped in.
#
# Usage:
#   trace_pairs.sh <binary> <err.json> <fastq-dir> <pairs.tsv> <out-dir>
#
# Environment:
#   BASE     baseline arm (default insertion), as given to run_member_order.sh
#   EXTRA    the floor run's extra dada flags, e.g. "--band 32 --kmer-size 7".
#            Must match, or the trace is of a different run; summarize.py
#            checks that each pair is reproduced.
#   SUFFIX   FASTQ suffix after the sample name (default .fastq.gz)
#   THREADS  default 4
#
# Writes <out-dir>/<sample>/<arm>.{json,trace.json,log} (":" -> "_"). Then:
#   python3 dev/rename_pairs/summarize.py <pairs.tsv> <out-dir>
set -euo pipefail

BIN="${1:?usage: trace_pairs.sh <binary> <err.json> <fastq-dir> <pairs.tsv> <out-dir>}"
ERR="${2:?missing err.json}"; FQDIR="${3:?missing fastq-dir}"
PAIRS="${4:?missing pairs.tsv}"; OUT="${5:?missing out-dir}"
BASE="${BASE:-insertion}"
SUFFIX="${SUFFIX:-.fastq.gz}"
THREADS="${THREADS:-4}"
read -r -a extra <<< "${EXTRA:-}"

# (sample, arm) jobs: column 1 is the sample, column 8 the arms (first used).
# Arm directory names map back to member orders: shuffle_2 -> shuffle:2.
jobs=$(awk -F'\t' -v base="$BASE" 'NR > 1 {
  split($8, a, ","); arm = a[1]; sub(/_/, ":", arm)
  print $1 "\t" base; print $1 "\t" arm
}' "$PAIRS" | sort -u)
[ -n "$jobs" ] || { echo "no pairs in $PAIRS" >&2; exit 1; }

while IFS=$'\t' read -r smp arm; do
  fq="$FQDIR/$smp$SUFFIX"
  [ -f "$fq" ] || { echo "ERROR: no FASTQ for $smp: $fq" >&2; exit 1; }
  o="$OUT/$smp"; tag=${arm/:/_}
  mkdir -p "$o"
  echo "==> $smp $arm"
  if ! DADA2RS_MEMBER_ORDER="$arm" "$BIN" dada "$fq" --error-model "$ERR" \
      --threads "$THREADS" --verbose --cluster-trace "$o/$tag.trace.json" \
      -o "$o/$tag.json" ${extra[@]+"${extra[@]}"} 2> "$o/$tag.log"; then
    echo "ERROR: dada failed for $smp $arm" >&2; tail -n 20 "$o/$tag.log" >&2; exit 1
  fi
  # Same guard as run_member_order.sh: an arm that silently did not take
  # effect would trace the baseline twice.
  if [ "$arm" != insertion ] && ! grep -q "member order=$arm (OVERRIDDEN)" "$o/$tag.log"; then
    echo "ERROR: $o/$tag.log does not confirm member order=$arm" >&2; exit 1
  fi
done <<< "$jobs"
echo "done: $OUT"
