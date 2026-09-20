#!/usr/bin/env bash
# probe_verbose_capture.sh -- is the sweep's prose capture losing lines?
#
# The ITS2 sweep's phase_split.txt is missing verbose lines that the code
# cannot explain: `kmer8 fill` appeared 2 times and `kmer8 pooled diversity`
# 29 times across the whole file, although they are two ADJACENT eprintln!
# calls inside one guard (src/dada.rs:595 and :601) and therefore execute the
# same number of times. Loss rates ranged from 0.7% (`phase times`) to 95%
# (the ASV summary).
#
# Everything testable has already been ruled out: no spliced lines, no
# truncation, every line starts with `[dada]` so the grep filtered nothing,
# both of tee's outputs are byte-identical, and line length does not predict
# which messages are lost. The remaining question is simply whether the
# `2>&1 | grep | tee` pipeline is responsible, and that needs one controlled
# run in the environment where it happens.
#
# This runs the SAME denoise twice in one job, on one node, over the same
# inputs, differing only in how stderr is captured:
#
#   A  direct    2> file                       (what the sweep does now)
#   B  piped     2>&1 | grep -E ... > file      (what the sweep used to do)
#
# Then it counts, in each, messages that must appear exactly once per sample.
# If A is complete and B is short, the pipeline was the cause. If both are
# short, the loss is inside the process and arm A holds the raw unfiltered
# stream to look at. The `--metrics-json` run count is printed alongside as
# ground truth for how many samples actually ran.
#
# Usage:
#   dev/probe_verbose_capture.sh <binary> <error-model.json> <out-dir> \
#       [threads] [sample-jobs] -- <input.fastq.gz>...
#
# Example (matching the ITS2 sweep's shape):
#   dev/probe_verbose_capture.sh ./target/release-native/dada2-rs \
#       sweep/base/errF.json /tmp/capture-probe 48 12 -- filtered/*_F.fastq.gz
set -euo pipefail

BIN="${1:?usage: probe_verbose_capture.sh <binary> <err.json> <out-dir> [threads] [jobs] -- <inputs>...}"
ERR="${2:?missing error model}"
OUT="${3:?missing out-dir}"
THREADS="${4:-48}"
JOBS="${5:-12}"
shift 5 || shift $#
[ "${1:-}" = "--" ] && shift
[ $# -gt 0 ] || { echo "no input files given (did you forget the -- separator?)" >&2; exit 2; }

mkdir -p "$OUT"
echo "==> $# input(s), --threads $THREADS --sample-jobs $JOBS"
echo

# Arm A: straight to a file. No pipe, no grep, no tee.
echo "==> arm A: direct redirect"
"$BIN" dada "$@" \
    --error-model "$ERR" --threads "$THREADS" --sample-jobs "$JOBS" --verbose \
    --metrics-json "$OUT/a.metrics.json" --metrics-attribution \
    --output-dir "$OUT/a.out" \
    > "$OUT/a.stdout" 2> "$OUT/a.stderr"

# Arm B: the old pipeline, verbatim.
echo "==> arm B: 2>&1 | grep"
"$BIN" dada "$@" \
    --error-model "$ERR" --threads "$THREADS" --sample-jobs "$JOBS" --verbose \
    --metrics-json "$OUT/b.metrics.json" --metrics-attribution \
    --output-dir "$OUT/b.out" 2>&1 \
  | grep -E "^\[(dada|derep)|maximum resident|Maximum resident|elapsed|real" \
  > "$OUT/b.stderr"

runs=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))['runs']))" \
        "$OUT/a.metrics.json")

echo
echo "Samples actually denoised, per --metrics-json: $runs"
echo "Every message below is emitted once per sample, so each should equal that."
echo
printf "%-26s %10s %10s\n" "message" "A direct" "B piped"
printf "%-26s %10s %10s\n" "-------" "--------" "-------"
for pat in "resident Raw footprint" "kmer8 fill" "kmer8 pooled diversity" \
           "compare attribution (of" "] phase times" "compare split (of" \
           "map parallel efficiency" "ASV(s) from" "] wrote"; do
  a=$(grep -cF "$pat" "$OUT/a.stderr" || true)
  b=$(grep -cF "$pat" "$OUT/b.stderr" || true)
  # `] wrote` also matches the single "wrote run metrics to ..." line, which is
  # once per RUN, not once per sample. Discount it so the row is comparable.
  if [ "$pat" = "] wrote" ]; then
    a=$(( a - $(grep -cF "wrote run metrics" "$OUT/a.stderr" || true) ))
    b=$(( b - $(grep -cF "wrote run metrics" "$OUT/b.stderr" || true) ))
  fi
  mark=""
  [ "$a" != "$runs" ] && mark="$mark  A-short"
  [ "$b" != "$runs" ] && mark="$mark  B-short"
  printf "%-26s %10s %10s%s\n" "$pat" "$a" "$b" "$mark"
done

echo
echo "  kmer8 fill and kmer8 pooled diversity are adjacent eprintln! calls in one"
echo "  guard. If they differ in either arm, lines are being dropped after the"
echo "  process emitted them -- that cannot come from execution."
echo
echo "  A complete, B short  -> the pipeline was the cause; the sweep fix stands."
echo "  both short           -> the loss is inside the process. $OUT/a.stderr is"
echo "                          the raw unfiltered stream; send that."
echo
echo "raw streams kept: $OUT/a.stderr  $OUT/b.stderr"
