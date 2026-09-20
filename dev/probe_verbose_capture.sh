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
# truncation, every line started with `[dada]` so the filter dropped nothing,
# both of tee's outputs were byte-identical, and line length does not predict
# which messages are lost. What remains is whether the `2>&1 | grep | tee`
# pipeline is responsible, which needs one controlled run where it happens.
#
# Runs the SAME denoise twice in one job, on one node, over the same inputs,
# differing only in how stderr is captured:
#
#   A  direct    2> file                       (what the sweep does now)
#   B  piped     2>&1 | grep -E ... > file     (what the sweep used to do)
#
# Then counts messages that are emitted once per run_dada call and compares
# both arms against the --metrics-json run count, which is ground truth for
# how many run_dada calls happened.
#
# RUN THIS IN THE MODE YOU ACTUALLY USE. `dada-pooled` denoises the merged
# table ONCE, so every message below should appear once, not once per sample;
# `dada` emits one set per sample. The minimizer work was done under full
# pooling, so `dada-pooled` is the mode that matters for those findings.
#
# Usage:
#   dev/probe_verbose_capture.sh <binary> <err.json> <out-dir> \
#       [subcommand] [threads] [sample-jobs] -- <input.fastq.gz>...
#
#   subcommand    dada (default) | dada-pooled | dada-pseudo
#   sample-jobs   ignored for dada-pooled, which has no such flag
#
# Examples:
#   # pooled, the shape the minimizer findings were produced in
#   dev/probe_verbose_capture.sh ./target/release-native/dada2-rs \
#       sweep/base/errF.json /tmp/probe dada-pooled 48 -- filtered/*_F.fastq.gz
#
#   # per-sample, the shape the ITS2 sweep ran
#   dev/probe_verbose_capture.sh ./target/release-native/dada2-rs \
#       sweep/base/errF.json /tmp/probe dada 48 12 -- filtered/*_F.fastq.gz
set -euo pipefail

BIN="${1:?usage: probe_verbose_capture.sh <binary> <err.json> <out-dir> [subcmd] [threads] [jobs] -- <inputs>...}"
ERR="${2:?missing error model}"
OUT="${3:?missing out-dir}"
CMD="${4:-dada}"
THREADS="${5:-48}"
JOBS="${6:-12}"
shift 6 2> /dev/null || shift $#
[ "${1:-}" = "--" ] && shift
[ $# -gt 0 ] || { echo "no input files given (did you forget the -- separator?)" >&2; exit 2; }

case "$CMD" in
  dada | dada-pooled | dada-pseudo) ;;
  *) echo "subcommand must be dada, dada-pooled or dada-pseudo (got '$CMD')" >&2; exit 2 ;;
esac

# `dada-pooled` denoises once over the merged table and has no --sample-jobs.
JOBFLAG=(--sample-jobs "$JOBS")
[ "$CMD" = "dada-pooled" ] && JOBFLAG=()

FILTER='^\[(dada|derep)|maximum resident|Maximum resident|elapsed|real'

mkdir -p "$OUT"
echo "==> $CMD over $# input(s), --threads $THREADS ${JOBFLAG[*]:-(no --sample-jobs)}"
echo

echo "==> arm A: direct redirect"
"$BIN" "$CMD" "$@" \
    --error-model "$ERR" --threads "$THREADS" ${JOBFLAG[@]+"${JOBFLAG[@]}"} --verbose \
    --metrics-json "$OUT/a.metrics.json" --metrics-attribution \
    --output-dir "$OUT/a.out" \
    > "$OUT/a.stdout" 2> "$OUT/a.stderr"

echo "==> arm B: 2>&1 | grep"
"$BIN" "$CMD" "$@" \
    --error-model "$ERR" --threads "$THREADS" ${JOBFLAG[@]+"${JOBFLAG[@]}"} --verbose \
    --metrics-json "$OUT/b.metrics.json" --metrics-attribution \
    --output-dir "$OUT/b.out" 2>&1 \
  | grep -E "$FILTER" > "$OUT/b.stderr"

runs=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))['runs']))" \
        "$OUT/a.metrics.json")

# dada-pseudo calls run_dada TWICE per sample (round 1 without priors, round 2
# with), but the metrics document records round 2 only -- round 1 goes through
# a path that does not serialize per-sample output. So the prose carries twice
# as many blocks as the JSON has entries, and that is expected, not loss.
expected=$runs
if [ "$CMD" = "dada-pseudo" ]; then
  expected=$((runs * 2))
  echo
  echo "NOTE dada-pseudo runs two rounds per sample and the JSON records round 2"
  echo "     only, so the prose is expected to carry 2x the JSON entry count."
fi

echo
echo "run_dada calls, per --metrics-json: $runs (expecting $expected prose block(s))"
echo "Each message below is emitted once per call, so each should equal that."
echo

# The ASV summary and the per-sample `wrote` line come from
# denoise_and_serialize, which `dada-pooled` does not use -- it has its own
# serialization. Expecting them under pooling would report a false loss.
MARKERS=("resident Raw footprint" "kmer8 fill" "kmer8 pooled diversity"
         "compare attribution (of" "] phase times" "compare split (of"
         "map parallel efficiency")
if [ "$CMD" = "dada-pooled" ]; then
  MARKERS+=("phase wall times" "derep split (of")
else
  MARKERS+=("ASV(s) from" "] wrote")
fi

printf "%-26s %10s %10s\n" "message" "A direct" "B piped"
printf "%-26s %10s %10s\n" "-------" "--------" "-------"
for pat in "${MARKERS[@]}"; do
  a=$(grep -cF "$pat" "$OUT/a.stderr" || true)
  b=$(grep -cF "$pat" "$OUT/b.stderr" || true)
  # `] wrote` also matches the single "wrote run metrics to ..." line, which is
  # once per RUN, not once per sample. Discount it so the row is comparable.
  if [ "$pat" = "] wrote" ]; then
    a=$((a - $(grep -cF "wrote run metrics" "$OUT/a.stderr" || true)))
    b=$((b - $(grep -cF "wrote run metrics" "$OUT/b.stderr" || true)))
  fi
  mark=""
  # The ASV summary and `wrote` are emitted once per SAMPLE, not once per
  # run_dada call, so under pseudo they track the JSON count, not the doubled one.
  want=$expected
  case "$pat" in "ASV(s) from" | "] wrote") want=$runs ;; esac
  [ "$a" != "$want" ] && mark="$mark  A-short"
  [ "$b" != "$want" ] && mark="$mark  B-short"
  printf "%-26s %10s %10s%s\n" "$pat" "$a" "$b" "$mark"
done

# Pooled only: the pipeline phase times exist in the JSON as pipeline.derep /
# merge / dada / output, and in the prose as one `[dada-pooled] phase wall
# times` line. Nothing has ever cross-checked them, because the ITS2 sweep ran
# per-sample and those fields do not exist there.
if [ "$CMD" = "dada-pooled" ]; then
  echo
  echo "PIPELINE PHASE TIMES (pooled only; never previously cross-checked)"
  python3 - "$OUT/a.stderr" "$OUT/a.metrics.json" <<'PY'
import json, re, sys
prose = open(sys.argv[1], errors="replace").read()
doc = json.load(open(sys.argv[2]))
m = re.search(
    r"phase wall times: derep=([\d.]+)s.*?merge=([\d.]+)s.*?"
    r"run_dada=([\d.]+)s.*?output=([\d.]+)s", prose)
if not m:
    print("    no `phase wall times` line found in arm A -- was --verbose on?")
    sys.exit(0)
pipe = doc.get("pipeline", {})
print(f"    {'field':10s} {'prose':>10s} {'json':>10s} {'diff':>9s}   (prose is 1dp)")
bad = 0
for i, key in enumerate(("derep", "merge", "dada", "output")):
    p = float(m.group(i + 1))
    j = pipe.get(key)
    if j is None:
        print(f"    {key:10s} {p:10.2f} {'absent':>10s}   <-- MISSING FROM JSON")
        bad += 1
        continue
    d = abs(p - j)
    flag = "" if d <= 0.05 else "   <-- MISMATCH"
    if flag:
        bad += 1
    print(f"    {key:10s} {p:10.2f} {j:10.2f} {d:9.2f}{flag}")
print("    (prose prints 1 decimal, so 0.05s tolerance)")
if bad:
    print("    *** prose and JSON disagree on the pipeline phases.")
PY
fi

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
