#!/usr/bin/env bash
# run_screen_sweep_pacbio.sh — sweep the minimizer screen against the k-mer screen
# on PacBio HiFi, on accuracy AND alignment work AND wall time.
#
# Takes ALREADY FILTERED/TRIMMED FASTQs (no primer removal, no filtering) so a
# prepped dataset can be used directly. Single-end throughout; there is no
# merge-pairs step on this platform.
#
# THE ONE DESIGN DECISION THAT MATTERS: every arm shares ONE error model,
# learned once with the k-mer screen.
#
# The screen is active inside learn-errors too -- build_trans_mat aligns each raw
# against its centre THROUGH the screen, so a shrouded pair contributes nothing to
# the transition counts. Letting each arm learn its own model therefore varies the
# error model AND the denoising at once, and the alignment counts that come out
# are not comparable. That confound produced a "13.3% fewer alignments" figure on
# this branch that was wrong by ~23 points once the model was held fixed.
# Holding the model fixed is what makes the sweep measure the screen.
#
# WHY THE CUTOFF IS SWEPT AND NOT JUST CALIBRATED. Calibrating for 100% recall of
# near-neighbour pairs sounds like the safe choice and is the expensive one: on
# PacBio it picked 0.50, which aligns ~10% MORE pairs than the k-mer screen and
# runs 13.5% SLOWER, while 0.42-0.45 aligns 7-11% FEWER, runs ~10% faster, and
# produces a bit-identical ASV set. The last 1.5% of recall bought nothing and
# cost ~20% of the alignment work. Calibrate to bracket the region; let cost
# choose within it.
#
# Usage:
#   run_screen_sweep_pacbio.sh <binary> <filtered-fastq-dir> <out-dir> [threads] [reps]
set -euo pipefail

BIN="${1:?usage: run_screen_sweep_pacbio.sh <binary> <filtered-fastq-dir> <out-dir> [threads] [reps]}"
DATA="${2:?missing filtered-fastq-dir}"
OUT="${3:?missing out-dir}"
THREADS="${4:-16}"
REPS="${5:-3}"

HERE="$(cd "$(dirname "$0")" && pwd)"

# PacBio full-length 16S defaults (match dev/concordance/run_pacbio.sh).
BAND="${BAND:-32}"
KMER="${KMER:-7}"          # k-mer screen size; 7 is the dada2-rs PacBio recommendation
ERRFUN="${ERRFUN:-pacbio}"
NBASES="${NBASES:-1000000000}"

# Sweep grid. Cutoffs bracket the k-mer screen's own pass rate rather than
# chasing 100% recall -- see the header.
# k=8 only: it beat k=9 on ASV set agreement at every cutoff on the 362-sample
# MiSeq run, and k=9 doubles the arm count for no benefit measured anywhere.
KS="${KS:-8}"
# Extended past the matched-pass point (~0.50 on HiFi) so the read-retention
# ZERO CROSSING is inside the grid rather than interpolated off its edge. On
# Illumina the crossing sat ~0.02 ABOVE matched-pass (0.62-0.65 matched, crossing
# 0.636), and it is the smoothest calibration signal available -- monotone, with a
# true zero -- where ASV churn is discrete and count L1 is a flat-bottomed U.
CUTS="${CUTS:-0.40 0.42 0.45 0.48 0.50 0.52 0.55 0.60}"

# Cutoffs to TIME. The accuracy sweep wants the whole grid; timing does not, and
# timing all of it is where this script spends most of its wall clock:
# REPS x (2 + |KS|x|CUTS|) full denoising passes, and the loose-cutoff arms are
# the expensive ones -- measured at 312% of baseline on soil 16S and 818% of
# baseline ALIGNMENTS on pooled ITS2. Timing a cutoff nobody would deploy costs
# more than timing the one they would.
#
# Default: the arms actually worth a wall-clock number -- near the k-mer screen's
# own pass rate, where alignment work is matched and the screen is the only
# variable. Set TIME_CUTS="$CUTS" to time everything, or "" to skip timing.
TIME_CUTS="${TIME_CUTS:-0.48 0.52}"

# COST WARNING. The timing and phase-split sections each run one full denoising
# pass per (arm x rep). On a large pooled run a single pass is enormous, so use
# REPS=1 and one or two TIME_CUTS there; the accuracy grid is cheap by comparison
# (one pass per arm, and arms are cached across re-runs).

# Bound kdist-calibrate: it aligns every sampled pair UNBANDED (deliberately --
# a band would truncate the divergence of distant pairs, which is the quantity
# being calibrated), and emits one CSV row per pair. On 1.5 kb reads that is the
# expensive part of this script.
CAL_PAIRS="${CAL_PAIRS:-300000}"
CAL_UNIQUES="${CAL_UNIQUES:-3000}"

mkdir -p "$OUT"/{derep,models,arms,.timing,.verbose}

shopt -s nullglob
fq=("$DATA"/*.fastq.gz "$DATA"/*.fq.gz "$DATA"/*.fastq "$DATA"/*.fq)
[ ${#fq[@]} -gt 0 ] || { echo "no FASTQ in $DATA" >&2; exit 1; }
echo "==> ${#fq[@]} pre-filtered samples"

echo "==> [1/5] derep"
for f in "${fq[@]}"; do
  b=$(basename "$f"); b=${b%%.*}
  [ -f "$OUT/derep/$b.json" ] || "$BIN" derep "$f" -o "$OUT/derep/$b.json" > /dev/null
done

echo "==> [2/5] calibrate both screens (same pairs, same seed)"
if [ ! -f "$OUT/models/cal_kmer.csv" ]; then
  "$BIN" kdist-calibrate "$OUT"/derep/*.json --k "$KMER" --per-sample \
      --max-uniques "$CAL_UNIQUES" --max-pairs "$CAL_PAIRS" \
      --threads "$THREADS" -o "$OUT/models/cal_kmer.csv" > /dev/null
fi
for K in $KS; do
  [ -f "$OUT/models/cal_mini_k$K.csv" ] && continue
  "$BIN" kdist-calibrate "$OUT"/derep/*.json --k "$KMER" --per-sample \
      --max-uniques "$CAL_UNIQUES" --max-pairs "$CAL_PAIRS" \
      --screen-backend minimizer --minimizer-k "$K" \
      --threads "$THREADS" -o "$OUT/models/cal_mini_k$K.csv" > /dev/null
done
cal_args=("kmer=$OUT/models/cal_kmer.csv")
for K in $KS; do cal_args+=("mini_k$K=$OUT/models/cal_mini_k$K.csv"); done
python3 "$HERE/analyze_kdist_curves.py" "${cal_args[@]}" | tee "$OUT/models/calibration.txt"

echo
echo "==> [3/5] learn-errors ONCE (k-mer screen) -- shared by every arm"
[ -f "$OUT/models/err.json" ] || \
  "$BIN" learn-errors "${fq[@]}" --nbases "$NBASES" --errfun "$ERRFUN" \
      --band "$BAND" --kmer-size "$KMER" --threads "$THREADS" \
      -o "$OUT/models/err.json" > /dev/null

echo "==> [4/5] denoise: baseline, control, and the grid (all on the shared model)"
run_arm() {  # name, extra args...
  local name="$1"; shift
  local d="$OUT/arms/$name"
  [ -d "$d" ] && { echo "    $name (cached)"; return; }
  mkdir -p "$d"
  echo "    $name"
  "$BIN" dada "${fq[@]}" --error-model "$OUT/models/err.json" \
      --band "$BAND" --kmer-size "$KMER" --output-dir "$d" \
      --threads "$THREADS" "$@" > /dev/null
}
run_arm kmer
run_arm kmerctl                       # control channel: identical config, run twice
for K in $KS; do
  for C in $CUTS; do
    run_arm "mini_k${K}_c${C}" --screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C"
  done
done

echo
echo "==> [5/5] accuracy + work"
python3 - "$OUT/arms" <<'PY'
import glob, json, os, sys
root = sys.argv[1]
def load(d):
    asv, al = {}, 0
    stats = []
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        j = json.load(open(p)); s = j["stats"]
        stats.append((s["nalign"], s["nshroud"]))
        for a in j["asvs"]:
            asv[(os.path.basename(p), a["sequence"])] = a["abundance"]
    # `dada-pooled` runs ONCE and writes the SAME global stats into every
    # per-sample output, so summing them multiplies the true count by the sample
    # number. Identical stats across every file is the tell (a per-sample run
    # essentially never produces that), and then one value is the whole run.
    if stats and len(set(stats)) == 1 and len(stats) > 1:
        al = stats[0][0] - stats[0][1]
    else:
        al = sum(a - b for a, b in stats)
    return asv, al
base, bal = load(os.path.join(root, "kmer"))
print(f"{'arm':>18s} {'ASVs':>7s} {'only_k':>7s} {'only_m':>7s} {'abund L1':>10s} {'aligned':>13s} {'vs kmer':>9s} {'reads vs base':>13s}")
print(f"{'kmer (baseline)':>18s} {len(base):7d} {'-':>7s} {'-':>7s} {'-':>10s} {bal:13,d} {'100.0%':>9s}")
for d in sorted(glob.glob(os.path.join(root, "*"))):
    n = os.path.basename(d)
    if n == "kmer": continue
    a, al = load(d)
    sh = set(base) & set(a)
    l1 = sum(abs(base[k] - a[k]) for k in sh)
    tot = sum(base.values())
    # Read retention vs baseline: monotone in cutoff, crosses zero where the
    # minimizer recovers as many reads as the k-mer screen. The smoothest
    # calibration signal available -- churn is discrete and L1 is a flat U.
    reads = sum(a.values()) - sum(base.values())
    print(f"{n:>18s} {len(a):7d} {len(set(base)-set(a)):7d} {len(set(a)-set(base)):7d} "
          f"{100*l1/max(tot,1):9.4f}% {al:13,d} {100*al/bal:8.1f}% {reads:+12,d}")
print("\nkmerctl MUST be identical to the baseline (0 / 0 / 0.0000%).")
print("If it is not, the run is nondeterministic and nothing below is a result.")
PY

echo
echo "==> wall time (arms interleaved across $REPS reps)"
: > "$OUT/timings.tsv"
declare -a T=("kmer::" "kmerctl::")
for K in $KS; do for C in $TIME_CUTS; do T+=("mini_k${K}_c${C}:$K:$C"); done; done
echo "    ${#T[@]} arms x $REPS reps = $(( ${#T[@]} * REPS )) denoising passes"
echo "    (narrow with TIME_CUTS=; the accuracy grid above is unaffected)"
for rep in $(seq 1 "$REPS"); do
  for spec in "${T[@]}"; do
    name="${spec%%:*}"; rest="${spec#*:}"; K="${rest%%:*}"; C="${rest#*:}"
    extra=(); [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C")
    t0=$(python3 -c 'import time;print(time.time())')
    "$BIN" dada "${fq[@]}" --error-model "$OUT/models/err.json" --band "$BAND" \
        --kmer-size "$KMER" --threads "$THREADS" ${extra[@]+"${extra[@]}"} \
        --output-dir "$OUT/.timing" > /dev/null 2>&1
    t1=$(python3 -c 'import time;print(time.time())')
    printf "%s\t%s\t%s\n" "$name" "$rep" "$(python3 -c "print(f'{$t1-$t0:.2f}')")" >> "$OUT/timings.tsv"
  done
done

echo
echo "==> phase split + resource use (one --verbose pass per arm; NOT timed reps)"
echo "Captures the FULL verbose block, not just the screen/align lines: parallel"
echo "efficiency, resident footprint, and the later phases (shuffle, p_update)."
echo "That matters because the screen's share is not fixed -- it ranges 0.9% to"
echo "76.5% across datasets -- and thread-scaling-and-placement.md found the"
echo "optimal thread count is PREDICTED by the screen/align split. A backend that"
echo "changes the split changes the right thread count with it, and shrinking"
echo "b_compare raises the serial phases' share (Amdahl), which caps the speedup"
echo "and idles cores. Peak RSS is recorded per arm for the same reason: the"
echo "screen structures are the largest resident allocation."
# Peak RSS via /usr/bin/time is OPTIONAL. It is absent or flag-incompatible on
# many cluster nodes, and making it a hard dependency of this block meant a
# missing binary silently killed the entire phase-split capture -- which is what
# happened on the first pooled ITS2 run. dada's own `resident Raw footprint` line
# reports the screen structures' size regardless, which is the number that matters
# here.
TIMER=()
if command -v /usr/bin/time > /dev/null 2>&1; then
  if /usr/bin/time -l true > /dev/null 2>&1; then TIMER=(/usr/bin/time -l)
  elif /usr/bin/time -v true > /dev/null 2>&1; then TIMER=(/usr/bin/time -v)
  fi
fi
[ ${#TIMER[@]} -eq 0 ] && echo "    (note: /usr/bin/time unavailable; peak RSS omitted, verbose block still captured)"
{
  for spec in "${T[@]}"; do
    name="${spec%%:*}"; rest="${spec#*:}"; K="${rest%%:*}"; C="${rest#*:}"
    [ "$name" = "kmerctl" ] && continue
    extra=(); [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C")
    echo "===== $name"
    ${TIMER[@]+"${TIMER[@]}"} "$BIN" dada "${fq[@]}" \
        --error-model "$OUT/models/err.json" --band "$BAND" --kmer-size "$KMER" \
        --threads "$THREADS" --verbose \
        ${extra[@]+"${extra[@]}"} --output-dir "$OUT/.verbose" 2>&1 \
      | grep -E "^\[dada\]|maximum resident|Maximum resident|elapsed|real" \
      || echo "    (no output)"
  done
} | tee "$OUT/phase_split.txt"

python3 - "$OUT/timings.tsv" <<'PY'
import statistics, sys
from collections import defaultdict
d = defaultdict(list)
for line in open(sys.argv[1]):
    n, r, t = line.split(); d[n].append(float(t))
base = statistics.median(d["kmer"])
print(f"\n{'arm':>18s} {'median s':>10s} {'min':>8s} {'max':>8s} {'spread':>8s} {'vs kmer':>9s}")
for n in sorted(d):
    ts = d[n]; med = statistics.median(ts)
    print(f"{n:>18s} {med:10.2f} {min(ts):8.2f} {max(ts):8.2f} "
          f"{(max(ts)-min(ts))/med*100:7.1f}% {100*med/base:8.1f}%")
if "kmerctl" in d:
    noise = abs(statistics.median(d["kmerctl"]) - base) / base * 100
    print(f"\nCONTROL CHANNEL (k-mer vs k-mer): {noise:.1f}%.")
    print("Any arm within that of the baseline is inside the noise floor and is")
    print("NOT a speed result. Report the control alongside every number.")
PY
