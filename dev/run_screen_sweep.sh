#!/usr/bin/env bash
# run_screen_sweep.sh — sweep the minimizer screen's (k, cutoff) grid against the
# k-mer screen, on three axes at once: ASV agreement, alignment WORK, and wall time.
#
# Why all three. The minimizer screen has a recall/cost knob, and the two ends
# pull against each other:
#
#   - too tight  -> shrouds genuine neighbours, raws never reach their parent
#                   cluster and birth spurious low-abundance ASVs (fragmentation)
#   - too loose  -> retains everything, but aligns MORE pairs than the k-mer
#                   screen, and b_compare is align-dominated (77-88%), so the
#                   run gets slower rather than faster
#
# On the 20-sample MiSeq SOP, cutoff 0.42 sits at the first failure (retains only
# 60.8% of pairs below 10% divergence) and 0.72 at the second (31.7% MORE
# alignments than k-mer). The useful region is between them, and it has to be
# found empirically per platform because the distance distribution moves with
# read length. See docs/findings/minimizer-screening.md.
#
# METHODOLOGY NOTES, learned the hard way on this branch:
#
#   * The control arm is not optional. The k-mer backend run TWICE establishes
#     the noise floor for both ASV agreement (it is exactly zero) and wall time
#     (it is not). A difference smaller than the control channel is not a result.
#   * Timing arms are INTERLEAVED, not run back-to-back, so thermal drift and
#     page-cache warming do not align with the arm boundary.
#   * ASV agreement is measured on SET IDENTITY and the full sample x ASV COUNT
#     MATRIX, never on n_asv, which stays flat while the sets churn underneath.
#
# Usage:
#   run_screen_sweep.sh <dada2rs-binary> <sop-style-dir> <out-dir> [threads] [reps]
#
#   <sop-style-dir> holds paired <sample>F.fastq.gz / <sample>R.fastq.gz (the
#   layout dev/concordance/run_illumina.sh expects). For the MiSeq SOP, symlink
#   the _R1_001/_R2_001 names across first.
#
#   PREFILTERED=1 when the inputs are ALREADY trimmed and filtered. Required for
#   ITS2, whose primers sit behind heterogeneity spacers and must come off with
#   cutadapt first (docs/findings/reading-the-prep.md). Without it this script
#   re-runs filter-and-trim with --trunc-len, and FIXED-LENGTH truncation on a
#   length-variable amplicon destroys the length variation being studied -- which
#   for ITS2 is the whole point of the test.
set -euo pipefail

BIN="${1:?usage: run_screen_sweep.sh <binary> <data-dir> <out-dir> [threads] [reps]}"
DATA="${2:?missing data-dir}"
OUT="${3:?missing out-dir}"
THREADS="${4:-8}"
REPS="${5:-3}"

HERE="$(cd "$(dirname "$0")" && pwd)"
RUN="$HERE/concordance/run_illumina.sh"

# Grid. Cutoffs span the fragmentation end to the over-alignment end; k is the
# sketch size, where 8 zeroed out every missed pair below 10 substitutions on the
# SOP and 11 (the shipped default) did not.
KS="${KS:-8}"
# Deliberately spans BELOW the k-mer screen's own 0.42 and ABOVE the zero-churn
# point, so the failure modes at both ends are visible rather than inferred: too
# tight shrouds genuine neighbours and fragments clusters, too loose aligns far
# more pairs for nothing. A grid that only covers the good region cannot show
# where the edges are.
CUTS="${CUTS:-0.40 0.45 0.50 0.55 0.60 0.62 0.63 0.64 0.65 0.70 0.75 0.80}"

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
TIME_CUTS="${TIME_CUTS:-0.62 0.64}"

# COST WARNING. The timing and phase-split sections each run one full denoising
# pass per (arm x rep), and on a large pooled run a single pass is enormous --
# pooled soil 16S is 11.8 BILLION comparisons. A 3-rep x 12-arm timing section is
# 36 such passes and will not finish. For pooled runs use REPS=1 and one or two
# TIME_CUTS; the accuracy grid above is cheap by comparison because it is one pass
# per arm and its arms are cached across re-runs.

# Reuse an existing error model for EVERY arm instead of learning one per arm.
#
# The screen is active inside learn-errors (build_trans_mat aligns each raw
# against its centre THROUGH the screen), so per-arm models make the arms differ
# in two ways at once and their ALIGNMENT COUNTS stop being comparable. The
# 362-sample MiSeq sweep ran without this and its work columns carry that caveat.
# Point ERR_DIR at a completed k-mer run to hold the model fixed:
#     ERR_DIR=/path/to/base bash dev/run_screen_sweep.sh ...
# ASV comparisons are valid either way -- those compare complete configurations.
ERR_DIR="${ERR_DIR:-}"

# Error-fitting function, passed through to run_illumina.sh and held FIXED across
# every arm. `loess` is the default; NovaSeq / binned-quality runs (soil 16S,
# ITS2) may want `binned-qual`. Holding it fixed is what keeps this a screen
# comparison rather than a screen-and-errfun comparison.
ERRFUN="${ERRFUN:-loess}"
# Companion flags for errfuns that need them, e.g.
#   ERRFUN=binned-qual ERRFUN_ARGS='--binned-quals 2,12,23,37'
ERRFUN_ARGS="${ERRFUN_ARGS:-}"

# Denoising mode, passed to run_illumina.sh: unset/false = per-sample `dada`,
# pseudo = `dada-pseudo`, true = `dada-pooled`.
#
# POOL=true is the configuration in which a SCREEN comparison is most
# informative. The screen is memory-bound, so its per-comparison cost scales with
# the working set rather than the vector size: the same k=5 vector costs 44
# ns/comp on a 1,979-unique sample (0.9% of busy time) and 841 ns/comp on a
# 272,574-unique pooled run (16.7%). A per-sample sweep will therefore show a
# sub-1% screen share and be dominated entirely by alignment count.
POOL="${POOL:-false}"

# Decouple the learn-errors screen from the denoising screen, per
# docs/findings/kdist-cutoff-decoupling.md: for the k-mer screen, keeping
# learn-errors lenient while tightening dada was both safe and fast, while
# tightening both churned real ASVs. Set LEARN_CUTOFF to a lenient value (e.g.
# 0.70) and the sweep's cutoffs then apply to DENOISING only.
LEARN_CUTOFF="${LEARN_CUTOFF:-}"

mkdir -p "$OUT"/.verbose "$OUT"

aligned_count() {  # sum (nalign - nshroud) over a run's forward dada outputs
  python3 - "$1" <<'PY'
import glob, json, os, sys
stats = []
for p in glob.glob(os.path.join(sys.argv[1], "dada_fwd", "*.json")):
    s = json.load(open(p))["stats"]
    stats.append((s["nalign"], s["nshroud"]))
# `dada-pooled` runs ONCE and writes the SAME global stats into every per-sample
# output, so summing them multiplies the true count by the sample number.
# Identical stats across every file is the tell -- a per-sample run essentially
# never produces that -- and then one value is the whole run.
if stats and len(set(stats)) == 1 and len(stats) > 1:
    na, ns = stats[0]
else:
    na = sum(a for a, _ in stats); ns = sum(b for _, b in stats)
print(f"{na} {ns} {na - ns}")
PY
}

echo "==> baseline: k-mer screen (production default)"
[ -d "$OUT/base" ] || PREFILTERED="${PREFILTERED:-}" ERR_DIR="${ERR_DIR:-}" ERRFUN="$ERRFUN" ERRFUN_ARGS="$ERRFUN_ARGS" POOL="$POOL" LEARN_CUTOFF="$LEARN_CUTOFF" \
  bash "$RUN" "$BIN" "$DATA" "$OUT/base" "$THREADS" > "$OUT/base.log" 2>&1

echo "==> control: k-mer screen AGAIN (establishes the ASV noise floor)"
[ -d "$OUT/control" ] || PREFILTERED="${PREFILTERED:-}" ERR_DIR="${ERR_DIR:-}" ERRFUN="$ERRFUN" ERRFUN_ARGS="$ERRFUN_ARGS" POOL="$POOL" LEARN_CUTOFF="$LEARN_CUTOFF" \
  bash "$RUN" "$BIN" "$DATA" "$OUT/control" "$THREADS" > "$OUT/control.log" 2>&1
echo "    control vs baseline (MUST be identical, or nothing below is interpretable):"
# compare_seqtab_matrix exits 1 when the tables differ, which for the CONTROL is
# a finding to print, not a reason to abort.
{ python3 "$HERE/compare_seqtab_matrix.py" \
    "$OUT/base/seqtab.nochim.json" "$OUT/control/seqtab.nochim.json" \
    --label-a base --label-b control || true; } | tail -2 | sed 's/^/      /'

echo
echo "==> sweeping the (k, cutoff) grid"
for K in $KS; do
  for C in $CUTS; do
    d="$OUT/mini_k${K}_c${C}"
    [ -d "$d" ] && { echo "    k=$K cutoff=$C (cached)"; continue; }
    echo "    k=$K cutoff=$C"
    SCREEN_BACKEND=minimizer MINIMIZER_K="$K" SCREEN_CUTOFF="$C" \
      ERR_DIR="${ERR_DIR:-}" PREFILTERED="${PREFILTERED:-}" ERRFUN="$ERRFUN" ERRFUN_ARGS="$ERRFUN_ARGS" POOL="$POOL" LEARN_CUTOFF="$LEARN_CUTOFF" \
      bash "$RUN" "$BIN" "$DATA" "$d" "$THREADS" > "$d.log" 2>&1
  done
done

echo
echo "==> error-model provenance"
python3 - "$OUT" <<'PY'
import glob, json, os, sys
root = sys.argv[1]
def flat(p):
    j = json.load(open(p))
    return [x for r in j["err_out"] for x in r]
base = os.path.join(root, "base", "errF.json")
if not os.path.exists(base):
    print("    (no errF.json; skipping)")
    raise SystemExit
b = flat(base)
same = diff = 0
for f in sorted(glob.glob(os.path.join(root, "*", "errF.json"))):
    if flat(f) == b:
        same += 1
    else:
        diff += 1
print(f"    arms sharing the baseline error model: {same}; differing: {diff}")
if diff == 0:
    print("    => alignment counts across arms are COMPARABLE.")
else:
    print("    => alignment counts across arms are NOT comparable: each arm")
    print("       learned its own model, and the screen shapes the model.")
    print("       Re-run with ERR_DIR=<base-run> to isolate the screen.")
    print("       ASV columns below are unaffected.")
PY

echo
echo "======================= ACCURACY + WORK ======================="
read -r bna bns bal <<< "$(aligned_count "$OUT/base")"
printf "%-22s %8s %8s %12s %10s %12s %s\n" "arm" "ASVs" "churn" "aligned" "vs kmer" "reads vs base" "count-matrix L1"
printf "%-22s %8s %8s %12s %10s %12s %s\n" "k-mer (baseline)" \
    "$(python3 -c "import json;print(len(json.load(open('$OUT/base/seqtab.nochim.json'))['sequences']))")" \
    "-" "$bal" "100.0%" "-" "-"
for K in $KS; do
  for C in $CUTS; do
    d="$OUT/mini_k${K}_c${C}"
    [ -d "$d" ] || continue
    read -r na ns al <<< "$(aligned_count "$d")"
    n=$(python3 -c "import json;print(len(json.load(open('$d/seqtab.nochim.json'))['sequences']))")
    # Read retention vs baseline: the calibration signal that is monotone in
    # cutoff and crosses zero where the minimizer recovers as many reads as the
    # k-mer screen (~0.64 on every screen-dominated dataset measured).
    rd=$(python3 -c "
import json
a=json.load(open('$OUT/base/seqtab.nochim.json')); b=json.load(open('$d/seqtab.nochim.json'))
ta=sum(sum(r) for r in a['counts']); tb=sum(sum(r) for r in b['counts'])
print(f'{tb-ta:+d}')")
    # `|| true` on both: these pipelines end in grep, and under `set -o pipefail`
    # a grep that matches nothing returns 1, which `set -e` turns into an abort.
    # A zero-churn arm produces exactly that -- so without this the script dies on
    # its own best result.
    churn=$(python3 "$HERE/compare_asvs.py" --baseline k="$OUT/base/seqtab.nochim.json" \
              --compare m="$d/seqtab.nochim.json" 2>/dev/null \
              | grep -o "churn=[0-9]*" | head -1 | cut -d= -f2 || true)
    l1=$(python3 "$HERE/compare_seqtab_matrix.py" "$OUT/base/seqtab.nochim.json" \
              "$d/seqtab.nochim.json" --label-a k --label-b m 2>/dev/null \
           | grep -o "L1 = [0-9]* reads ([0-9.]*%" | grep -o "([0-9.]*%" | tr -d '(' || true)
    printf "%-22s %8s %8s %12s %9s%% %12s %s\n" "mini k=$K cut=$C" "$n" "${churn:-?}" "$al" \
        "$(python3 -c "print(f'{100*$al/$bal:.1f}')")" "${rd:-?}" "${l1:-?}"
  done
done

echo
echo "========================== WALL TIME =========================="
echo "Arms interleaved across $REPS reps; the k-mer arm appears twice so the"
echo "control channel carries the same thermal/cache exposure as the test arms."
echo "Timing the dada step only -- filter/merge/chimera are identical across arms."
case "$POOL" in
  true)   DADA_CMD="dada-pooled" ;;
  pseudo) DADA_CMD="dada-pseudo" ;;
  *)      DADA_CMD="dada" ;;
esac
echo "Denoising mode: $DADA_CMD (POOL=$POOL)"

# Under PREFILTERED the run wrote no filtered/ dir -- it used the inputs
# directly, so the timing arms must too.
if [ -n "${PREFILTERED:-}" ]; then
  filtF=("$DATA"/*F.fastq.gz)
else
  filtF=("$OUT"/base/filtered/*_F_filt.fastq.gz)
fi
if [ ${#filtF[@]} -eq 0 ]; then
  echo "    (no forward reads found for timing; skipping)" >&2
  REPS=0
fi
declare -a ARMS=("kmer::" "kmerctl::")
for K in $KS; do for C in $TIME_CUTS; do ARMS+=("mini_k${K}_c${C}:$K:$C"); done; done

: > "$OUT/timings.tsv"
echo "    ${#ARMS[@]} arms x $REPS reps = $(( ${#ARMS[@]} * REPS )) denoising passes"
echo "    (narrow with TIME_CUTS=; the accuracy grid above is unaffected)"
for rep in $(seq 1 "$REPS"); do
  for arm in "${ARMS[@]}"; do
    name="${arm%%:*}"; rest="${arm#*:}"; K="${rest%%:*}"; C="${rest#*:}"
    extra=()
    [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C")
    t0=$(python3 -c 'import time;print(time.time())')
    # Same denoising mode as the arms, or the timings describe a different
    # workload than the accuracy table above.
    "$BIN" $DADA_CMD "${filtF[@]}" --error-model "$OUT/base/errF.json" \
        --output-dir "$OUT/.timing" --threads "$THREADS" \
        ${extra[@]+"${extra[@]}"} > /dev/null 2>&1
    t1=$(python3 -c 'import time;print(time.time())')
    printf "%s\t%s\t%s\n" "$name" "$rep" \
        "$(python3 -c "print(f'{$t1-$t0:.2f}')")" >> "$OUT/timings.tsv"
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
  for spec in "${ARMS[@]}"; do
    name="${spec%%:*}"; rest="${spec#*:}"; K="${rest%%:*}"; C="${rest#*:}"
    [ "$name" = "kmerctl" ] && continue
    extra=(); [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C")
    echo "===== $name"
    ${TIMER[@]+"${TIMER[@]}"} "$BIN" $DADA_CMD "${filtF[@]}" \
        --error-model "$OUT/base/errF.json" --threads "$THREADS" --verbose \
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
    name, rep, t = line.split()
    d[name].append(float(t))
base = statistics.median(d["kmer"]) if "kmer" in d else None
ctl = d.get("kmerctl")
print(f"\n{'arm':>22s} {'median s':>10s} {'min':>8s} {'max':>8s} {'spread':>8s} {'vs kmer':>9s}")
for name, ts in d.items():
    med = statistics.median(ts)
    spread = (max(ts) - min(ts)) / med * 100
    rel = f"{100 * med / base:.1f}%" if base else "-"
    print(f"{name:>22s} {med:10.2f} {min(ts):8.2f} {max(ts):8.2f} {spread:7.1f}% {rel:>9s}")
if ctl and base:
    noise = abs(statistics.median(ctl) - base) / base * 100
    print(f"\nCONTROL CHANNEL: k-mer vs k-mer differs by {noise:.1f}%.")
    print("Any arm closer to the baseline than this is INSIDE the noise floor")
    print("and must not be reported as a speed difference.")
PY
