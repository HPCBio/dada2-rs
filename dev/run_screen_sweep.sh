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

# Winnowing windows to sweep. w controls SAMPLING DENSITY -- the sketch keeps
# ~2/(w+1) of positions -- and it has never been varied on this branch, so every
# result here is a single w=5 point.
#
# The reason to sweep it on Illumina is not tuning, it is a test of why the
# platforms differ. Sketch size tracks read length: 475 entries/raw on 1.5 kb
# HiFi, but only **74** on 231 bp ITS2. The screen distance is a sum of minima
# over those entries, so its granularity is 1/74 against 1/475 and its relative
# sampling error is ~2.5x larger -- and HiFi is exactly concordant while Illumina
# churns. If churn falls monotonically as w drops (raising the entry count), the
# platform gap is estimator variance from a small sketch, not the error-profile
# story currently in the findings page.
#
# w=1 is the informative endpoint: every k-mer becomes its own window minimum, so
# the sketch degenerates to the full k-mer multiset (bar consecutive duplicates
# inside homopolymers, which emit once by design). That arm isolates WINNOWING
# NOISE from the choice of k -- if w=1 churns zero, all churn is sampling.
WS="${WS:-5}"
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

# Cutoffs to ALSO time with the inverted index DISABLED, as separate `_noidx`
# arms. Timing and phase-split only -- the accuracy grid here is the expensive
# part (full pipeline per arm), and index exactness is already established on
# PacBio, where two `_noidx` arms reproduced their index-on twins cell for cell.
#
# OPT-IN here, unlike the PacBio script, but the question is more open on this
# platform, not less. On per-sample PacBio the index WON: the per-pair merge-join
# cost 2379 ns/comp against the k-mer sweep's 2460 -- both memory-latency-bound,
# so the serial scatter was worth paying. Two things differ on a pooled Illumina
# run and both push the other way. Sketches are ~4x smaller on 300 bp reads than
# on 1.5 kb HiFi, making the merge-join proportionally cheaper; and a single pool
# at 48 threads has no cross-sample concurrency to hide serial time behind, where
# per-sample PacBio ran 12 samples at once. Pooled ITS2 measured setup at 23.7%
# of compare for exactly that reason. So set NOIDX_CUTS on a pooled run.
NOIDX_CUTS="${NOIDX_CUTS:-}"

# Cutoffs to ALSO time with the index left to `decide_index` (the shipped
# behaviour), as `_auto` arms. Once the selection rule exists this is the arm
# that matters: `_c<C>` forces it on, `_noidx` forces it off, only `_auto`
# measures what a user gets. Defaults to TIME_CUTS.
AUTO_CUTS="${AUTO_CUTS:-$TIME_CUTS}"

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

# Arm directories and timing rows are cached by NAME, and the name encodes the
# cutoff but not the denoising mode or the error function. Pointing a POOL=true
# run at a per-sample output directory would therefore reuse every per-sample arm
# as though it were pooled -- silently, and with the accuracy table and timings
# both looking plausible. Stamp the configuration and refuse to mix.
STAMP="$OUT/.sweep_mode"
# POOL, not DADA_CMD: the two scripts define DADA_CMD at different points and
# this guard must not depend on that ordering (it broke the ITS run once).
WANT="pool=$POOL errfun=$ERRFUN"
mkdir -p "$OUT"
if [ -f "$STAMP" ]; then
  HAVE="$(cat "$STAMP")"
  if [ "$HAVE" != "$WANT" ]; then
    echo "ERROR: $OUT was built with a different configuration." >&2
    echo "  existing: $HAVE" >&2
    echo "  requested: $WANT" >&2
    echo "  Cached arms and timings are keyed by name only, so reusing this" >&2
    echo "  directory would mix the two. Use a new out-dir. To keep the" >&2
    echo "  expensive cached inputs, copy them over first:" >&2
    echo "    mkdir -p <new-out> && cp -a $OUT/models $OUT/derep <new-out>/" >&2
    exit 1
  fi
else
  printf '%s' "$WANT" > "$STAMP"
fi

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
# Arm name. The `_w<N>` segment appears ONLY for a non-default w, so every arm
# already on disk from a w=5 run keeps its name and stays cached -- adding the w
# axis must not silently invalidate a completed grid.
arm_name() {  # k, w, cutoff
  if [ "$2" = "5" ]; then echo "mini_k$1_c$3"; else echo "mini_k$1_w$2_c$3"; fi
}

echo "==> sweeping the (k, w, cutoff) grid"
for K in $KS; do
  for W in $WS; do
    for C in $CUTS; do
      d="$OUT/$(arm_name "$K" "$W" "$C")"
      [ -d "$d" ] && { echo "    k=$K w=$W cutoff=$C (cached)"; continue; }
      echo "    k=$K w=$W cutoff=$C"
      SCREEN_BACKEND=minimizer MINIMIZER_K="$K" MINIMIZER_W="$W" SCREEN_CUTOFF="$C" \
        ERR_DIR="${ERR_DIR:-}" PREFILTERED="${PREFILTERED:-}" ERRFUN="$ERRFUN" ERRFUN_ARGS="$ERRFUN_ARGS" POOL="$POOL" LEARN_CUTOFF="$LEARN_CUTOFF" \
        bash "$RUN" "$BIN" "$DATA" "$d" "$THREADS" > "$d.log" 2>&1
    done
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
# spec = name:minimizer-k:cutoff:w:index  (empty k = k-mer arm; index 0 = no index)
declare -a ARMS=("kmer::::" "kmerctl::::")
for K in $KS; do for W in $WS; do
  # Forced ON, not defaulted: `decide_index` now chooses, so an unqualified arm
  # would change meaning between binaries while keeping its name -- and
  # timings.tsv keys on (arm, rep), so old and new replicates would be pooled.
  for C in $TIME_CUTS;  do ARMS+=("$(arm_name "$K" "$W" "$C"):$K:$C:$W:1"); done
  for C in $NOIDX_CUTS; do ARMS+=("$(arm_name "$K" "$W" "$C")_noidx:$K:$C:$W:0"); done
  for C in $AUTO_CUTS;  do ARMS+=("$(arm_name "$K" "$W" "$C")_auto:$K:$C:$W:"); done
done; done

# Append, do not truncate: timing passes are the expensive part (a pooled pass is
# billions of comparisons), so raising REPS on a re-run should ADD replicates
# rather than discard the ones already paid for. Each (arm, rep) is skipped if
# already recorded.
touch "$OUT/timings.tsv"
echo "    ${#ARMS[@]} arms x $REPS reps = $(( ${#ARMS[@]} * REPS )) denoising passes"
echo "    (narrow with TIME_CUTS=; the accuracy grid above is unaffected)"
for rep in $(seq 1 "$REPS"); do
  for arm in "${ARMS[@]}"; do
    IFS=: read -r name K C W IDX <<< "$arm"
    extra=()
    [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --minimizer-w "$W" --kdist-cutoff "$C")
    env=(); [ -n "$IDX" ] && env=(env "DADA2RS_MINIMIZER_INDEX=$IDX")
    # Already timed on an earlier invocation? Skip it, so raising REPS adds
    # replicates instead of redoing the ones already paid for.
    if awk -F'\t' -v n="$name" -v r="$rep" '$1==n && $2==r{f=1} END{exit !f}' \
         "$OUT/timings.tsv" 2>/dev/null; then
      echo "    $name rep $rep (cached)"
      continue
    fi
    t0=$(python3 -c 'import time;print(time.time())')
    # Same denoising mode as the arms, or the timings describe a different
    # workload than the accuracy table above.
    ${env[@]+"${env[@]}"} "$BIN" $DADA_CMD "${filtF[@]}" \
        --error-model "$OUT/base/errF.json" \
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
    IFS=: read -r name K C W IDX <<< "$spec"
    [ "$name" = "kmerctl" ] && continue
    extra=(); [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --minimizer-w "$W" --kdist-cutoff "$C")
    env=(); [ -n "$IDX" ] && env=(env "DADA2RS_MINIMIZER_INDEX=$IDX")
    echo "===== $name"
    ${env[@]+"${env[@]}"} ${TIMER[@]+"${TIMER[@]}"} "$BIN" $DADA_CMD "${filtF[@]}" \
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
