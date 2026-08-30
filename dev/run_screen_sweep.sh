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
KS="${KS:-8 9}"
CUTS="${CUTS:-0.55 0.58 0.60 0.62 0.65 0.70 0.72}"

mkdir -p "$OUT"

aligned_count() {  # sum (nalign - nshroud) over a run's forward dada outputs
  python3 - "$1" <<'PY'
import glob, json, os, sys
na = ns = 0
for p in glob.glob(os.path.join(sys.argv[1], "dada_fwd", "*.json")):
    s = json.load(open(p))["stats"]
    na += s["nalign"]; ns += s["nshroud"]
print(f"{na} {ns} {na - ns}")
PY
}

echo "==> baseline: k-mer screen (production default)"
[ -d "$OUT/base" ] || bash "$RUN" "$BIN" "$DATA" "$OUT/base" "$THREADS" > "$OUT/base.log" 2>&1

echo "==> control: k-mer screen AGAIN (establishes the ASV noise floor)"
[ -d "$OUT/control" ] || bash "$RUN" "$BIN" "$DATA" "$OUT/control" "$THREADS" > "$OUT/control.log" 2>&1
echo "    control vs baseline (MUST be identical, or nothing below is interpretable):"
python3 "$HERE/compare_seqtab_matrix.py" \
    "$OUT/base/seqtab.nochim.json" "$OUT/control/seqtab.nochim.json" \
    --label-a base --label-b control | tail -2 | sed 's/^/      /'

echo
echo "==> sweeping the (k, cutoff) grid"
for K in $KS; do
  for C in $CUTS; do
    d="$OUT/mini_k${K}_c${C}"
    [ -d "$d" ] && { echo "    k=$K cutoff=$C (cached)"; continue; }
    echo "    k=$K cutoff=$C"
    SCREEN_BACKEND=minimizer MINIMIZER_K="$K" SCREEN_CUTOFF="$C" \
      bash "$RUN" "$BIN" "$DATA" "$d" "$THREADS" > "$d.log" 2>&1
  done
done

echo
echo "======================= ACCURACY + WORK ======================="
read -r bna bns bal <<< "$(aligned_count "$OUT/base")"
printf "%-22s %8s %8s %12s %10s %s\n" "arm" "ASVs" "churn" "aligned" "vs kmer" "count-matrix L1"
printf "%-22s %8s %8s %12s %10s %s\n" "k-mer (baseline)" \
    "$(python3 -c "import json;print(len(json.load(open('$OUT/base/seqtab.nochim.json'))['sequences']))")" \
    "-" "$bal" "100.0%" "-"
for K in $KS; do
  for C in $CUTS; do
    d="$OUT/mini_k${K}_c${C}"
    [ -d "$d" ] || continue
    read -r na ns al <<< "$(aligned_count "$d")"
    n=$(python3 -c "import json;print(len(json.load(open('$d/seqtab.nochim.json'))['sequences']))")
    churn=$(python3 "$HERE/compare_asvs.py" --baseline k="$OUT/base/seqtab.nochim.json" \
              --compare m="$d/seqtab.nochim.json" 2>/dev/null | grep -o "churn=[0-9]*" | head -1 | cut -d= -f2)
    l1=$(python3 "$HERE/compare_seqtab_matrix.py" "$OUT/base/seqtab.nochim.json" \
              "$d/seqtab.nochim.json" --label-a k --label-b m 2>/dev/null \
           | grep -o "L1 = [0-9]* reads ([0-9.]*%" | grep -o "([0-9.]*%" | tr -d '(')
    printf "%-22s %8s %8s %12s %9s%% %s\n" "mini k=$K cut=$C" "$n" "${churn:-?}" "$al" \
        "$(python3 -c "print(f'{100*$al/$bal:.1f}')")" "${l1:-?}"
  done
done

echo
echo "========================== WALL TIME =========================="
echo "Arms interleaved across $REPS reps; the k-mer arm appears twice so the"
echo "control channel carries the same thermal/cache exposure as the test arms."
echo "Timing the dada step only -- filter/merge/chimera are identical across arms."

filtF=("$OUT"/base/filtered/*_F_filt.fastq.gz)
declare -a ARMS=("kmer::" "kmerctl::")
for K in $KS; do for C in $CUTS; do ARMS+=("mini_k${K}_c${C}:$K:$C"); done; done

: > "$OUT/timings.tsv"
for rep in $(seq 1 "$REPS"); do
  for arm in "${ARMS[@]}"; do
    name="${arm%%:*}"; rest="${arm#*:}"; K="${rest%%:*}"; C="${rest#*:}"
    extra=()
    [ -n "$K" ] && extra=(--screen-backend minimizer --minimizer-k "$K" --kdist-cutoff "$C")
    t0=$(python3 -c 'import time;print(time.time())')
    "$BIN" dada "${filtF[@]}" --error-model "$OUT/base/errF.json" \
        --output-dir "$OUT/.timing" --threads "$THREADS" \
        ${extra[@]+"${extra[@]}"} > /dev/null 2>&1
    t1=$(python3 -c 'import time;print(time.time())')
    printf "%s\t%s\t%s\n" "$name" "$rep" \
        "$(python3 -c "print(f'{$t1-$t0:.2f}')")" >> "$OUT/timings.tsv"
  done
done

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
