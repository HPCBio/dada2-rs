#!/usr/bin/env bash
# run_shuffle_ab.sh
# ---------------------------------------------------------------------------
# A/B the two-level bounded shuffle build (#124) against a baseline build,
# through the standard bench_pooled.py harness so wall time, CPU, effective
# cores and peak RSS all come from the same instrumented source as every other
# perf result in docs/results.md.
#
# WHY a wrapper: the change is expected to be *byte-identical* and to move only
# wall time, and its payoff is workload-conditional on cluster count — the same
# axis that flipped #87's sign between two amplicons of one run. So a useful run
# has to produce three things together, or it cannot be interpreted:
#
#   1. wall/RSS deltas          (compare_bench.py)
#   2. proof the ASVs did not move  (compare_asvs.py — churn MUST be 0)
#   3. the realised prune rate  (the `shuffle build prune` verbose line)
#
# (3) is the one that makes a null result readable. The bound only pays if it
# prunes; if `examined %` comes back near 100 the phase never got cheaper and
# any wall delta is noise or overhead, not evidence about the design.
#
# Usage:
#   bash run_shuffle_ab.sh <platform> <input-dir> <baseline-bin> <branch-bin> [outdir] [threads]
#
#   <platform>     : illumina | pacbio  (passed straight to bench_pooled.py)
#   <baseline-bin> : dada2-rs built from the comparison point (e.g. main)
#   <branch-bin>   : dada2-rs built from this branch
#
# Build BOTH binaries the same way (target/release, not release-native, if the
# numbers are to be compared across nodes) — the build target moves wall time by
# more than this change is expected to.
#
# Env overrides:
#   REPS=1                 # runs per arm; >=3 gives compare_bench.py a median +- spread
#   HOT_SWEEP=""           # e.g. "8 32 128" — extra branch arms at each
#                          #   DADA2RS_SHUFFLE_HOT_CLUSTERS value. Free axis: it
#                          #   changes only the bound's tightness, never the
#                          #   partition, so every arm must still show churn=0.
#   PRIMER_FWD=/PRIMER_REV # required for pacbio (raw, primered reads)
#   EXTRA=""               # extra args forwarded to bench_pooled.py verbatim
#
# Examples:
#   # MiSeq, 3 reps per arm:
#   REPS=3 bash run_shuffle_ab.sh illumina /path/to/miseq_raw \
#       build/main/dada2-rs build/branch/dada2-rs out_shuf 32
#
#   # PacBio HiFi + hot-cluster sweep:
#   PRIMER_FWD=AGRGTTYGATYMTGGCTCAG PRIMER_REV=RGYTACCTTGTTACGACTT \
#   HOT_SWEEP="8 32 128" bash run_shuffle_ab.sh pacbio /path/to/hifi_raw \
#       build/main/dada2-rs build/branch/dada2-rs out_shuf_pb 32
set -euo pipefail

PLATFORM="${1:?usage: run_shuffle_ab.sh <platform> <input-dir> <baseline-bin> <branch-bin> [outdir] [threads]}"
INPUT="${2:?missing input-dir}"
BASE_BIN="${3:?missing baseline binary}"
BRANCH_BIN="${4:?missing branch binary}"
OUTDIR="${5:-shuffle_ab_out}"
THREADS="${6:-8}"

REPS="${REPS:-1}"
HOT_SWEEP="${HOT_SWEEP:-}"
EXTRA="${EXTRA:-}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$HERE/benchmark/bench_pooled.py"

platform_args=()
if [ "$PLATFORM" = "pacbio" ]; then
    platform_args+=(--primer-fwd "${PRIMER_FWD:?pacbio needs PRIMER_FWD}")
    platform_args+=(--primer-rev "${PRIMER_REV:?pacbio needs PRIMER_REV}")
fi

mkdir -p "$OUTDIR"

# One arm = one (label, binary, hot-cluster override) triple, run REPS times
# into its own directory so compare_bench.py can take the median.
run_arm() {
    local label="$1" bin="$2" hot="${3:-}"
    local dirs=()
    for rep in $(seq 1 "$REPS"); do
        local d="$OUTDIR/${label}_rep${rep}"
        echo "=== arm $label rep $rep -> $d" >&2
        # --verbose is what puts the `shuffle build prune` / `shuffle scan`
        # lines into the per-step logs; without it this run cannot be read.
        env ${hot:+DADA2RS_SHUFFLE_HOT_CLUSTERS="$hot"} \
            python3 "$BENCH" "$PLATFORM" "$INPUT" \
                --outdir "$d" --threads "$THREADS" --verbose \
                --dada2rs "$bin" ${platform_args[@]+"${platform_args[@]}"} $EXTRA \
            > "$d.stdout" 2>&1 || { echo "ARM $label rep $rep FAILED — see $d.stdout" >&2; tail -20 "$d.stdout" >&2; exit 1; }
        dirs+=("$d")
    done
    # Comma-separated list = reps, which compare_bench.py reduces to a median.
    # This is run_arm's ONLY stdout — progress goes to stderr above, or it would
    # be swallowed by the caller's command substitution.
    (IFS=,; echo "${dirs[*]}")
}

BASE_DIRS="$(run_arm base "$BASE_BIN")"
BRANCH_DIRS="$(run_arm branch "$BRANCH_BIN")"

compare_args=(--baseline "base=$BASE_DIRS" --compare "branch=$BRANCH_DIRS")
asv_args=(--baseline "base=$(echo "$BASE_DIRS" | cut -d, -f1)/rust/dada_fwd")
asv_args+=(--compare "branch=$(echo "$BRANCH_DIRS" | cut -d, -f1)/rust/dada_fwd")

for hot in $HOT_SWEEP; do
    d="$(run_arm "hot$hot" "$BRANCH_BIN" "$hot")"
    compare_args+=(--compare "hot$hot=$d")
    asv_args+=(--compare "hot$hot=$(echo "$d" | cut -d, -f1)/rust/dada_fwd")
done

echo
echo "############################################################"
echo "# 1. wall / cpu / cores / peak RSS"
echo "############################################################"
python3 "$HERE/benchmark/compare_bench.py" "${compare_args[@]}" | tee "$OUTDIR/compare_bench.txt"

echo
echo "############################################################"
echo "# 2. ASV identity — churn MUST be 0 in every arm"
echo "############################################################"
python3 "$HERE/compare_asvs.py" "${asv_args[@]}" | tee "$OUTDIR/compare_asvs.txt"

echo
echo "############################################################"
echo "# 3. shuffle scan split / prune rate (from the denoise-step logs)"
echo "############################################################"
{
    for d in "$OUTDIR"/*_rep1; do
        [ -d "$d" ] || continue
        echo "--- $(basename "$d")"
        grep -hE "shuffle scan split|shuffle build prune|shuffle scan time|phase times" \
            "$d"/rust/dada_*.log 2>/dev/null || echo "  (no shuffle lines — was --verbose honoured?)"
    done
} | tee "$OUTDIR/shuffle_lines.txt"

echo
echo "Wrote $OUTDIR/{compare_bench,compare_asvs,shuffle_lines}.txt"
echo "Read them together: a wall delta is only interpretable next to the prune rate,"
echo "and only valid if ASV churn is 0."
