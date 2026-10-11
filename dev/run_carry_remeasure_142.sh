#!/usr/bin/env bash
# run_carry_remeasure_142.sh — re-measure the #139 carry's compare penalty on
# soil 16S (#142), with an R pool=TRUE reference built from the SAME reads.
#
# Every number in #142 predates the screen-timer change (447c1f3), the LOESS
# `interpolate` default (a42fd33) and the abundance-sorted pool (b703d64), and
# its 16S input was never length-normalised. This script rebuilds the
# measurement from one filtering pass, so the dada2-rs A/B arms and the R
# reference all read identical files.
#
# Input: the output of dev/cutadapt_v3v4.sh run WITH EQUALIZE_R1/_R2 set
# (<sample>_1.fastq.gz / <sample>_2.fastq.gz).
#
# Usage:
#   run_carry_remeasure_142.sh <dada2rs-binary> <cutadapt-out-dir> <out-dir>
#
# Stages (STAGES, space-separated, run in this order when listed):
#   prep    filter-and-trim once into <out>/filtered/<sample>{F,R}.fastq.gz,
#           then pin a training set (<out>/train.txt)
#   learn   dada2-rs learn-errors on the pinned set -> <out>/err/err{F,R}.json
#   ab      carry vs DADA2RS_SHUFFLE_NO_CARRY=1, both reads, REPS reps each,
#           sequential, one binary, NUMA policy fixed -> <out>/ab/
#   r       R reference: write_reference.R illumina --pool=true on the same
#           filtered reads and training set -> <out>/r/ (long; own job)
#   parity  dada2-rs dada-pooled under R's own error models, default and
#           DADA2RS_GAPLESS_X86=1, diffed against R's per-direction tables
#           -> <out>/parity/   (needs `r`)
#
# Env:
#   STAGES="prep learn ab"   default; `r` and `parity` are opt-in
#   TRUNC_LEN_F / TRUNC_LEN_R  REQUIRED for prep: the EQUALIZE_R1/_R2 values.
#                            Reads that truncQ cuts below this are dropped, so
#                            the pool stays one length per direction.
#   MAX_EE=2  TRUNC_Q=2  MAX_N=0   filter-and-trim, as run_illumina.sh
#   TRAIN_BASES=3e8          pinned training-set budget (pick_training_subset.py)
#   NBASES=1000000000000     learn-errors cap on BOTH sides (an integer: our
#                            --nbases rejects 1e12); must exceed the training
#                            set or R trains on a prefix of it
#   THREADS=64  REPS=2
#   NUMA_POLICY=interleave   see dev/numa_pin.sh
#   ERR_DIR=                 use err{F,R}.json from here for `ab` instead of
#                            <out>/err (e.g. <out>/parity/err for R's model)
#   R_THREADS=$THREADS       R's thread count (`r`)
#
# Run on an exclusive node and never time two arms concurrently (#281). Check
# the R side with the x86 caveat in mind: x86-64 R takes the gapless shortcut on
# some unequal-length pairs, which DADA2RS_GAPLESS_X86=1 reproduces (#277). With
# length-equalised input that should be rare, which the parity stage measures.
set -euo pipefail

BIN="${1:?usage: run_carry_remeasure_142.sh <binary> <cutadapt-out-dir> <out-dir>}"
IN="${2:?missing cutadapt-out-dir}"
OUT="${3:?missing out-dir}"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"

STAGES="${STAGES:-prep learn ab}"
MAX_EE="${MAX_EE:-2}"; TRUNC_Q="${TRUNC_Q:-2}"; MAX_N="${MAX_N:-0}"
TRAIN_BASES="${TRAIN_BASES:-3e8}"
NBASES="${NBASES:-1000000000000}"
THREADS="${THREADS:-64}"
REPS="${REPS:-2}"
R_THREADS="${R_THREADS:-$THREADS}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
FILT="$OUT/filtered"

has_stage() { [[ " $STAGES " == *" $1 "* ]]; }

# The binary's version tag, recorded with every stage so an arm built from a
# stale checkout cannot pass for this one (#132).
VERSION="$("$BIN" --version | awk '{print $NF}')"
mkdir -p "$OUT"
echo "binary: $BIN ($VERSION)" | tee -a "$OUT/versions.txt"

# The system a run came from. The re-measure is planned for a migrated cluster
# (new OS and SLURM, Ceph instead of GPFS). The leading account of #142 is
# #151's reused-buffer effect, which depends on glibc's mmap threshold and on
# page-fault behaviour, so glibc, kernel and THP are recorded with the numbers.
{
    echo "--- environment ($(date -u +%FT%TZ), stages: $STAGES)"
    echo "host: $(hostname)  kernel: $(uname -sr)"
    echo "glibc: $( (ldd --version 2>/dev/null || true) | head -1)"
    echo "thp: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null || echo n/a)"
    echo "numa: $( (numactl --hardware 2>/dev/null || true) | head -1)"
    echo "perf: $(command -v perf || echo absent)  paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo n/a)"
    echo "out fs: $( (df -PT "$OUT" 2>/dev/null || true) | awk 'NR==2 {print $2, $1}')"
    echo "slurm: job ${SLURM_JOB_ID:-none}, cpus ${SLURM_CPUS_ON_NODE:-n/a}"
} | tee -a "$OUT/versions.txt"

# Filtered inputs in one fixed order. Pooled ties break in first-seen order, so
# both tools must see the samples in the same order; R's list.files() sorts,
# and so does this (C locale; SRA run accessions collate the same either way).
filtered_list() {  # filtered_list F|R
    LC_ALL=C ls "$FILT"/*"$1".fastq.gz
}

# ---------------------------------------------------------------------------
if has_stage prep; then
    : "${TRUNC_LEN_F:?prep needs TRUNC_LEN_F (the EQUALIZE_R1 value)}"
    : "${TRUNC_LEN_R:?prep needs TRUNC_LEN_R (the EQUALIZE_R2 value)}"
    mkdir -p "$FILT"
    n=0
    for r1 in $(LC_ALL=C ls "$IN"/*_1.fastq.gz); do
        s=$(basename "$r1" _1.fastq.gz)
        r2="$IN/${s}_2.fastq.gz"
        [ -e "$r2" ] || { echo "missing mate for $s" >&2; exit 1; }
        "$BIN" filter-and-trim --fwd "$r1" --filt "$FILT/${s}F.fastq.gz" \
            --rev "$r2" --filt-rev "$FILT/${s}R.fastq.gz" \
            --trunc-len "$TRUNC_LEN_F" "$TRUNC_LEN_R" --max-n "$MAX_N" \
            --max-ee "$MAX_EE" "$MAX_EE" --trunc-q "$TRUNC_Q" --compress \
            > /dev/null
        n=$((n + 1))
    done
    echo "==> prep: $n pair(s) filtered into $FILT (trunc-len $TRUNC_LEN_F/$TRUNC_LEN_R)"
    python3 "$HERE/concordance/pick_training_subset.py" "$FILT" \
        --bases "$TRAIN_BASES" --suffix F.fastq.gz --out "$OUT/train.txt"
fi

# ---------------------------------------------------------------------------
if has_stage learn; then
    mkdir -p "$OUT/err"
    for d in F R; do
        train=()
        while IFS= read -r n; do train+=("$FILT/$n$d.fastq.gz"); done < "$OUT/train.txt"
        echo "==> learn-errors $d (${#train[@]} pinned samples)"
        "$BIN" learn-errors "${train[@]}" --nbases "$NBASES" --threads "$THREADS" \
            --verbose -o "$OUT/err/err$d.json" 2> "$OUT/err/learn$d.log" \
            || { tail -5 "$OUT/err/learn$d.log" >&2; exit 1; }
    done
fi

# ---------------------------------------------------------------------------
if has_stage ab; then
    ERR="${ERR_DIR:-$OUT/err}"
    # shellcheck source=dev/numa_pin.sh
    . "$HERE/numa_pin.sh"
    numa_pin_init "$THREADS"
    AB="$OUT/ab"; mkdir -p "$AB"
    echo "ab: error models from $ERR; numa: ${NUMA_PREFIX:-none}" | tee -a "$OUT/versions.txt"
    for rep in $(seq 1 "$REPS"); do
        # Alternate the arm order per rep so monotonic drift does not land on one arm.
        arms="carry nocarry"; [ $((rep % 2)) -eq 0 ] && arms="nocarry carry"
        for arm in $arms; do
            gate=(); want="shuffle carry=on"
            [ "$arm" = nocarry ] && { gate=(DADA2RS_SHUFFLE_NO_CARRY=1); want="shuffle carry=off"; }
            for d in F R; do
                o="$AB/${arm}_rep${rep}_$d"; mkdir -p "$o"
                inputs=(); while IFS= read -r f; do inputs+=("$f"); done < <(filtered_list "$d")
                echo "==> ab $arm rep $rep $d"
                env ${gate[@]+"${gate[@]}"} $NUMA_PREFIX "$BIN" dada-pooled "${inputs[@]}" \
                    --error-model "$ERR/err$d.json" -o "$o" --threads "$THREADS" \
                    --verbose --metrics-json "$o.metrics.json" 2> "$o.log" \
                    || { tail -5 "$o.log" >&2; exit 1; }
                # Positive control on the gate: refuse an arm that did not run
                # the configuration it is labelled with.
                grep -q "$want" "$o.metrics.json" \
                    || { echo "ERROR: $o did not run with '$want'" >&2; exit 1; }
            done
        done
    done
    python3 - "$AB" <<'PY' | tee "$AB/summary.tsv"
import glob, json, os, sys
ab = sys.argv[1]
cols = ["arm", "rep", "dir", "builds", "nclusters", "screened", "aligned",
        "run_dada", "compare", "map", "store", "shuffle",
        "screen_ns", "dp_ns", "wall"]
print("\t".join(cols))
for p in sorted(glob.glob(os.path.join(ab, "*.metrics.json"))):
    arm, rep, d = os.path.basename(p)[: -len(".metrics.json")].split("_")
    r = json.load(open(p))["runs"][0]
    c, ph, sh = r["compare"], r["phases"], r["shuffle"]
    row = [arm, rep, d, sh["builds"], r["run"]["nclusters"], c["screened"],
           c["aligned"], ph["loop_total"], ph["compare"], c["attribution"]["map"],
           c["attribution"]["store"], ph["shuffle"],
           c["split"]["screen_ns_per_screened"], c["split"]["dp_ns_per_aligned"],
           json.load(open(p))["wall_seconds"]]
    print("\t".join(f"{v:.1f}" if isinstance(v, float) else str(v) for v in row))
PY
    # The carry must not move a single ASV.
    for d in F R; do
        python3 "$HERE/compare_asvs.py" --baseline "nocarry=$AB/nocarry_rep1_$d" \
            --compare "carry=$AB/carry_rep1_$d" | tee "$AB/asvs_$d.txt"
    done
fi

# ---------------------------------------------------------------------------
if has_stage r; then
    mkdir -p "$OUT/r"
    echo "==> R reference (pool=TRUE) on $FILT"
    Rscript "$HERE/concordance/write_reference.R" illumina "$FILT" "$OUT/r/ref.csv" \
        --prefiltered --pool=true --train-samples="$OUT/train.txt" \
        --nbases="$NBASES" --threads="$R_THREADS" 2>&1 | tee "$OUT/r/log.txt"
    Rscript -e 'cat("R dada2", as.character(packageVersion("dada2")), R.version$platform, "\n")' \
        | tee -a "$OUT/versions.txt"
fi

# ---------------------------------------------------------------------------
if has_stage parity; then
    P="$OUT/parity"; mkdir -p "$P/err"
    for d in F R; do
        Rscript "$REPO/scripts/learnerrors_to_dada2rs.R" "$OUT/r/ref.err$d.rds" "$P/err/err$d.json"
        python3 "$HERE/concordance/rtab_to_seqtab.py" "$OUT/r/ref.dada$d.csv" "$P/r_dada$d.json"
    done
    # Same model and same reads on both sides, so a residual is denoising alone.
    # The x86 arm reproduces x86-64 R's gapless test (#277); if R ran on x86 it
    # is the arm that should match.
    for arm in default x86; do
        gate=(); [ "$arm" = x86 ] && gate=(DADA2RS_GAPLESS_X86=1)
        for d in F R; do
            o="$P/${arm}_$d"; mkdir -p "$o"
            inputs=(); while IFS= read -r f; do inputs+=("$f"); done < <(filtered_list "$d")
            echo "==> parity $arm $d"
            env ${gate[@]+"${gate[@]}"} "$BIN" dada-pooled "${inputs[@]}" \
                --error-model "$P/err/err$d.json" -o "$o" --threads "$THREADS" 2> "$o.log" \
                || { tail -5 "$o.log" >&2; exit 1; }
        done
    done
    for d in F R; do
        python3 "$HERE/compare_asvs.py" --baseline "R=$P/r_dada$d.json" \
            --compare "default=$P/default_$d" --compare "x86=$P/x86_$d" \
            | tee "$P/asvs_$d.txt"
    done
fi

echo "done: $OUT ($STAGES)"
