#!/usr/bin/env bash
# run_endmodel_sweep_289.sh — alignment end model x backend sweep (#289, #51).
#
# Runs the full bench_pooled.py pipeline once per arm and replicate:
#   backend  nw | wfa (affine, gap opening 0) | wfalin (gap-linear WFA)
#   model    ends-free (R) | free-3p | global
# Replicates are interleaved (rep 1 of every arm, then rep 2, ...) so drift over
# a long job spreads across arms instead of landing on one. Summarise with
# compare_endmodel_arms.py.
#
# Usage:
#   run_endmodel_sweep_289.sh <dada2rs-binary> <illumina|pacbio> <input-dir> <out-dir> [-- bench args...]
#
# Anything after `--` goes to bench_pooled.py unchanged (e.g. --trunc-len 240,160,
# --fwd-pattern _R1_001, --kmer-size 7, --pool false).
#
# Env:
#   BACKENDS="nw wfa wfalin"   arms to run; wfa/wfalin need a --features wfa build
#   END_MODELS=                default: pacbio "ends-free free-3p global",
#                              illumina "ends-free free-3p" (global only fits
#                              reads that all run through the reverse primer)
#   REPS=2  THREADS=22  NUMA=interleave  POOL=true
#   WFA_MAX_EDITS=             WFA cap; unset = the binary's default (50)
#
# Run on an exclusive node, from node-local scratch, and never time two arms
# concurrently (#281).
set -euo pipefail

BIN="${1:?usage: run_endmodel_sweep_289.sh <binary> <illumina|pacbio> <input-dir> <out-dir> [-- bench args]}"
PLATFORM="${2:?missing platform}"
IN="${3:?missing input-dir}"
OUT="${4:?missing out-dir}"
shift 4
[ "${1:-}" = "--" ] && shift
BENCH_EXTRA=("$@")

case "$PLATFORM" in
    pacbio) default_models="ends-free free-3p global" ;;
    illumina) default_models="ends-free free-3p" ;;
    *) echo "platform must be illumina or pacbio" >&2; exit 1 ;;
esac
BACKENDS="${BACKENDS:-nw wfa wfalin}"
END_MODELS="${END_MODELS:-$default_models}"
REPS="${REPS:-2}"
THREADS="${THREADS:-22}"
NUMA="${NUMA:-interleave}"
POOL="${POOL:-true}"
WFA_MAX_EDITS="${WFA_MAX_EDITS:-}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
mkdir -p "$OUT"

# Appended, so a rerun into the same directory (which skips finished arms)
# keeps the record of every invocation.
{
    echo "=== $(date -u +%FT%TZ)"
    echo "binary: $("$BIN" --version)"
    echo "platform: $PLATFORM  input: $IN"
    echo "backends: $BACKENDS  end models: $END_MODELS  reps: $REPS"
    echo "threads: $THREADS  numa: $NUMA  pool: $POOL  wfa_max_edits: ${WFA_MAX_EDITS:-default}"
    echo "bench extra: ${BENCH_EXTRA[*]:-}"
    echo "host: $(hostname)"
} | tee -a "$OUT/sweep_env.txt"

for rep in $(seq 1 "$REPS"); do
    for backend in $BACKENDS; do
        for model in $END_MODELS; do
            arm="$backend.$model"
            dir="$OUT/$arm/rep$rep"
            if [ -f "$dir/summary.csv" ]; then
                echo "skip $arm rep$rep (done)"
                continue
            fi
            args=(--dada2rs "$BIN" --outdir "$dir" --pool "$POOL" --threads "$THREADS"
                  --numa "$NUMA" --verbose --end-model "$model")
            envs=()
            case "$backend" in
                nw) ;;
                wfa) args+=(--align-backend wfa2) ;;
                wfalin) args+=(--align-backend wfa2); envs=(DADA2RS_WFA_LINEAR=1) ;;
                *) echo "unknown backend $backend" >&2; exit 1 ;;
            esac
            if [ "$backend" != nw ] && [ -n "$WFA_MAX_EDITS" ]; then
                args+=(--wfa-max-edits "$WFA_MAX_EDITS")
            fi
            echo "==> $arm rep$rep"
            start=$(date +%s)
            env ${envs[@]+"${envs[@]}"} python3 -I "$HERE/bench_pooled.py" "$PLATFORM" "$IN" \
                "${args[@]}" ${BENCH_EXTRA[@]+"${BENCH_EXTRA[@]}"} > "$OUT/$arm.rep$rep.log" 2>&1
            echo "    $(( $(date +%s) - start ))s"
        done
    done
done

echo "done; summarise with:"
echo "  python3 $HERE/compare_endmodel_arms.py $OUT"
