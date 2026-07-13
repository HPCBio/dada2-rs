#!/usr/bin/env bash
# run_band_sweep.sh
# ---------------------------------------------------------------------------
# Sweep the alignment band radius (--band) and measure its effect on the
# pooled ASV set and the pooled partition (which uniques are absorbed vs.
# split), holding the k-mer screen (--kmer-size / --kdist-cutoff) and every
# other parameter CONSTANT. The band is swept in BOTH error learning
# (errors-from-sample) and inference (dada-pooled), so the result is the full
# end-to-end effect of the BAND_SIZE default.
#
# WHY: banded alignment is a soft constraint, not a screen — a pair whose true
# optimal alignment needs more net gaps than the band is aligned SUB-OPTIMALLY
# (not rejected), and that alignment feeds the error-model transition tally.
# The kdist-calibration work found genuine error copies need band<=8 on both
# Illumina and PacBio HiFi; the large-band tail is truncation/off-target
# artifacts that no practical band reaches anyway. This sweep tests that
# empirically: does tightening the band change the ASVs or the partition?
#
# The pass criteria (evaluated by compare_bands.py against the FIRST band as
# baseline) are: (1) ASV set identical, (2) no MULTI-READ unique changes its
# assignment. Singleton churn is reported but benign.
#
# Usage:
#   bash run_band_sweep.sh <input> [outdir] [errfun] [kmer-size] [band-list]
#
#   <input> : a directory, a quoted glob, or a single derep/sample JSON (or
#             FASTQ). Multiple samples are denoised POOLED (dada-pooled).
#   band-list: space-separated; the FIRST is the baseline (use the shipped
#             default: 16 for Illumina, 32 for PacBio HiFi).
#
# Examples:
#   # Illumina MiSeq, baseline 16, tighten to 8 and 4:
#   bash run_band_sweep.sh "$prep/derep/R1" out_R1 loess 5 "16 8 4"
#
#   # PacBio HiFi, baseline 32, tighten to 16 and 8 (k=7):
#   bash run_band_sweep.sh /path/to/pacbio/derep out_pb pacbio 7 "32 16 8"
#
# Env overrides (held constant across all bands):
#   MAX_CONSIST=10  KDIST_CUTOFF=0.42  THREADS=1  GZIP=1
#   FILE_GLOB="*.json.gz *.json *.fastq.gz *.fastq"
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer an optimized build; fall back to debug.
DADA2RS=""
for cand in "${SCRIPT_DIR}/../target/release-native/dada2-rs" \
            "${SCRIPT_DIR}/../target/release/dada2-rs" \
            "${SCRIPT_DIR}/../target/debug/dada2-rs"; do
    [[ -x "$cand" ]] && { DADA2RS="$cand"; break; }
done
DADA2RS="${DADA2RS_BIN:-$DADA2RS}"

INPUT="${1:?Usage: run_band_sweep.sh <dir|glob|file> [outdir] [errfun] [kmer-size] [band-list]}"
OUTDIR="${2:-${SCRIPT_DIR}/band_sweep_out}"
ERRFUN="${3:-loess}"
KMER_SIZE="${4:-5}"
BANDLIST="${5:-16 8 4}"

MAX_CONSIST="${MAX_CONSIST:-10}"
KDIST_CUTOFF="${KDIST_CUTOFF:-0.42}"
THREADS="${THREADS:-1}"
GZIP="${GZIP:-1}"
FILE_GLOB="${FILE_GLOB:-*.json.gz *.json *.fastq.gz *.fastq}"

if [[ -z "$DADA2RS" || ! -x "$DADA2RS" ]]; then
    echo "ERROR: dada2-rs binary not found. Build it (cargo build --release) or set DADA2RS_BIN." >&2
    exit 1
fi
echo "==> binary: $DADA2RS"
"$DADA2RS" --version 2>/dev/null || true

# --- resolve inputs -------------------------------------------------------
declare -a INPUTS=()
read -r -a GLOB_PATS <<< "$FILE_GLOB"
if [[ -d "$INPUT" ]]; then
    shopt -s nullglob
    for pat in "${GLOB_PATS[@]}"; do for f in "$INPUT"/$pat; do INPUTS+=("$f"); done; done
    shopt -u nullglob
elif [[ -f "$INPUT" ]]; then
    INPUTS+=("$INPUT")
else
    shopt -s nullglob
    for f in $INPUT; do INPUTS+=("$f"); done
    shopt -u nullglob
fi
[[ ${#INPUTS[@]} -gt 1 ]] && { IFS=$'\n' INPUTS=($(printf '%s\n' "${INPUTS[@]}" | sort)); unset IFS; }
[[ ${#INPUTS[@]} -eq 0 ]] && { echo "ERROR: no inputs for '$INPUT' (glob: $FILE_GLOB)" >&2; exit 1; }

echo "==> ${#INPUTS[@]} sample(s); band list: [$BANDLIST]  (baseline = first)"

# R1/R2 mixing guard: pooling forward+reverse gives meaningless ASVs.
n_r1=0; n_r2=0
for f in "${INPUTS[@]}"; do
    b="$(basename "$f")"
    [[ "$b" == *_R1[._]* ]] && n_r1=$((n_r1+1))
    [[ "$b" == *_R2[._]* ]] && n_r2=$((n_r2+1))
done
if [[ $n_r1 -gt 0 && $n_r2 -gt 0 ]]; then
    echo "WARNING: inputs mix R1 ($n_r1) and R2 ($n_r2); denoise each orientation separately." >&2
    echo "         Continuing in 3s (Ctrl-C to abort)..." >&2; sleep 3
fi

mkdir -p "$OUTDIR"

# --- build the derep set (JSON passthrough; FASTQ -> derep) ----------------
DEREP_DIR="${OUTDIR}/derep"; mkdir -p "$DEREP_DIR"
declare -a DEREPS=()
for f in "${INPUTS[@]}"; do
    case "$f" in
        *.json|*.json.gz) DEREPS+=("$f") ;;
        *)
            base="$(basename "$f")"; base="${base%.gz}"; base="${base%.fastq}"
            dj="${DEREP_DIR}/${base}.derep.json"
            [[ -s "$dj" ]] || { echo "==> derep $f"; "$DADA2RS" derep "$f" -o "$dj"; }
            DEREPS+=("$dj") ;;
    esac
done

# --- per-band: learn errors + dada-pooled (+ pooled record) ----------------
ext="json"; [[ "$GZIP" == "1" ]] && ext="json.gz"
declare -a POOLED_ARGS=() ERR_ARGS=()
for BAND in $BANDLIST; do
    echo ""
    echo "================  band = $BAND  ================"
    ERR_JSON="${OUTDIR}/errors_band${BAND}.json"
    DADA_OUT="${OUTDIR}/dada_band${BAND}"
    POOLED_REC="${OUTDIR}/pooled_band${BAND}.${ext}"
    mkdir -p "$DADA_OUT"

    echo "==> errors-from-sample (band=$BAND, k=$KMER_SIZE, cutoff=$KDIST_CUTOFF, errfun=$ERRFUN)"
    "$DADA2RS" errors-from-sample "${DEREPS[@]}" \
        --errfun "$ERRFUN" --band "$BAND" \
        --kmer-size "$KMER_SIZE" --kdist-cutoff "$KDIST_CUTOFF" \
        --max-consist "$MAX_CONSIST" --threads "$THREADS" \
        -o "$ERR_JSON" 2> "${OUTDIR}/learn_band${BAND}.log"

    echo "==> dada-pooled (band=$BAND)"
    gzip_flag=(); [[ "$GZIP" == "1" ]] && gzip_flag=(--gzip)
    "$DADA2RS" dada-pooled "${DEREPS[@]}" \
        --error-model "$ERR_JSON" --output-dir "$DADA_OUT" \
        --pooled-record "$POOLED_REC" "${gzip_flag[@]}" \
        --band "$BAND" --kmer-size "$KMER_SIZE" --kdist-cutoff "$KDIST_CUTOFF" \
        --threads "$THREADS" 2> "${OUTDIR}/dada_band${BAND}.log"

    echo "    errors -> $ERR_JSON"
    echo "    pooled -> $POOLED_REC"
    if [[ ${#POOLED_ARGS[@]} -eq 0 ]]; then
        POOLED_ARGS+=(--baseline "${BAND}=${POOLED_REC}")
        ERR_ARGS+=(--errors "${BAND}=${ERR_JSON}")
    else
        POOLED_ARGS+=(--compare "${BAND}=${POOLED_REC}")
        ERR_ARGS+=(--errors "${BAND}=${ERR_JSON}")
    fi
done

# --- compare ---------------------------------------------------------------
echo ""
echo "================  COMPARE  ================"
python3 "${SCRIPT_DIR}/compare_bands.py" "${POOLED_ARGS[@]}" "${ERR_ARGS[@]}" \
    --json "${OUTDIR}/band_report.json" || true
echo ""
echo "Done. Outputs in ${OUTDIR}/ (report: band_report.json)"
