#!/usr/bin/env bash
# run_cutoff_sweep.sh
# ---------------------------------------------------------------------------
# Sweep the k-mer pre-alignment screen cutoff (--kdist-cutoff) and measure its
# effect on the pooled ASV set and partition, holding --band, --kmer-size and
# every other parameter CONSTANT. The cutoff is swept in BOTH error learning
# (learn-errors) and inference (dada-pooled), so the result is the full
# end-to-end effect of the KDIST_CUTOFF default (0.42).
#
# WHY: unlike the band (a soft alignment constraint), the k-mer cutoff is a HARD
# screen — a pair whose k-mer distance exceeds it is never aligned, so it can
# never be linked as an error copy. Lowering the cutoff prunes more pairs
# (faster), but below the real error-copy distance it starts screening out
# genuine error copies, which then fail to be corrected (over-splitting / more
# failures). The kdist-calibration work put real error copies far below 0.42
# (median kdist ~0.02, clear-error-copy ceiling ~0.10-0.16 on both platforms),
# so there is large nominal headroom; this sweep tests whether tightening the
# cutoff actually changes the ASVs.
#
# Companion to run_band_sweep.sh; same structure, same comparison tool
# (compare_bands.py). Pass criteria vs the FIRST cutoff (baseline): ASV-set
# identity AND reads-moved-% under the ceiling.
#
# Usage:
#   bash run_cutoff_sweep.sh <input> [outdir] [errfun] [kmer-size] [band] [cutoff-list]
#
#   <input>      : directory, quoted glob, or a single derep/sample JSON (or FASTQ).
#   band         : held constant (use the platform default: 16 Illumina, 32 HiFi).
#   cutoff-list  : space-separated; FIRST is the baseline (0.42 = shipped default).
#
# Examples:
#   # Illumina MiSeq, band 16, tighten cutoff 0.42 -> 0.30 -> 0.20:
#   bash run_cutoff_sweep.sh "$prep/derep/R1" out_R1 loess 5 16 "0.42 0.30 0.20"
#
#   # PacBio HiFi, band 32, k=7:
#   bash run_cutoff_sweep.sh /path/to/pacbio/derep out_pb pacbio 7 32 "0.42 0.30 0.20"
#
# Env overrides (held constant across all cutoffs):
#   MAX_CONSIST=10  THREADS=1  GZIP=1  NBASES=100000000
#   MAX_READS_MOVED_PCT=0.5   # OVERALL pass ceiling on reads reassigned
#   LEARN_CUTOFF=""           # decouple learn cutoff from inference (see below)
#   FILE_GLOB="*.json.gz *.json *.fastq.gz *.fastq"
#
# DECOUPLED MODE (LEARN_CUTOFF set): learn the error model ONCE at the looser
# LEARN_CUTOFF, then sweep only the INFERENCE cutoff reusing that fixed model.
# This isolates the inference-screen effect (partition / fragmentation) from the
# error-model destabilization that joint tightening incurs — i.e. "looser KDIST
# for error modeling, tighter for inference". Example:
#   LEARN_CUTOFF=0.42 bash run_cutoff_sweep.sh "$prep/derep/R1" out_dec loess 5 16 "0.42 0.30 0.20"
#
# The cutoff is a HARD screen, so tightening it reassigns far more reads than the
# band (mostly singletons screened out). The pass ceiling on reads-moved defaults
# to 0.5% here (vs the band sweep's tighter implicit 0.01%); ASV-set identity is
# always required regardless. Read the ASV-identity line as the primary verdict.
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DADA2RS=""
for cand in "${SCRIPT_DIR}/../target/release-native/dada2-rs" \
            "${SCRIPT_DIR}/../target/release/dada2-rs" \
            "${SCRIPT_DIR}/../target/debug/dada2-rs"; do
    [[ -x "$cand" ]] && { DADA2RS="$cand"; break; }
done
DADA2RS="${DADA2RS_BIN:-$DADA2RS}"

INPUT="${1:?Usage: run_cutoff_sweep.sh <dir|glob|file> [outdir] [errfun] [kmer-size] [band] [cutoff-list]}"
OUTDIR="${2:-${SCRIPT_DIR}/cutoff_sweep_out}"
ERRFUN="${3:-loess}"
KMER_SIZE="${4:-5}"
BAND="${5:-16}"
CUTOFFLIST="${6:-0.42 0.30 0.20}"

MAX_CONSIST="${MAX_CONSIST:-10}"
THREADS="${THREADS:-1}"
GZIP="${GZIP:-1}"
NBASES="${NBASES:-100000000}"
MAX_READS_MOVED_PCT="${MAX_READS_MOVED_PCT:-0.5}"
# DECOUPLE the learn cutoff from the inference cutoff. Empty = joint (learn AND
# infer at each swept cutoff, the default). Set LEARN_CUTOFF (e.g. 0.42) to learn
# the error model ONCE at that looser cutoff and reuse it for every inference
# cutoff — isolating the inference-screen effect (partition / fragmentation) from
# the error-model destabilization that joint tightening also incurs.
LEARN_CUTOFF="${LEARN_CUTOFF:-}"
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

echo "==> ${#INPUTS[@]} sample(s); band=$BAND (fixed); cutoff list: [$CUTOFFLIST]  (baseline = first)"

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

# --- optional: learn ONCE at a fixed (looser) cutoff -----------------------
ext="json"; [[ "$GZIP" == "1" ]] && ext="json.gz"
SHARED_ERR=""
if [[ -n "$LEARN_CUTOFF" ]]; then
    SHARED_ERR="${OUTDIR}/errors_learn${LEARN_CUTOFF}.json"
    echo ""
    echo "==> DECOUPLED: learn-errors ONCE at cutoff=$LEARN_CUTOFF (reused for all inference cutoffs)"
    "$DADA2RS" learn-errors "${DEREPS[@]}" \
        --errfun "$ERRFUN" --band "$BAND" \
        --kmer-size "$KMER_SIZE" --kdist-cutoff "$LEARN_CUTOFF" \
        --nbases "$NBASES" \
        --max-consist "$MAX_CONSIST" --threads "$THREADS" \
        --verbose -o "$SHARED_ERR" 2> "${OUTDIR}/learn_shared.log"
fi

# --- per-cutoff: (learn) + dada-pooled (+ pooled record) -------------------
declare -a POOLED_ARGS=() ERR_ARGS=()
for CUT in $CUTOFFLIST; do
    echo ""
    echo "================  cutoff = $CUT  ================"
    tag="cut${CUT}"
    ERR_JSON="${OUTDIR}/errors_${tag}.json"
    DADA_OUT="${OUTDIR}/dada_${tag}"
    POOLED_REC="${OUTDIR}/pooled_${tag}.${ext}"
    mkdir -p "$DADA_OUT"

    if [[ -n "$SHARED_ERR" ]]; then
        ERR_JSON="$SHARED_ERR"
        echo "==> using shared error model (learned at cutoff=$LEARN_CUTOFF)"
    else
        echo "==> learn-errors (cutoff=$CUT, band=$BAND, k=$KMER_SIZE, errfun=$ERRFUN, nbases=$NBASES)"
        "$DADA2RS" learn-errors "${DEREPS[@]}" \
            --errfun "$ERRFUN" --band "$BAND" \
            --kmer-size "$KMER_SIZE" --kdist-cutoff "$CUT" \
            --nbases "$NBASES" \
            --max-consist "$MAX_CONSIST" --threads "$THREADS" \
            --verbose -o "$ERR_JSON" 2> "${OUTDIR}/learn_${tag}.log"
    fi

    echo "==> dada-pooled (inference cutoff=$CUT)"
    gzip_flag=(); [[ "$GZIP" == "1" ]] && gzip_flag=(--gzip)
    # --verbose emits the "ALIGN: N aligns, M shrouded" line to the log.
    "$DADA2RS" dada-pooled "${DEREPS[@]}" \
        --error-model "$ERR_JSON" --output-dir "$DADA_OUT" \
        --pooled-record "$POOLED_REC" "${gzip_flag[@]}" \
        --band "$BAND" --kmer-size "$KMER_SIZE" --kdist-cutoff "$CUT" \
        --threads "$THREADS" --verbose 2> "${OUTDIR}/dada_${tag}.log"

    echo "    errors -> $ERR_JSON"
    echo "    pooled -> $POOLED_REC"
    if [[ ${#POOLED_ARGS[@]} -eq 0 ]]; then
        POOLED_ARGS+=(--baseline "${CUT}=${POOLED_REC}")
    else
        POOLED_ARGS+=(--compare "${CUT}=${POOLED_REC}")
    fi
    ERR_ARGS+=(--errors "${CUT}=${ERR_JSON}")
done

# --- compare ---------------------------------------------------------------
echo ""
echo "================  COMPARE  ================"
python3 "${SCRIPT_DIR}/compare_bands.py" "${POOLED_ARGS[@]}" "${ERR_ARGS[@]}" \
    --max-reads-moved-pct "$MAX_READS_MOVED_PCT" \
    --json "${OUTDIR}/cutoff_report.json" || true
echo ""
echo "Done. Outputs in ${OUTDIR}/ (report: cutoff_report.json)"
