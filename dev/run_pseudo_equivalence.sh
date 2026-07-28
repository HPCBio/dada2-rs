#!/usr/bin/env bash
# run_pseudo_equivalence.sh
# ---------------------------------------------------------------------------
# Test whether `dada-pseudo` is equivalent to a hand-rolled two-round run:
#
#   arm A (built-in) : dada-pseudo
#   arm B (manual)   : dada  ->  make-sequence-table  ->  seq-table-to-fasta
#                      --prevalence N  ->  dada --prior
#
# In dada2-rs these two SHOULD be identical. `dada-pseudo`'s round 2 is
# `mark_priors` + `dada_uniques`, which is exactly what `dada --prior` does:
# same prior-marking function, same core, and one error model resolved once and
# used for both rounds. Arm B reproduces the prior selection with
# `seq-table-to-fasta --prevalence`, which calls the very same
# `select_sequences` that dada-pseudo uses internally.
#
# So a difference here is a bug on our side. Issue #100.
#
# Why this matters beyond an internal check: the manual two-round form is what
# you need to denoise samples in parallel across a cluster (round 1 per sample,
# collect priors centrally, round 2 per sample). This verifies that workflow
# rather than assuming it. It is also the baseline for the separate, still
# UNCONFIRMED question in #100 of whether R's `pool="pseudo"` re-estimates its
# error model between rounds -- which, if true, would make R's built-in and R's
# manual form differ even though ours do not.
#
# Both arms are fed the SAME dereplicated JSONs and the SAME error model, so the
# only thing that can differ is the pseudo/manual plumbing itself.
#
# Usage:
#   bash run_pseudo_equivalence.sh <input> [outdir] [errfun] [band] [prevalence]
#
#   <input> may be:
#     - a directory     : all *.fastq / *.fastq.gz (or *.json / *.json.gz) inside
#     - a glob/file list: quote it, e.g. "data/*.fastq.gz"
#
#   At least 2 samples are required: with one sample nothing can reach
#   prevalence 2, the prior set is empty, and both arms trivially reduce to a
#   plain per-sample run.
#
# Examples:
#   bash run_pseudo_equivalence.sh "data/dada2/sam[12]F.fastq.gz" ./pseudo_eq_out loess 16
#   bash run_pseudo_equivalence.sh /path/to/fastq_dir ./pseudo_eq_out loess 16 2
#
# Defaults: outdir=./pseudo_equivalence_out  errfun=loess  band=16  prevalence=2
#           (prevalence 2 mirrors R DADA2's PSEUDO_PREVALENCE; the abundance
#            disjunct is left off, mirroring R's PSEUDO_ABUNDANCE=Inf)
#
# Environment overrides (held CONSTANT across both arms):
#   THREADS=1  KDIST_CUTOFF=0.42  MAX_CONSIST=10  NBASES=100000000
#   FILE_GLOB="*.fastq.gz *.fastq *.json *.json.gz"
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer an optimized build; fall back across profiles. DADA2RS_BIN overrides.
DADA2RS=""
for cand in "${SCRIPT_DIR}/../target/release-native/dada2-rs" \
            "${SCRIPT_DIR}/../target/release/dada2-rs" \
            "${SCRIPT_DIR}/../target/debug/dada2-rs"; do
    [[ -x "$cand" ]] && { DADA2RS="$cand"; break; }
done
DADA2RS="${DADA2RS_BIN:-$DADA2RS}"

INPUT="${1:?Usage: run_pseudo_equivalence.sh <dir|glob> [outdir] [errfun] [band] [prevalence]}"
OUTDIR="${2:-${SCRIPT_DIR}/pseudo_equivalence_out}"
ERRFUN="${3:-loess}"
BAND="${4:-16}"
PREVALENCE="${5:-2}"

THREADS="${THREADS:-1}"
KDIST_CUTOFF="${KDIST_CUTOFF:-0.42}"
MAX_CONSIST="${MAX_CONSIST:-10}"
NBASES="${NBASES:-100000000}"
FILE_GLOB="${FILE_GLOB:-*.fastq.gz *.fastq *.json *.json.gz}"

if [[ ! -x "$DADA2RS" ]]; then
    echo "ERROR: binary not found/executable at $DADA2RS" >&2
    echo "       Run 'cargo build --release' first." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Step 0: resolve <input> into an array of sample files.
# ---------------------------------------------------------------------------
declare -a INPUTS=()
read -r -a GLOB_PATS <<< "$FILE_GLOB"
if [[ -d "$INPUT" ]]; then
    shopt -s nullglob
    for pat in "${GLOB_PATS[@]}"; do
        for f in "$INPUT"/$pat; do INPUTS+=("$f"); done
    done
    shopt -u nullglob
elif [[ -f "$INPUT" ]]; then
    INPUTS+=("$INPUT")
else
    shopt -s nullglob
    for f in $INPUT; do INPUTS+=("$f"); done
    shopt -u nullglob
fi
if [[ ${#INPUTS[@]} -gt 1 ]]; then
    IFS=$'\n' INPUTS=($(printf '%s\n' "${INPUTS[@]}" | sort)); unset IFS
fi

if [[ ${#INPUTS[@]} -lt 2 ]]; then
    echo "ERROR: need at least 2 samples (got ${#INPUTS[@]}); with one sample the" >&2
    echo "       prior set is empty and both arms reduce to a plain dada run." >&2
    exit 1
fi
echo "==> ${#INPUTS[@]} sample(s)"

# ---------------------------------------------------------------------------
# Step 1: dereplicate once. BOTH arms consume these exact files, so the arms
# cannot differ because of a re-read or a re-derep.
# ---------------------------------------------------------------------------
DEREP_DIR="${OUTDIR}/derep"
mkdir -p "$DEREP_DIR"
declare -a DEREPS=()

for f in "${INPUTS[@]}"; do
    case "$f" in
        *.json|*.json.gz)
            DEREPS+=("$f")
            ;;
        *)
            base="$(basename "$f")"
            base="${base%.gz}"; base="${base%.fastq}"; base="${base%.fq}"
            dj="${DEREP_DIR}/${base}.derep.json"
            if [[ ! -s "$dj" ]]; then
                echo "==> derep $base"
                "$DADA2RS" derep "$f" -o "$dj" >/dev/null
            fi
            DEREPS+=("$dj")
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Step 2: learn the error model ONCE. Shared by both arms and, within arm B, by
# both rounds -- matching dada-pseudo, which resolves its params a single time.
# ---------------------------------------------------------------------------
ERR_JSON="${OUTDIR}/err.json"
if [[ ! -s "$ERR_JSON" ]]; then
    echo "==> learn-errors"
    "$DADA2RS" learn-errors "${DEREPS[@]}" \
        --errfun "$ERRFUN" --band "$BAND" \
        --kdist-cutoff "$KDIST_CUTOFF" \
        --nbases "$NBASES" \
        --max-consist "$MAX_CONSIST" --threads "$THREADS" \
        -o "$ERR_JSON" 2> "${OUTDIR}/learn_errors.log"
fi

COMMON=(--error-model "$ERR_JSON" --band "$BAND"
        --kdist-cutoff "$KDIST_CUTOFF" --threads "$THREADS")

# ---------------------------------------------------------------------------
# Arm A: the built-in two-round path, dumping the prior set it derived.
# ---------------------------------------------------------------------------
echo "==> arm A: dada-pseudo"
rm -rf "${OUTDIR}/A_pseudo"
"$DADA2RS" dada-pseudo "${DEREPS[@]}" \
    "${COMMON[@]}" \
    --pseudo-prevalence "$PREVALENCE" \
    --priors-out "${OUTDIR}/priors.A.fa" \
    -o "${OUTDIR}/A_pseudo" > "${OUTDIR}/A_pseudo.log" 2>&1

# ---------------------------------------------------------------------------
# Arm B: the manual form. Round 1 is a plain per-sample run; priors come from
# the round-1 sequence table via the same selection rule; round 2 flags them.
# ---------------------------------------------------------------------------
echo "==> arm B round 1: dada (no priors)"
rm -rf "${OUTDIR}/B_round1"
"$DADA2RS" dada "${DEREPS[@]}" "${COMMON[@]}" \
    --output-dir "${OUTDIR}/B_round1" > "${OUTDIR}/B_round1.log" 2>&1

echo "==> arm B: sequence table + prior selection (prevalence >= ${PREVALENCE})"
"$DADA2RS" make-sequence-table "${OUTDIR}"/B_round1/*.json \
    -o "${OUTDIR}/B_seqtab.r1.json" > /dev/null
"$DADA2RS" seq-table-to-fasta "${OUTDIR}/B_seqtab.r1.json" \
    --prevalence "$PREVALENCE" \
    -o "${OUTDIR}/priors.B.fa" > /dev/null

echo "==> arm B round 2: dada --prior"
rm -rf "${OUTDIR}/B_round2"
"$DADA2RS" dada "${DEREPS[@]}" "${COMMON[@]}" \
    --prior "${OUTDIR}/priors.B.fa" \
    --output-dir "${OUTDIR}/B_round2" > "${OUTDIR}/B_round2.log" 2>&1

# ---------------------------------------------------------------------------
# Arm C (diagnostic): only meaningful if the two prior sets differ. Re-runs
# round 2 with ARM A's priors, which separates "prior selection differs" from
# "denoising given identical priors differs".
# ---------------------------------------------------------------------------
if ! "${SCRIPT_DIR}/compare_pseudo_equivalence.py" \
        --priors-only "${OUTDIR}/priors.A.fa" "${OUTDIR}/priors.B.fa" >/dev/null 2>&1; then
    echo "==> prior sets DIFFER -> arm C: round 2 with arm A's priors"
    rm -rf "${OUTDIR}/C_round2_priorsA"
    "$DADA2RS" dada "${DEREPS[@]}" "${COMMON[@]}" \
        --prior "${OUTDIR}/priors.A.fa" \
        --output-dir "${OUTDIR}/C_round2_priorsA" > "${OUTDIR}/C_round2.log" 2>&1
    ARM_C=("--arm-c" "${OUTDIR}/C_round2_priorsA")
else
    ARM_C=()
fi

# ---------------------------------------------------------------------------
# Compare.
# ---------------------------------------------------------------------------
echo
"${SCRIPT_DIR}/compare_pseudo_equivalence.py" \
    --priors-a "${OUTDIR}/priors.A.fa" \
    --priors-b "${OUTDIR}/priors.B.fa" \
    --pseudo "${OUTDIR}/A_pseudo" \
    --manual "${OUTDIR}/B_round2" \
    --round1 "${OUTDIR}/B_round1" \
    ${ARM_C[@]+"${ARM_C[@]}"} | tee "${OUTDIR}/report.txt"
