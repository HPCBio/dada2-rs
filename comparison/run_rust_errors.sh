#!/usr/bin/env bash
# run_rust_errors.sh
# Filter, dereplicate, and learn an error model from a single forward FASTQ
# using dada2-rs, with parameters matching the DADA2 MiSeq SOP tutorial
# (http://benjjneb.github.io/dada2/tutorial.html) and R DADA2 defaults.
#
# Usage:
#   bash run_rust_errors.sh [fastq.gz] [outdir]
#
# Defaults to F3D0 R1 and ./comparison_out/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DADA2RS="${SCRIPT_DIR}/../target/release/dada2-rs"

FASTQ="${1:-/Users/cjfields/projects/hpcbio/tada/debugging/fastq_quals/test/miseq/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq.gz}"
OUTDIR="${2:-${SCRIPT_DIR}/comparison_out}"

if [[ ! -f "$DADA2RS" ]]; then
    echo "ERROR: binary not found at $DADA2RS" >&2
    echo "       Run 'cargo build --release' first." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

SAMPLE="$(basename "$FASTQ" .fastq.gz)"
FILTERED="${OUTDIR}/${SAMPLE}_filtered.fastq.gz"
DEREP_JSON="${OUTDIR}/${SAMPLE}_derep.json"
ERRORS_JSON="${OUTDIR}/${SAMPLE}_errors_rust.json"

# ---------------------------------------------------------------------------
# Step 1: Filter and trim (MiSeq SOP forward-read parameters)
#
#   truncLen=240   — truncate to 240 bp (discard if shorter after truncQ trim)
#   maxN=0         — discard any read containing an N
#   maxEE=2        — discard reads with expected errors > 2
#   truncQ=2       — truncate at first quality score <= 2
#   minLen=20      — Rust default; keeps reads that pass truncQ before truncLen
#
# Note: rm.phix is skipped here (requires a phiX FASTA). The R comparison
# script also sets rm.phix=FALSE so both pipelines see identical reads.
# ---------------------------------------------------------------------------
echo "==> Filter and trim: $FASTQ"
"$DADA2RS" filter-and-trim \
    --fwd    "$FASTQ" \
    --filt   "$FILTERED" \
    --trunc-len 240 \
    --max-n  0 \
    --max-ee 2 \
    --trunc-q 2 \
    --compress \
    --verbose
echo "    -> $FILTERED"

# ---------------------------------------------------------------------------
# Step 2: Dereplicate the filtered reads
# ---------------------------------------------------------------------------
echo "==> Dereplicating: $FILTERED"
"$DADA2RS" derep "$FILTERED" -o "$DEREP_JSON"
echo "    -> $DEREP_JSON"

# ---------------------------------------------------------------------------
# Step 3: Learn the error model
#
# Parameters set to match R DADA2 defaults exactly:
#   OMEGA_A=1e-40, OMEGA_C=1e-40, OMEGA_P=1e-4
#   MIN_FOLD=1, MIN_HAMMING=1, MIN_ABUNDANCE=1
#   MAX_CONSIST=10, GREEDY=TRUE, USE_QUALS=TRUE
#   errorEstimationFunction = loessErrfun
# ---------------------------------------------------------------------------
echo "==> Learning error model (Rust)..."
"$DADA2RS" errors-from-sample "$DEREP_JSON" \
    --omega-a    1e-40 \
    --omega-c    1e-40 \
    --omega-p    1e-4  \
    --min-fold   1.0   \
    --min-hamming 1    \
    --min-abund  1     \
    --max-consist 10   \
    --errfun     loess \
    -o "$ERRORS_JSON"  \
    --verbose
echo "    -> $ERRORS_JSON"

echo ""
echo "Done. Now run the R comparison:"
echo "  Rscript comparison/compare_errors.R $FILTERED $ERRORS_JSON"
