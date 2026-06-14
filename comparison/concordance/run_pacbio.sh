#!/usr/bin/env bash
# run_pacbio.sh — dada2-rs single-end PacBio pipeline for the concordance
# guardrail. Produces seqtab.nochim.json from a small set of raw (primered)
# PacBio FASTQs, which compare_to_reference.py diffs against the static R
# reference CSV.
#
# Parameters MUST match write_reference.R (pacbio branch): same primers, length
# filters, PacBio errfun, BAND_SIZE=32, pool=FALSE. NOTE: --kmer-size 5 is used
# deliberately to match R DADA2's fixed KMER_SIZE=5, so the comparison is
# apples-to-apples (dada2-rs defaults to k=7 for PacBio speed, but that is a
# screening-only difference; the reference is k=5).
#
# Usage: run_pacbio.sh <binary> <data-dir> <out-dir> <primer_fwd> <primer_rev> [threads]
#   <data-dir> holds raw <sample>.fastq.gz (primered, single-end).
set -euo pipefail

BIN="${1:?usage: run_pacbio.sh <binary> <data-dir> <out-dir> <primer_fwd> <primer_rev> [threads]}"
DATA="${2:?missing data-dir}"
OUT="${3:?missing out-dir}"
PRIMER_FWD="${4:?missing primer_fwd}"
PRIMER_REV="${5:?missing primer_rev}"
THREADS="${6:-2}"

# --- Parameters (keep in sync with write_reference.R pacbio branch) ---
MIN_LEN=1000
MAX_LEN=1600
MAX_EE=2
TRUNC_Q=0
MAX_N=0
BAND=32
# Default k=5 matches R's fixed KMER_SIZE (apples-to-apples vs the reference).
# Override with PACBIO_KMER=7 to spot-check dada2-rs's recommended PacBio setting
# against the same R(k=5) reference — the k-mer screen is a prefilter, so this
# should give the same ASVs (see issue #15).
KMER="${PACBIO_KMER:-5}"
MAX_MISMATCH=2
NBASES=200000000

mkdir -p "$OUT"/{filtered,dada}

reads=("$DATA"/*.fastq.gz)
if [ ! -e "${reads[0]}" ]; then
  echo "run_pacbio.sh: no *.fastq.gz in $DATA" >&2
  exit 1
fi

filts=()
for f in "${reads[@]}"; do
  name=$(basename "$f" .fastq.gz)
  ff="$OUT/filtered/${name}_filt.fastq.gz"
  echo "==> remove-primers + filter $name"
  "$BIN" remove-primers "$f" --fout "$ff" \
      --primer-fwd "$PRIMER_FWD" --primer-rev "$PRIMER_REV" \
      --max-mismatch "$MAX_MISMATCH" --trim-fwd --trim-rev --orient \
      --min-len "$MIN_LEN" --max-len "$MAX_LEN" --max-n "$MAX_N" \
      --max-ee "$MAX_EE" --trunc-q "$TRUNC_Q" --compress \
      -o "$OUT/primers_${name}.json"
  filts+=("$ff")
done

echo "==> learn-errors (pacbio errfun, k=$KMER)"
"$BIN" learn-errors "${filts[@]}" --nbases "$NBASES" --errfun pacbio \
    --band "$BAND" --kmer-size "$KMER" --threads "$THREADS" -o "$OUT/err.json"

echo "==> dada (per-sample)"
"$BIN" dada "${filts[@]}" --error-model "$OUT/err.json" \
    --output-dir "$OUT/dada" --band "$BAND" --kmer-size "$KMER" --threads "$THREADS"

echo "==> make-sequence-table"
"$BIN" make-sequence-table "$OUT"/dada/*.json -o "$OUT/seqtab.json"

echo "==> remove-bimera-denovo"
"$BIN" remove-bimera-denovo "$OUT/seqtab.json" --method consensus \
    --threads "$THREADS" -o "$OUT/seqtab.nochim.json"

echo "==> done: $OUT/seqtab.nochim.json"
