#!/usr/bin/env bash
# run_illumina.sh — dada2-rs paired-end Illumina pipeline for the concordance
# guardrail. Produces a chimera-filtered sequence table (seqtab.nochim.json) from
# a small set of paired FASTQs, which compare_to_reference.py then diffs against
# the static R DADA2 reference CSV.
#
# The parameters here MUST match write_reference.R exactly (same truncLen, maxEE,
# truncQ, pool=FALSE), or the comparison is apples-to-oranges. Keep them in sync.
#
# Usage: run_illumina.sh <dada2rs-binary> <data-dir> <out-dir> [threads]
#   <data-dir> holds paired files named <sample>F.fastq.gz / <sample>R.fastq.gz
#              (e.g. sam1F.fastq.gz + sam1R.fastq.gz).
set -euo pipefail

BIN="${1:?usage: run_illumina.sh <binary> <data-dir> <out-dir> [threads]}"
DATA="${2:?missing data-dir}"
OUT="${3:?missing out-dir}"
THREADS="${4:-2}"

# Optional alignment backend (nw|wfa2). When set, it is threaded through every
# alignment-using subcommand (learn-errors, dada, remove-bimera-denovo) so the
# concordance guardrail can run the whole pipeline with WFA. Unset = default
# (nw), leaving existing behavior unchanged. The `+"${...}"` form keeps the
# empty-array expansion safe under `set -u`.
ALIGN_BACKEND="${ALIGN_BACKEND:-}"
backend_arg=()
[ -n "$ALIGN_BACKEND" ] && backend_arg=(--align-backend "$ALIGN_BACKEND")

# --- Parameters (keep in sync with write_reference.R) ---
TRUNC_LEN_F=240
TRUNC_LEN_R=160
MAX_EE=2
TRUNC_Q=2
MAX_N=0
NBASES=20000000   # learn-errors subsampling cap; small data uses all reads anyway

mkdir -p "$OUT"/{filtered,dada_fwd,dada_rev}

fwds=("$DATA"/*F.fastq.gz)
if [ ! -e "${fwds[0]}" ]; then
  echo "run_illumina.sh: no *F.fastq.gz in $DATA" >&2
  exit 1
fi

filtFs=(); filtRs=()
for f in "${fwds[@]}"; do
  base=$(basename "$f")
  name=${base%F.fastq.gz}
  r="$DATA/${name}R.fastq.gz"
  [ -e "$r" ] || { echo "run_illumina.sh: missing reverse mate for $f" >&2; exit 1; }
  ff="$OUT/filtered/${name}_F_filt.fastq.gz"
  fr="$OUT/filtered/${name}_R_filt.fastq.gz"
  echo "==> filter-and-trim $name"
  "$BIN" filter-and-trim --fwd "$f" --filt "$ff" --rev "$r" --filt-rev "$fr" \
      --trunc-len "$TRUNC_LEN_F" "$TRUNC_LEN_R" --max-n "$MAX_N" \
      --max-ee "$MAX_EE" "$MAX_EE" --trunc-q "$TRUNC_Q" --compress
  filtFs+=("$ff"); filtRs+=("$fr")
done

echo "==> learn-errors (fwd, rev)"
"$BIN" learn-errors "${filtFs[@]}" --nbases "$NBASES" --errfun loess \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} -o "$OUT/errF.json"
"$BIN" learn-errors "${filtRs[@]}" --nbases "$NBASES" --errfun loess \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} -o "$OUT/errR.json"

# Per-sample denoising (pool=FALSE analog) — matches R dada() default.
echo "==> dada (fwd, rev; per-sample)"
"$BIN" dada "${filtFs[@]}" --error-model "$OUT/errF.json" \
    --output-dir "$OUT/dada_fwd" --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"}
"$BIN" dada "${filtRs[@]}" --error-model "$OUT/errR.json" \
    --output-dir "$OUT/dada_rev" --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"}

echo "==> merge-pairs"
"$BIN" merge-pairs \
    --fwd-dada "$OUT"/dada_fwd/*.json --rev-dada "$OUT"/dada_rev/*.json \
    --fwd-fastq "${filtFs[@]}" --rev-fastq "${filtRs[@]}" \
    --threads "$THREADS" -o "$OUT/merged.json"

echo "==> make-sequence-table"
"$BIN" make-sequence-table "$OUT/merged.json" -o "$OUT/seqtab.json"

echo "==> remove-bimera-denovo"
"$BIN" remove-bimera-denovo "$OUT/seqtab.json" --method consensus \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} -o "$OUT/seqtab.nochim.json"

echo "==> done: $OUT/seqtab.nochim.json"
