#!/usr/bin/env bash
# run_nbases_ladder.sh — is the learn-errors model converged, and does it depend
# on WHICH samples it saw?
#
# Two questions the small MiSeq SOP subset cannot answer, because it holds only
# ~33 Mbp total and because `learn-errors` accumulates WHOLE FILES (it never
# subsamples reads within a file, see `--nbases` help and issue #68). On that
# subset a `--nbases` ladder silently varies *how many samples* were read as much
# as *how much data*, so depth and composition are confounded.
#
# This script separates them:
#
#   ARM A (depth)        every sample pre-subsampled to the SAME read depth, then
#                        --nbases laddered. Each rung adds depth uniformly rather
#                        than adding whole samples.
#   ARM B (composition)  --nbases held FIXED, only the sample draw varies
#                        (--randomize with several seeds).
#
# Compare the two arms with dev/compare_error_models.py. If arm B's spread is
# comparable to arm A's, sample choice matters as much as data volume, and
# raising --nbases is not the lever it appears to be.
#
# Usage:
#   run_nbases_ladder.sh <dada2rs-binary> <fastq-dir> <out-dir> [threads]
#
#   <fastq-dir> holds the RAW FASTQs (*.fastq.gz) for ONE read orientation --
#   forward only, say. A mixed R1/R2 set is detected and rejected, because
#   learn-errors would silently fit one model across two different error
#   profiles. They are filtered here with the MiSeq SOP parameters; edit
#   TRUNC_LEN/MAX_EE below for other data.
set -euo pipefail

BIN="${1:?usage: run_nbases_ladder.sh <binary> <fastq-dir> <out-dir> [threads]}"
DATA="${2:?missing fastq-dir}"
OUT="${3:?missing out-dir}"
THREADS="${4:-16}"

# Filtering (MiSeq SOP forward-read defaults).
TRUNC_LEN=240
MAX_EE=2
TRUNC_Q=2
MAX_N=0

# Per-sample read depth for ARM A. Every sample is cut to exactly this many
# reads, so a --nbases rung means "more depth from the same samples" rather than
# "more samples". Set to the smallest per-sample depth in the run, or lower.
EQUAL_DEPTH="${EQUAL_DEPTH:-20000}"

# ARM A ladder. Deliberately extends past 1e8: the KDIST cutoff study found the
# model still moving at 1e8, and nothing has yet tested beyond it.
LADDER="${LADDER:-20000000 100000000 500000000 1000000000 5000000000}"

# ARM B seeds, at a fixed budget.
SEEDS="${SEEDS:-1 2 3 4 5}"
FIXED_NBASES="${FIXED_NBASES:-100000000}"

ERRFUN="${ERRFUN:-loess}"

mkdir -p "$OUT"/{filt,equal,models_depth,models_draw}

echo "==> [1/4] filter-and-trim"
shopt -s nullglob
raws=("$DATA"/*.fastq.gz)
if [ ${#raws[@]} -eq 0 ]; then
  echo "run_nbases_ladder.sh: no *.fastq.gz in $DATA" >&2
  exit 1
fi
# Refuse an obviously mixed forward/reverse set.
#
# learn-errors fits ONE error model over whatever it is given, and forward and
# reverse reads have materially different error profiles -- a model fit across
# both is not slightly wrong, it is uninterpretable. Nothing downstream reports
# this, and on a cluster the cost is hours of runtime for a meaningless model.
#
# Deliberately a detector rather than a hardcoded R1 glob: naming conventions
# vary too much for a fixed pattern to be anything but brittle, so this accepts
# ANY set of files and only rejects one where both orientations are plainly
# present. Set ALLOW_MIXED_READS=1 if the match is a false positive (e.g. a
# sample legitimately named "..._R2_site").
fwd=0; rev=0
for f in "${raws[@]}"; do
  case "$(basename "$f")" in
    *_R1_*|*_R1.*|*_1.fastq*|*_1.fq*) fwd=1 ;;
  esac
  case "$(basename "$f")" in
    *_R2_*|*_R2.*|*_2.fastq*|*_2.fq*) rev=1 ;;
  esac
done
if [ "$fwd" = 1 ] && [ "$rev" = 1 ] && [ "${ALLOW_MIXED_READS:-0}" != 1 ]; then
  cat >&2 <<'MIXED'
run_nbases_ladder.sh: input looks like BOTH forward and reverse reads.

learn-errors fits a single error model over every file it is given, and R1/R2
have different error profiles, so the result would be uninterpretable -- and
nothing downstream would tell you. Point this at one orientation only:

    mkdir r1 && ln -s /path/to/raw/*_R1_*.fastq.gz r1/
    run_nbases_ladder.sh <binary> r1 <outdir> <threads>

Set ALLOW_MIXED_READS=1 to override if this is a false positive.
MIXED
  exit 1
fi

for f in "${raws[@]}"; do
  b=$(basename "$f" .fastq.gz)
  [ -f "$OUT/filt/$b.fastq.gz" ] && continue
  "$BIN" filter-and-trim --fwd "$f" --filt "$OUT/filt/$b.fastq.gz" \
      --trunc-len "$TRUNC_LEN" --max-n "$MAX_N" --max-ee "$MAX_EE" \
      --trunc-q "$TRUNC_Q" --compress > /dev/null
done
echo "    filtered: $(ls "$OUT"/filt/*.fastq.gz | wc -l) samples"

echo "==> [2/4] equalise per-sample depth to $EQUAL_DEPTH reads (ARM A only)"
for f in "$OUT"/filt/*.fastq.gz; do
  b=$(basename "$f")
  [ -f "$OUT/equal/$b" ] && continue
  # Truncate to the first EQUAL_DEPTH reads (4 lines each). Deliberately not
  # `dada2-rs sample`, which dereplicates to JSON rather than writing FASTQ.
  # Taking a prefix is fine for Illumina, where read order carries no quality
  # ordering; a sample shallower than EQUAL_DEPTH passes through whole.
  #
  # pipefail is disabled for this one pipeline on purpose: `head` closes the pipe
  # as soon as it has its lines, `gzip -cd` then dies of SIGPIPE, and pipefail
  # would abort the script on what is the normal, intended outcome.
  set +o pipefail
  gzip -cd "$f" | head -n $(( EQUAL_DEPTH * 4 )) | gzip -c > "$OUT/equal/$b"
  set -o pipefail
done
echo "    equalised: $(ls "$OUT"/equal/*.fastq.gz | wc -l) samples"

echo "==> [3/4] ARM A — depth ladder on equal-depth samples"
for N in $LADDER; do
  o="$OUT/models_depth/err_$N.json"
  [ -f "$o" ] && { echo "    nbases=$N (cached)"; continue; }
  echo "    nbases=$N"
  "$BIN" learn-errors "$OUT"/equal/*.fastq.gz --nbases "$N" --errfun "$ERRFUN" \
      --threads "$THREADS" -o "$o" > /dev/null
done

echo "==> [4/4] ARM B — sample-draw sensitivity at fixed nbases=$FIXED_NBASES"
for S in $SEEDS; do
  o="$OUT/models_draw/err_seed$S.json"
  [ -f "$o" ] && { echo "    seed=$S (cached)"; continue; }
  echo "    seed=$S"
  "$BIN" learn-errors "$OUT"/filt/*.fastq.gz --nbases "$FIXED_NBASES" \
      --errfun "$ERRFUN" --randomize --seed "$S" --threads "$THREADS" -o "$o" > /dev/null
done

echo
echo "================ ARM A: depth (equal-depth samples) ================"
args=()
for N in $LADDER; do args+=("$(( N / 1000000 ))Mb=$OUT/models_depth/err_$N.json"); done
python3 "$(dirname "$0")/compare_error_models.py" "${args[@]}"

echo
echo "================ ARM B: sample draw (fixed nbases) ================="
args=()
for S in $SEEDS; do args+=("seed$S=$OUT/models_draw/err_seed$S.json"); done
python3 "$(dirname "$0")/compare_error_models.py" "${args[@]}"

echo
cat <<'NOTE'
================ how to read this =================================
ARM A tells you whether the model has converged: the step-to-step
"max fold" should fall toward 1.00x. If the LAST step is as large as
the previous one, the model is still moving and the top rung is not a
converged model.

ARM B tells you how much of any movement is really sample choice. If
ARM B's spread is comparable to ARM A's, then raising --nbases is not
buying convergence, it is buying a different sample mix (issue #68).

CHECK `nq` IN BOTH TABLES FIRST. It is max_q + 1 and varies with the
data. Comparing models with different nq without aligning by quality
column manufactures huge spurious differences that land exactly on the
nq boundary. compare_error_models.py aligns and truncates correctly,
but a differing nq is still worth noticing.
NOTE
