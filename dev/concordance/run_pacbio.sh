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

BIN="${1:?usage: run_pacbio.sh <binary> <data-dir> <out-dir> [primer_fwd primer_rev] [threads]}"
DATA="${2:?missing data-dir}"
OUT="${3:?missing out-dir}"

# PREFILTERED=1: <data-dir> already holds primer-stripped, length/EE-filtered
# reads, so remove-primers and the filter are skipped and the primers are not
# needed. Mirrors run_illumina.sh's PREFILTERED and write_reference.R's
# --prefiltered. This is what lets ONE filtering pass feed both tools, rather
# than each filtering separately and the comparison silently carrying that
# difference too.
PREFILTERED="${PREFILTERED:-}"

PRIMER_FWD="${4:-}"
PRIMER_REV="${5:-}"
if [ -z "$PREFILTERED" ]; then
  : "${PRIMER_FWD:?missing primer_fwd (or set PREFILTERED=1)}"
  : "${PRIMER_REV:?missing primer_rev (or set PREFILTERED=1)}"
fi
# Positional 6 still works; THREADS= is the form to use with PREFILTERED=1,
# where there are no primer arguments to hold its place.
THREADS="${THREADS:-${6:-2}}"

# Optional alignment backend (nw|wfa2), threaded through the alignment-using
# subcommands (learn-errors, dada, remove-bimera-denovo) so the concordance
# guardrail can run the pipeline with WFA. Unset = default (nw). The
# `+"${...}"` form keeps empty-array expansion safe under `set -u`.
ALIGN_BACKEND="${ALIGN_BACKEND:-}"
backend_arg=()
[ -n "$ALIGN_BACKEND" ] && backend_arg=(--align-backend "$ALIGN_BACKEND")

# Optional pre-alignment screen backend (kmer|minimizer), experimental. Applied
# ONLY to the screening subcommands -- unlike ALIGN_BACKEND, remove-bimera-denovo
# has no screen and rejects the flag. Unset = default (kmer).
SCREEN_BACKEND="${SCREEN_BACKEND:-}"
screen_arg=()
if [ -n "$SCREEN_BACKEND" ]; then
  screen_arg=(--screen-backend "$SCREEN_BACKEND")
  [ -n "${MINIMIZER_K:-}" ] && screen_arg+=(--minimizer-k "$MINIMIZER_K")
  [ -n "${MINIMIZER_W:-}" ] && screen_arg+=(--minimizer-w "$MINIMIZER_W")
fi

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
# learnErrors' subsampling budget. MUST be raised past the training set when
# PACBIO_TRAIN_SAMPLES pins one -- otherwise dada2-rs trains on a PREFIX of the
# manifest while R (given a matching --nbases) uses all of it, and the arms
# differ in training data as well as in whatever is under test. write_reference.R
# refuses that on the R side; nothing can refuse it here, because this script
# hands learn-errors an explicit file list and cannot know the manifest's size
# without reading every FASTQ. So it is checked below instead.
NBASES="${NBASES:-200000000}"

# Pin the error-model training set to these samples (one name per line, matching
# the filtered stems without the _filt suffix), then denoise everything --
# the same contract as write_reference.R's --train-samples. Produce the manifest
# with pick_training_subset.py --suffix _filt.fastq.gz so both sides agree.
PACBIO_TRAIN_SAMPLES="${PACBIO_TRAIN_SAMPLES:-}"

# Full pooling (R pool=TRUE) instead of per-sample dada. Must match the R
# reference's --pool, or the comparison is measuring the pooling mode.
POOL="${POOL:-false}"

# VERBOSE=1: per-step progress to stderr. A pooled PacBio run can sit in
# `dada-pooled` for hours with nothing on the terminal, which is
# indistinguishable from a hang; `--verbose` prints the bud loop's progress
# line so you can see it advancing and estimate the remaining work.
# METRICS=1 additionally writes <step>.metrics.json per denoising step --
# machine-readable phase timings, which is what to attach when a run is slow
# rather than a wall-clock number (see docs/findings on the screen vs align
# split).
VERBOSE="${VERBOSE:-}"
verbose_arg=()
[ -n "$VERBOSE" ] && verbose_arg=(--verbose)
METRICS="${METRICS:-}"

mkdir -p "$OUT"/{filtered,dada}

reads=("$DATA"/*.fastq.gz)
if [ ! -e "${reads[0]}" ]; then
  echo "run_pacbio.sh: no *.fastq.gz in $DATA" >&2
  exit 1
fi

# `names` runs parallel to `filts` so the training pin can resolve a manifest
# entry to a file without assuming how the file is named.
filts=()
names=()
if [ -n "$PREFILTERED" ]; then
  echo "==> PREFILTERED: using ${#reads[@]} input(s) as-is (no remove-primers, no filter)"
  for f in "${reads[@]}"; do
    name=$(basename "$f" .fastq.gz)
    # Strip the suffix run_pacbio.sh itself writes, so a manifest built from
    # either raw or filtered names resolves. Matches write_reference.R, which
    # strips `_filt` only under --prefiltered.
    name=${name%_filt}
    filts+=("$f")
    names+=("$name")
  done
else
  for f in "${reads[@]}"; do
    name=$(basename "$f" .fastq.gz)
    ff="$OUT/filtered/${name}_filt.fastq.gz"
    echo "==> remove-primers + filter $name"
    "$BIN" remove-primers "$f" --fout "$ff" \
        --primer-fwd "$PRIMER_FWD" --primer-rev "$PRIMER_REV" \
        --max-mismatch "$MAX_MISMATCH" --trim-fwd --trim-rev --orient \
        --min-len "$MIN_LEN" --max-len "$MAX_LEN" --max-n "$MAX_N" \
        --max-ee "$MAX_EE" --trunc-q "$TRUNC_Q" --compress \
        ${verbose_arg[@]+"${verbose_arg[@]}"} \
        -o "$OUT/primers_${name}.json"
    filts+=("$ff")
    names+=("$name")
  done
fi

trains=("${filts[@]}")
if [ -n "$PACBIO_TRAIN_SAMPLES" ]; then
  trains=()
  while IFS= read -r want; do
    want=$(printf '%s' "$want" | tr -d '[:space:]')   # tolerate CRLF manifests
    [ -z "$want" ] && continue
    found=""
    for i in "${!names[@]}"; do
      if [ "${names[$i]}" = "$want" ]; then found="${filts[$i]}"; break; fi
    done
    if [ -z "$found" ]; then
      echo "run_pacbio.sh: --train-samples name not present in $DATA: $want" >&2
      exit 1
    fi
    trains+=("$found")
  done < "$PACBIO_TRAIN_SAMPLES"
  echo "==> training on ${#trains[@]} of ${#filts[@]} samples (pinned via $PACBIO_TRAIN_SAMPLES)"
  # Guard the mirror image of the trap write_reference.R catches: a pinned
  # manifest bigger than NBASES means learn-errors silently takes a prefix.
  train_bases=$(gzip -cd "${trains[@]}" | awk 'NR%4==2 {n+=length($0)} END {print n+0}')
  echo "==> training set: $train_bases bases; nbases = $NBASES"
  if [ "$train_bases" -gt "$NBASES" ]; then
    echo "run_pacbio.sh: NBASES ($NBASES) is smaller than the pinned training set" >&2
    echo "  ($train_bases bases), so learn-errors would train on only a prefix." >&2
    echo "  Re-run with NBASES=$(( train_bases * 11 / 10 )) or larger." >&2
    exit 1
  fi
fi

# ERR_DIR: reuse an error model instead of learning one, mirroring
# run_illumina.sh. The k-mer screen is active inside build_trans_mat, so
# PACBIO_KMER changes the LEARNED MODEL as well as the denoising -- a k=7 run
# is therefore not comparable to an R reference on the model axis, because R's
# KMER_SIZE is a compile-time 5. This lets the two be decoupled: learn once at
# k=5 to match R, then denoise at k=7 for the ~7x screen saving.
ERR_DIR="${ERR_DIR:-}"
if [ -n "$ERR_DIR" ]; then
  echo "==> reusing error model from $ERR_DIR (skipping learn-errors)"
  cp "$ERR_DIR/err.json" "$OUT/err.json"
else
echo "==> learn-errors (pacbio errfun, k=$KMER)"
"$BIN" learn-errors "${trains[@]}" --nbases "$NBASES" --errfun pacbio \
    --band "$BAND" --kmer-size "$KMER" --threads "$THREADS" \
    ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} \
    ${verbose_arg[@]+"${verbose_arg[@]}"} -o "$OUT/err.json"
fi

if [ "$POOL" = "true" ]; then
  echo "==> dada-pooled (full pooling)"
  metrics_arg=()
  [ -n "$METRICS" ] && metrics_arg=(--metrics-json "$OUT/dada-pooled.metrics.json")
  "$BIN" dada-pooled "${filts[@]}" --error-model "$OUT/err.json" \
      -o "$OUT/dada" --band "$BAND" --kmer-size "$KMER" --threads "$THREADS" \
      ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} \
      ${verbose_arg[@]+"${verbose_arg[@]}"} ${metrics_arg[@]+"${metrics_arg[@]}"}
else
  echo "==> dada (per-sample)"
  metrics_arg=()
  [ -n "$METRICS" ] && metrics_arg=(--metrics-json "$OUT/dada.metrics.json")
  "$BIN" dada "${filts[@]}" --error-model "$OUT/err.json" \
      --output-dir "$OUT/dada" --band "$BAND" --kmer-size "$KMER" --threads "$THREADS" \
      ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} \
      ${verbose_arg[@]+"${verbose_arg[@]}"} ${metrics_arg[@]+"${metrics_arg[@]}"}
fi

echo "==> make-sequence-table"
"$BIN" make-sequence-table "$OUT"/dada/*.json -o "$OUT/seqtab.json"

echo "==> remove-bimera-denovo"
"$BIN" remove-bimera-denovo "$OUT/seqtab.json" --method consensus \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} \
    ${verbose_arg[@]+"${verbose_arg[@]}"} -o "$OUT/seqtab.nochim.json"

echo "==> done: $OUT/seqtab.nochim.json"
