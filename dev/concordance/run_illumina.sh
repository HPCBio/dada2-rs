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
# POOL=pseudo switches the two denoising steps from per-sample `dada` to
# `dada-pseudo` (R `dada(pool="pseudo")`). Unset/false = per-sample, the original
# behaviour. The reference CSV must have been generated with the SAME pool mode --
# write_reference.R takes a matching argument.
POOL="${POOL:-false}"

ALIGN_BACKEND="${ALIGN_BACKEND:-}"
backend_arg=()
[ -n "$ALIGN_BACKEND" ] && backend_arg=(--align-backend "$ALIGN_BACKEND")

# Optional pre-alignment screen backend (kmer|minimizer), experimental. Applied
# ONLY to the screening subcommands (learn-errors, dada, dada-pseudo) -- unlike
# ALIGN_BACKEND, remove-bimera-denovo has no screen and rejects the flag. Unset =
# default (kmer), leaving existing behaviour unchanged.
SCREEN_BACKEND="${SCREEN_BACKEND:-}"
screen_arg=()
if [ -n "$SCREEN_BACKEND" ]; then
  screen_arg=(--screen-backend "$SCREEN_BACKEND")
  [ -n "${MINIMIZER_K:-}" ] && screen_arg+=(--minimizer-k "$MINIMIZER_K")
  [ -n "${MINIMIZER_W:-}" ] && screen_arg+=(--minimizer-w "$MINIMIZER_W")
  # The minimizer metric induces a different distance DISTRIBUTION than the
  # frequency vector, so KDIST_CUTOFF does not transfer between backends even
  # though the two formulas have the same shape -- at 0.42 the sketch passes
  # ~9% of pairs where the frequency vector passes ~28%. Applied to BOTH
  # learn-errors and dada, since both stages screen.
  [ -n "${SCREEN_CUTOFF:-}" ] && screen_arg+=(--kdist-cutoff "$SCREEN_CUTOFF")
fi

# --- Parameters (keep in sync with write_reference.R) ---
TRUNC_LEN_F=240
TRUNC_LEN_R=160
MAX_EE=2
TRUNC_Q=2
MAX_N=0
NBASES=20000000   # learn-errors subsampling cap; small data uses all reads anyway

mkdir -p "$OUT"/{filtered,dada_fwd,dada_rev,control_persample}

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

# Optional: reuse an existing error model instead of learning one. Lets the
# denoising-stage screen effect be isolated from the error-model effect, which
# are otherwise confounded -- the screen is active inside build_trans_mat, so
# changing it changes the learned model as well as the denoising.
ERR_DIR="${ERR_DIR:-}"
if [ -n "$ERR_DIR" ]; then
  echo "==> reusing error models from $ERR_DIR (skipping learn-errors)"
  cp "$ERR_DIR/errF.json" "$OUT/errF.json"
  cp "$ERR_DIR/errR.json" "$OUT/errR.json"
else
echo "==> learn-errors (fwd, rev)"
"$BIN" learn-errors "${filtFs[@]}" --nbases "$NBASES" --errfun loess \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} -o "$OUT/errF.json"
"$BIN" learn-errors "${filtRs[@]}" --nbases "$NBASES" --errfun loess \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} -o "$OUT/errR.json"
fi

if [ "$POOL" = "pseudo" ]; then
  # Pseudo-pooling. NOTE: dada-pseudo takes -o for its output DIRECTORY, whereas
  # `dada` uses --output-dir (-o there means a single-sample output FILE).
  echo "==> dada-pseudo (fwd, rev)"
  "$BIN" dada-pseudo "${filtFs[@]}" --error-model "$OUT/errF.json" \
      -o "$OUT/dada_fwd" --priors-out "$OUT/priors_fwd.fasta" \
      --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"}
  "$BIN" dada-pseudo "${filtRs[@]}" --error-model "$OUT/errR.json" \
      -o "$OUT/dada_rev" --priors-out "$OUT/priors_rev.fasta" \
      --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"}

  # Positive control: did the priors actually CHANGE anything?
  #
  # A non-empty prior set is NOT sufficient. The 2-sample dev/concordance/data/
  # illumina fixture selects 8 priors that alter nothing, so round 2 equals round 1,
  # the pseudo path goes untested, and every comparison still passes -- a green
  # result that means nothing. That trap has bitten twice (see
  # docs/findings/pseudo-pooling-priors-vs-error-model.md), so this checks against
  # an actual per-sample run rather than trusting a prior count. Costs one extra
  # dada pass on a small fixture.
  echo "==> positive control: per-sample run for comparison"
  "$BIN" dada "${filtFs[@]}" --error-model "$OUT/errF.json" \
      --output-dir "$OUT/control_persample" --threads "$THREADS" \
      ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"} > /dev/null
  python3 - "$OUT/dada_fwd" "$OUT/control_persample" <<'PYCTL'
import glob, json, os, sys
def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "*.json"))):
        j = json.load(open(p))
        out[j["sample"]] = {a["sequence"] for a in j["asvs"]}
    return out
pseudo, ctrl = load(sys.argv[1]), load(sys.argv[2])
if pseudo == ctrl:
    sys.exit("run_illumina.sh: POOL=pseudo output is IDENTICAL to a per-sample run.\n"
             "  The priors were inert, so the pseudo path is untested and any PASS\n"
             "  below is meaningless. The fixture needs >=2 samples sharing\n"
             "  low-abundance variants.")
gained = sum(len(pseudo[s] - ctrl[s]) for s in ctrl if s in pseudo)
print(f"==> control OK: pseudo rescued {gained} ASV(s) a per-sample run did not call")
PYCTL
else
  # Per-sample denoising (pool=FALSE analog) — matches R dada() default.
  echo "==> dada (fwd, rev; per-sample)"
  "$BIN" dada "${filtFs[@]}" --error-model "$OUT/errF.json" \
      --output-dir "$OUT/dada_fwd" --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"}
  "$BIN" dada "${filtRs[@]}" --error-model "$OUT/errR.json" \
      --output-dir "$OUT/dada_rev" --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} ${screen_arg[@]+"${screen_arg[@]}"}
fi

echo "==> merge-pairs"
# merge-pairs pairs its four lists POSITIONALLY, so the dada JSONs must be
# enumerated in the same order as filtFs/filtRs -- NOT globbed independently.
# Globbing breaks on the MiSeq SOP and any similar naming: the FASTQ glob sorts
# `F3D1F.fastq.gz` against `F3D141F.fastq.gz` while the JSON glob sorts
# `F3D1_F_filt.json` against `F3D141_F_filt.json`, and the inserted `_` flips
# their relative order. The 2-sample fixture never exposed it; at 20 samples it
# silently mis-pairs 10 of them, and only a downstream map-length check turns it
# into a visible error.
dadaFs=(); dadaRs=()
for ff in "${filtFs[@]}"; do
  dadaFs+=("$OUT/dada_fwd/$(basename "$ff" .fastq.gz).json")
done
for fr in "${filtRs[@]}"; do
  dadaRs+=("$OUT/dada_rev/$(basename "$fr" .fastq.gz).json")
done
"$BIN" merge-pairs \
    --fwd-dada "${dadaFs[@]}" --rev-dada "${dadaRs[@]}" \
    --fwd-fastq "${filtFs[@]}" --rev-fastq "${filtRs[@]}" \
    --threads "$THREADS" -o "$OUT/merged.json"

echo "==> make-sequence-table"
"$BIN" make-sequence-table "$OUT/merged.json" -o "$OUT/seqtab.json"

echo "==> remove-bimera-denovo"
"$BIN" remove-bimera-denovo "$OUT/seqtab.json" --method consensus \
    --threads "$THREADS" ${backend_arg[@]+"${backend_arg[@]}"} -o "$OUT/seqtab.nochim.json"

echo "==> done: $OUT/seqtab.nochim.json"
