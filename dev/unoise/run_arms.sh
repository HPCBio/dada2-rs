#!/usr/bin/env bash
# run_arms.sh — the #323 arm grid on one dataset: DADA2's divisive partition
# against the experimental single pass, with and without the indel-aware λ
# (#291), and the UNOISE birth rule.
#
# Every arm denoises the same reads with the same error model. learn-errors with
# --indel-rates-out learns the substitution model exactly as without it
# (byte-identical JSON) and the indel rates after it, so one learn serves the
# arms with and without indels. Only the dada step differs; merge, table and
# chimera removal are the same for every arm.
#
# Arms (ARMS, space-separated; default all eight):
#   divisive            dada as released
#   sp                  --partition single-pass (p-value birth rule)
#   sp_reassign         ... --reassign (one DADA2 shuffle after the pass)
#   unoise              ... --birth-rule unoise (no indel variant: d counts
#                       gaps by definition)
#   <arm>_indel         the first three with --indel-rates from learn-errors
#   vsearch             UNOISE3 itself, as vsearch runs it, as an external
#                       control: merge (paired), --fastq_maxee filter, derep,
#                       --cluster_unoise, --uchime3_denovo, then every read
#                       mapped to the ZOTUs at 97% identity for the table.
#                       Its chimera removal and table differ from the other
#                       arms by design; it uses no error model.
#
# Usage:
#   run_arms.sh <binary> <out-dir> <manifest.tsv>
#
# The manifest is tab-separated, no header: `sample  fwd.fastq.gz  [rev.fastq.gz]`.
# Two columns is single-end (PacBio); three is paired (Illumina), denoised per
# direction and merged. Inputs must already be primer-trimmed and filtered. Row
# order is training order: learn-errors stops at --nbases, so the first rows
# select the training set.
#
# Environment:
#   THREADS           default 8
#   MODE              pooled (default; dada-pooled) or per-sample (dada)
#   ARMS              see above
#   MIN_ABUND         --min-abund for EVERY arm, including unoise. Unset: the
#                     DADA2 arms keep their default (1) and unoise uses
#                     UNOISE_MIN_ABUND. Set it (e.g. 8) for a matched comparison.
#   UNOISE_MIN_ABUND  unoise's --min-abund when MIN_ABUND is unset; default 8,
#                     vsearch --cluster_unoise's --minsize default
#   EXTRA             flags for learn-errors and dada, e.g. PacBio's
#                     "--band 32 --homo-gap-p -1 --kmer-size 7"
#   LEARN_EXTRA       flags for learn-errors only, e.g. "--errfun pacbio" or
#                     "--errfun binned-qual --binned-quals 2,11,25,37"
#   DADA_EXTRA        flags for every arm's dada step
#   INDEL_EXTRA       flags for the _indel arms only, e.g. "--indel-end-zone 0"
#   MERGE_EXTRA       flags for merge-pairs, e.g. ITS2's "--min-overlap 20"
#   CHIMERA_EXTRA     flags for remove-bimera-denovo (consensus) in every
#                     DADA2-side arm, e.g. "--min-fold-parent-over-abundance 3.5"
#                     for PacBio; the vsearch arm uses UCHIME3 regardless
#   LEARN_DIR         where the error models and indel rates live; default
#                     <out-dir>/learn. Learning is skipped when every file is
#                     already there, so a dev/indels/run_indel_ab.sh learn/
#                     directory can be reused as is.
#   REFERENCE         truth FASTA (mock alleles); scores each arm with
#                     reference-eval
#   REF_EXTRA         flags for reference-eval, e.g. "--max-diffs 0"
#   FORCE             1 re-runs arms that already finished
#   VSEARCH           vsearch binary for the vsearch arm; default `vsearch`
#   VSEARCH_MAXEE     the vsearch arm's --fastq_maxee before denoising;
#                     default 1, the UNOISE3 workflow's
#   VSEARCH_ALPHA     --unoise_alpha; default 2
#   VSEARCH_MERGE_EXTRA  flags for --fastq_mergepairs, e.g.
#                     "--fastq_allowmergestagger" for ITS2
#   The vsearch arm's --minsize follows the unoise arm's --min-abund.
#
# Writes <out-dir>/learn/, then per arm <out-dir>/<arm>/: per-sample JSON per
# direction, <dir>.metrics.json (phase timings), time_<dir>.txt (wall time and
# peak memory), seqtab.json, nochim.json, nochim.fa and, with REFERENCE,
# refeval.json plus per-ASV and per-reference TSVs. Then:
#   python3 dev/unoise/summarize_arms.py <out-dir>
set -eEuo pipefail
# Most steps send stderr to a log, so name the step that failed.
trap 'echo "run_arms.sh: failed (line $LINENO): $BASH_COMMAND; see the logs under $OUT" >&2' ERR

BIN="${1:?usage: run_arms.sh <binary> <out-dir> <manifest.tsv>}"
OUT="${2:?missing out-dir}"
MAN="${3:?missing manifest.tsv}"
THREADS="${THREADS:-8}"
MODE="${MODE:-pooled}"
ARMS="${ARMS:-divisive divisive_indel sp sp_indel sp_reassign sp_reassign_indel unoise vsearch}"
MIN_ABUND="${MIN_ABUND:-}"
UNOISE_MIN_ABUND="${UNOISE_MIN_ABUND:-8}"
LEARN_DIR="${LEARN_DIR:-$OUT/learn}"
REFERENCE="${REFERENCE:-}"
read -r -a extra <<< "${EXTRA:-}"
read -r -a learn_extra <<< "${LEARN_EXTRA:-}"
read -r -a dada_extra <<< "${DADA_EXTRA:-}"
read -r -a indel_extra <<< "${INDEL_EXTRA:-}"
read -r -a merge_extra <<< "${MERGE_EXTRA:-}"
read -r -a chimera_extra <<< "${CHIMERA_EXTRA:-}"
read -r -a ref_extra <<< "${REF_EXTRA:-}"
VSEARCH="${VSEARCH:-vsearch}"
VSEARCH_MAXEE="${VSEARCH_MAXEE:-1}"
VSEARCH_ALPHA="${VSEARCH_ALPHA:-2}"
read -r -a vs_merge_extra <<< "${VSEARCH_MERGE_EXTRA:-}"

case "$MODE" in
  pooled) dada_cmd=dada-pooled ;;
  per-sample) dada_cmd=dada ;;
  *) echo "MODE must be pooled or per-sample, not '$MODE'" >&2; exit 2 ;;
esac
[ -z "$REFERENCE" ] || [ -r "$REFERENCE" ] || { echo "REFERENCE not readable: $REFERENCE" >&2; exit 2; }

samples=(); fwd=(); rev=()
while IFS=$'\t' read -r s f r; do
  [ -z "$s" ] && continue
  samples+=("$s"); fwd+=("$f"); [ -n "${r:-}" ] && rev+=("$r")
done < "$MAN"
n=${#samples[@]}
[ "$n" -gt 0 ] || { echo "no samples in $MAN" >&2; exit 2; }
paired=0
if [ "${#rev[@]}" -gt 0 ]; then
  [ "${#rev[@]}" -eq "$n" ] || { echo "$MAN: some rows lack a reverse file" >&2; exit 2; }
  paired=1
fi
dirs=(fwd); [ "$paired" = 1 ] && dirs=(fwd rev)
names=$(IFS=,; echo "${samples[*]}")
# (Plain loops and copies throughout: macOS ships bash 3.2, without mapfile.)
set_inputs() { if [ "$1" = fwd ]; then in=("${fwd[@]}"); else in=("${rev[@]}"); fi; }

# The per-sample JSON a run writes for input $1 of direction $2. dada-pooled
# names it by --sample-names; dada by the input's stem, as fastq_stem() does.
out_json() {
  local i=$1 d=$2 o=$3 f
  if [ "$MODE" = pooled ]; then echo "$o/$d/${samples[$i]}.json"; return; fi
  if [ "$d" = fwd ]; then f=${fwd[$i]}; else f=${rev[$i]}; fi
  f=$(basename "$f")
  for suf in .fastq.gz .fq.gz .fastq .fq; do f=${f%"$suf"}; done
  echo "$o/$d/$f.json"
}

# Wall time and peak memory: GNU time on Linux, BSD time on macOS; untimed
# (with a note in the log) where /usr/bin/time is missing.
timed() {
  local log=$1; shift
  if /usr/bin/time -v true 2> /dev/null; then /usr/bin/time -v -o "$log" "$@"
  elif /usr/bin/time -l true 2> /dev/null; then /usr/bin/time -l "$@" 2> >(tee "$log" >&2)
  else echo "no /usr/bin/time: untimed" > "$log"; "$@"; fi
}

# The flags that define an arm, and what its output JSON must record for it.
arm_flags() {
  local arm=$1 d=$2 base=${1%_indel}
  flags=(); want_partition=None; want_rule=None; want_reassign=False; want_indel=False
  case "$base" in
    divisive) ;;
    sp) flags=(--partition single-pass); want_partition=single-pass; want_rule=pvalue ;;
    sp_reassign)
      flags=(--partition single-pass --reassign)
      want_partition=single-pass; want_rule=pvalue; want_reassign=True ;;
    unoise)
      flags=(--partition single-pass --birth-rule unoise)
      want_partition=single-pass; want_rule=unoise ;;
    vsearch) ;;
    *) echo "unknown arm '$arm'" >&2; exit 2 ;;
  esac
  if [ "$arm" != "$base" ]; then
    case "$base" in unoise|vsearch)
      echo "$base has no _indel arm: d counts gaps already" >&2; exit 2 ;;
    esac
    flags+=(--indel-rates "$LEARN_DIR/$d.indel_rates.tsv" ${indel_extra[@]+"${indel_extra[@]}"})
    want_indel=True
  fi
  want_min_abund=1
  if [ -n "$MIN_ABUND" ]; then want_min_abund=$MIN_ABUND
  elif [ "$base" = unoise ] || [ "$base" = vsearch ]; then want_min_abund=$UNOISE_MIN_ABUND; fi
  flags+=(--min-abund "$want_min_abund")
}

for a in $ARMS; do arm_flags "$a" fwd; done   # reject a bad ARMS before any work
case " $ARMS " in *" vsearch "*)
  command -v "$VSEARCH" > /dev/null || {
    echo "vsearch arm: '$VSEARCH' not found; set VSEARCH or drop vsearch from ARMS" >&2; exit 2; }
esac

mkdir -p "$OUT" "$LEARN_DIR"
cp "$MAN" "$OUT/manifest.tsv"
{
  echo "binary: $BIN ($("$BIN" --version))"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(hostname)"
  case " $ARMS " in *" vsearch "*) echo "vsearch: $("$VSEARCH" --version 2>&1 | sed -n 1p)" ;; esac
  for v in THREADS MODE ARMS MIN_ABUND UNOISE_MIN_ABUND EXTRA LEARN_EXTRA DADA_EXTRA \
           INDEL_EXTRA MERGE_EXTRA CHIMERA_EXTRA LEARN_DIR REFERENCE REF_EXTRA VSEARCH_MAXEE \
           VSEARCH_ALPHA VSEARCH_MERGE_EXTRA; do
    echo "$v: ${!v:-}"
  done
} > "$OUT/config.txt"

for d in "${dirs[@]}"; do
  if [ -s "$LEARN_DIR/$d.err.json" ] && [ -s "$LEARN_DIR/$d.indel_rates.tsv" ]; then
    echo "==> learn-errors $d: reusing $LEARN_DIR"
    continue
  fi
  echo "==> learn-errors $d ($n input(s))"
  set_inputs "$d"
  timed "$LEARN_DIR/time_$d.txt" "$BIN" learn-errors "${in[@]}" --threads "$THREADS" \
    ${extra[@]+"${extra[@]}"} ${learn_extra[@]+"${learn_extra[@]}"} \
    -o "$LEARN_DIR/$d.err.json" --indel-rates-out "$LEARN_DIR/$d.indel_rates.tsv" \
    --indel-history "$LEARN_DIR/$d.indel_history.tsv" --verbose 2> "$LEARN_DIR/$d.log"
  grep "indel iter" "$LEARN_DIR/$d.log" | tail -1 || true
done

# A dada2-rs arm: dada per direction, merge, table, chimera removal.
run_dada_arm() {
  local arm=$1 o=$2
  for d in "${dirs[@]}"; do
    arm_flags "$arm" "$d"
    mkdir -p "$o/$d"
    echo "==> $arm $dada_cmd $d"
    set_inputs "$d"
    names_arg=(); [ "$MODE" = pooled ] && names_arg=(--sample-names "$names")
    timed "$o/time_$d.txt" "$BIN" "$dada_cmd" "${in[@]}" ${names_arg[@]+"${names_arg[@]}"} \
      --error-model "$LEARN_DIR/$d.err.json" --output-dir "$o/$d" --threads "$THREADS" \
      --metrics-json "$o/$d.metrics.json" ${extra[@]+"${extra[@]}"} \
      ${dada_extra[@]+"${dada_extra[@]}"} "${flags[@]}" --verbose 2> "$o/$d.log"
    # Refuse an arm whose output does not record the settings it was asked for:
    # arms that silently ran alike would look exactly like a null result.
    first=$(out_json 0 "$d" "$o")
    got=$(python3 -I - "$first" <<'PY'
import json, sys
j = json.load(open(sys.argv[1]))
p = j.get("params") or {}
print(p.get("partition"), p.get("birth_rule"), bool(p.get("reassign")),
      "indel_lambda" in json.dumps(j), p.get("min_abund"))
PY
)
    want="$want_partition $want_rule $want_reassign $want_indel $want_min_abund"
    [ "$got" = "$want" ] || { echo "ERROR: $first records '$got', expected '$want'" >&2; exit 1; }
  done
  fj=(); rj=()
  for ((i = 0; i < n; i++)); do
    fj+=("$(out_json "$i" fwd "$o")")
    [ "$paired" = 1 ] && rj+=("$(out_json "$i" rev "$o")")
  done
  if [ "$paired" = 1 ]; then
    echo "==> $arm merge-pairs"
    "$BIN" merge-pairs --fwd-dada "${fj[@]}" --rev-dada "${rj[@]}" --fwd-fastq "${fwd[@]}" \
      --rev-fastq "${rev[@]}" --sample-names "$names" --threads "$THREADS" \
      ${merge_extra[@]+"${merge_extra[@]}"} -o "$o/merged.json" 2> "$o/merge.log"
    "$BIN" make-sequence-table "$o/merged.json" -o "$o/seqtab.json" 2> /dev/null
  else
    "$BIN" make-sequence-table "${fj[@]}" -o "$o/seqtab.json" 2> /dev/null
  fi
  echo "==> $arm remove-bimera-denovo"
  "$BIN" remove-bimera-denovo "$o/seqtab.json" --method consensus --threads "$THREADS" \
    ${chimera_extra[@]+"${chimera_extra[@]}"} -o "$o/nochim.json" 2> "$o/chimera.log"
  "$BIN" seq-table-to-fasta "$o/nochim.json" -o "$o/nochim.fa"
  "$BIN" seq-table-to-fasta "$o/seqtab.json" -o "$o/seqtab.fa"
}

# The vsearch arm: the UNOISE3 workflow. Reads are quality-filtered for
# denoising, but every read is mapped back for the table, as Edgar recommends.
run_vsearch_arm() {
  local o=$1 vd=$1/vsearch i s
  arm_flags vsearch fwd
  mkdir -p "$vd/reads"
  : > "$vd/all.fq"
  for ((i = 0; i < n; i++)); do
    s=${samples[$i]}
    if [ "$paired" = 1 ]; then
      "$VSEARCH" --fastq_mergepairs "${fwd[$i]}" --reverse "${rev[$i]}" \
        --fastqout "$vd/reads/$s.fq" --sample "$s" --relabel "$s." \
        --fastq_qmax 93 --fastq_qmaxout 93 --threads "$THREADS" \
        ${vs_merge_extra[@]+"${vs_merge_extra[@]}"} 2>> "$vd/merge.log"
    else
      "$VSEARCH" --fastq_filter "${fwd[$i]}" --fastqout "$vd/reads/$s.fq" \
        --sample "$s" --relabel "$s." --fastq_qmax 93 2>> "$vd/relabel.log"
    fi
    cat "$vd/reads/$s.fq" >> "$vd/all.fq"
  done
  rm -rf "$vd/reads"
  "$VSEARCH" --fastq_filter "$vd/all.fq" --fastq_qmax 93 --fastaout "$vd/all.fa" \
    --fasta_width 0 2> "$vd/all.log"
  "$VSEARCH" --fastq_filter "$vd/all.fq" --fastq_qmax 93 --fastq_maxee "$VSEARCH_MAXEE" \
    --fastaout "$vd/filt.fa" --fasta_width 0 2> "$vd/filter.log"
  rm -f "$vd/all.fq"
  "$VSEARCH" --derep_fulllength "$vd/filt.fa" --sizeout --relabel Uniq \
    --output "$vd/uniques.fa" --fasta_width 0 2> "$vd/derep.log"
  echo "==> vsearch --cluster_unoise (--minsize $want_min_abund, alpha $VSEARCH_ALPHA)"
  timed "$o/time_fwd.txt" "$VSEARCH" --cluster_unoise "$vd/uniques.fa" \
    --minsize "$want_min_abund" --unoise_alpha "$VSEARCH_ALPHA" --sizeout \
    --centroids "$vd/centroids.fa" --fasta_width 0 --threads "$THREADS" 2> "$vd/unoise.log"
  "$VSEARCH" --uchime3_denovo "$vd/centroids.fa" --nonchimeras "$vd/zotus.fa" \
    --fasta_width 0 2> "$vd/uchime3.log"
  echo "==> vsearch --usearch_global (every read against the ZOTUs, 97%)"
  "$VSEARCH" --usearch_global "$vd/all.fa" --db "$vd/zotus.fa" --id 0.97 --strand plus \
    --otutabout "$vd/otutab.tsv" --threads "$THREADS" 2> "$vd/map.log"
  python3 -I "$(dirname "$0")/vsearch_to_seqtab.py" --centroids "$vd/centroids.fa" \
    --zotus "$vd/zotus.fa" --otutab "$vd/otutab.tsv" --samples "$names" --out-dir "$o"
  rm -f "$vd/all.fa"
}

for arm in $ARMS; do
  o="$OUT/$arm"
  if [ -e "$o/.done" ] && [ "${FORCE:-0}" != 1 ]; then
    echo "==> $arm: done already (FORCE=1 re-runs)"
    continue
  fi
  rm -rf "$o"; mkdir -p "$o"
  if [ "$arm" = vsearch ]; then run_vsearch_arm "$o"; else run_dada_arm "$arm" "$o"; fi
  if [ -n "$REFERENCE" ]; then
    echo "==> $arm reference-eval"
    "$BIN" reference-eval --asvs "$o/nochim.fa" --reference "$REFERENCE" \
      --threads "$THREADS" ${ref_extra[@]+"${ref_extra[@]}"} -o "$o/refeval.json" \
      --per-asv "$o/refeval_asv.tsv" --per-ref "$o/refeval_ref.tsv"
    # The same over the pre-chimera table, with the chimera 2x2.
    "$BIN" reference-eval --asvs "$o/seqtab.fa" --reference "$REFERENCE" \
      --non-chimeric "$o/nochim.fa" --threads "$THREADS" ${ref_extra[@]+"${ref_extra[@]}"} \
      -o "$o/refeval_prechim.json"
  fi
  touch "$o/.done"
done
echo "done: $OUT"
echo "next: python3 dev/unoise/summarize_arms.py $OUT"
