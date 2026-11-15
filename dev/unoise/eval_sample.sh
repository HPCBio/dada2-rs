#!/usr/bin/env bash
# eval_sample.sh — score one sample of a finished run_arms.sh grid (#323)
# against a truth set: typically the mock sample of a run whose other samples
# have no truth, such as MiSeq SOP's `Mock`.
#
# Usage:
#   eval_sample.sh <binary> <out-dir> <sample> <truth.fasta>
#
# For each finished arm (<arm>/.done), the sample's ASVs (count > 0 in its
# column) are written from nochim.json and seqtab.json to
# <arm>/sample_<sample>/, as FASTA with ;size= read counts, and scored as
# run_arms.sh scores the whole table:
#   refeval.json, refeval_asv.tsv, refeval_ref.tsv   the chimera-free ASVs
#   refeval_prechim.json                             pre-chimera, with the
#                                                    chimera x class 2x2
# The vsearch arm has no per-sample table before chimera removal, so it gets
# only the first.
#
# The truth set must be (mock, region, primers)-specific: trimmed to the same
# primer-free region as the reads. Then:
#   python3 dev/unoise/summarize_arms.py <out-dir> --sample <sample>
#
# Environment: THREADS (default 8), REF_EXTRA (flags for reference-eval).
set -eEuo pipefail
trap 'echo "eval_sample.sh: failed (line $LINENO): $BASH_COMMAND" >&2' ERR

BIN="${1:?usage: eval_sample.sh <binary> <out-dir> <sample> <truth.fasta>}"
OUT="${2:?missing out-dir}"
SAMPLE="${3:?missing sample}"
TRUTH="${4:?missing truth.fasta}"
THREADS="${THREADS:-8}"
read -r -a ref_extra <<< "${REF_EXTRA:-}"
[ -r "$TRUTH" ] || { echo "truth not readable: $TRUTH" >&2; exit 2; }

# Writes <table>'s ASVs present in the sample as FASTA; exit 3 when the table
# has no column for the sample, 1 on any other problem.
sample_fasta() {
  python3 - "$1" "$SAMPLE" "$2" <<'PY'
import json, sys
table, sample, out = sys.argv[1:]
with open(table) as fh:
    j = json.load(fh)
if sample not in j["samples"]:
    sys.exit(3)
row = j["counts"][j["samples"].index(sample)]
recs = [(s.upper(), c) for s, c in zip(j["sequences"], row) if c > 0]
recs.sort(key=lambda r: -r[1])
with open(out, "w") as fh:
    for i, (s, c) in enumerate(recs, 1):
        fh.write(f">ASV{i};size={c}\n{s}\n")
print(f"{len(recs)} ASVs, {sum(c for _, c in recs)} reads", file=sys.stderr)
PY
}

found=0
for o in "$OUT"/*/; do
  o=${o%/}; arm=$(basename "$o")
  [ -e "$o/.done" ] || continue
  d="$o/sample_$SAMPLE"
  mkdir -p "$d"
  rc=0; sample_fasta "$o/nochim.json" "$d/nochim.fa" 2> "$d/extract.log" || rc=$?
  if [ "$rc" = 3 ]; then echo "$arm: no sample '$SAMPLE' in nochim.json" >&2; exit 2; fi
  [ "$rc" = 0 ] || { cat "$d/extract.log" >&2; exit 1; }
  found=1
  echo "==> $arm: $SAMPLE $(cat "$d/extract.log")"
  "$BIN" reference-eval --asvs "$d/nochim.fa" --reference "$TRUTH" --threads "$THREADS" \
    ${ref_extra[@]+"${ref_extra[@]}"} -o "$d/refeval.json" \
    --per-asv "$d/refeval_asv.tsv" --per-ref "$d/refeval_ref.tsv" 2> "$d/refeval.log"
  rc=0; sample_fasta "$o/seqtab.json" "$d/seqtab.fa" 2>> "$d/extract.log" || rc=$?
  if [ "$rc" = 0 ]; then
    "$BIN" reference-eval --asvs "$d/seqtab.fa" --reference "$TRUTH" \
      --non-chimeric "$d/nochim.fa" --threads "$THREADS" ${ref_extra[@]+"${ref_extra[@]}"} \
      -o "$d/refeval_prechim.json" 2>> "$d/refeval.log"
  elif [ "$rc" = 3 ]; then
    rm -f "$d/refeval_prechim.json"   # vsearch: no per-sample pre-chimera table
  else
    cat "$d/extract.log" >&2; exit 1
  fi
done
[ "$found" = 1 ] || { echo "$OUT: no finished arms" >&2; exit 2; }
