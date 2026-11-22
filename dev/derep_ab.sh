#!/usr/bin/env bash
# derep_ab.sh — A/B: pooled denoising from FASTQ vs from pre-dereplicated JSON.
# ---------------------------------------------------------------------------
# Question (#1): does a dedicated per-sample `derep` step, run before
# `dada-pooled`, shorten the pooled run's serial load front, and does it leave
# the denoised output unchanged?
#
# Arm A  dada-pooled on the FASTQs. The serial load front dereplicates every
#        sample in turn (`derep_detail.build`).
# Step   `derep` each FASTQ to compact `.json.gz`, DEREP_JOBS at a time. This is
#        the work a workflow would move into per-sample jobs.
# Arm B  dada-pooled on those JSONs. The load front now reads and parses JSON
#        (`derep_detail.read` / `.parse`).
#
# Both arms get the same inputs in the same order and the same sample names,
# because the pooled tie-break is input-order sensitive (#260). The per-sample
# outputs must match exactly apart from `input_file`, which names the arm's own
# input; `source_fastq` must still match.
#
# Usage: derep_ab.sh <dada2rs-binary> <error-model.json> <out-dir> <fastq>...
#
# Environment:
#   THREADS     dada-pooled --threads (default 8)
#   DEREP_JOBS  concurrent derep processes (default 1; 1 gives the serial cost)
#   DADA_ARGS   extra dada-pooled flags, whitespace-split (PacBio: band size,
#               homopolymer gap, etc. -- whatever the production run uses)
#   PREFIX      command prefix for every timed step, whitespace-split, e.g.
#               PREFIX="$NUMA_PREFIX" after sourcing dev/numa_pin.sh on the
#               cluster (see feedback on NUMA placement in that file)
#   STRIP       shell glob (not a regex) removed from the end of each FASTQ
#               basename to form the sample name, shortest match. It is tried on
#               the full name, then on the name minus .fastq.gz/.fq.gz/.fastq/.fq,
#               so STRIP='_R[12].fastq.gz', STRIP='_R[12]' and STRIP='_R[12]_*'
#               all work; quote it so the shell does not expand it. Default none
#               (the extension is still removed).
#               Paired data: pass one direction per run (R1s only, say).
#   SKIP_A=1    reuse an existing arm A (out-dir/a_fastq) instead of rerunning
#
# Writes out-dir/summary.tsv and prints it. Exits non-zero if the arms differ.
set -euo pipefail

BIN="${1:?usage: derep_ab.sh <binary> <error-model> <out-dir> <fastq>...}"
ERR="${2:?missing error model}"
OUT="${3:?missing out-dir}"
shift 3
[ $# -gt 0 ] || { echo "derep_ab.sh: no FASTQ inputs" >&2; exit 1; }

THREADS="${THREADS:-8}"
DEREP_JOBS="${DEREP_JOBS:-1}"
STRIP="${STRIP:-}"
# shellcheck disable=SC2206
dada_extra=(${DADA_ARGS:-})
# shellcheck disable=SC2206
prefix=(${PREFIX:-})

now() { python3 -c 'import time; print(f"{time.time():.3f}")'; }

mkdir -p "$OUT"/{a_fastq,b_json,derep,logs}

names=()
for f in "$@"; do
  n=$(basename "$f")
  # STRIP is a glob, tried first on the full basename, then (only if that did
  # not match) on the extension-less name. One strip at most, so `_R[12]*`
  # cannot also eat an `_R1` that is part of the sample name.
  if [ -n "$STRIP" ] && [ "${n%$STRIP}" != "$n" ]; then
    n=${n%$STRIP}
    for ext in .fastq.gz .fq.gz .fastq .fq; do n=${n%"$ext"}; done
  else
    for ext in .fastq.gz .fq.gz .fastq .fq; do n=${n%"$ext"}; done
    [ -n "$STRIP" ] && n=${n%$STRIP}
  fi
  names+=("$n")
done
dups=$(printf '%s\n' "${names[@]}" | sort | uniq -d)
[ -z "$dups" ] || { echo "derep_ab.sh: duplicate sample names: $dups" >&2; exit 1; }
sample_csv=$(IFS=,; echo "${names[*]}")

{
  echo "binary: $BIN ($("$BIN" --version 2>&1 | head -1))"
  echo "host: $(hostname)  threads: $THREADS  derep_jobs: $DEREP_JOBS  samples: $#"
  echo "dada_args: ${DADA_ARGS:-}  prefix: ${PREFIX:-}"
} | tee "$OUT/run_info.txt"

# --- Arm A: FASTQ in -----------------------------------------------------------
if [ -z "${SKIP_A:-}" ]; then
  echo "==> arm A: dada-pooled on FASTQ"
  t0=$(now)
  ${prefix[@]+"${prefix[@]}"} "$BIN" dada-pooled --error-model "$ERR" \
    --output-dir "$OUT/a_fastq" --threads "$THREADS" --sample-names "$sample_csv" \
    --metrics-json "$OUT/a_fastq.metrics.json" --verbose \
    ${dada_extra[@]+"${dada_extra[@]}"} "$@" 2> "$OUT/logs/a_fastq.log"
  echo "$(now) $t0" | awk '{printf "%.3f\n", $1-$2}' > "$OUT/a_fastq.wall"
fi

# --- Derep step: one process per sample ---------------------------------------
echo "==> derep: $# samples, $DEREP_JOBS at a time"
: > "$OUT/derep/times.tsv"
derep_one() {
  local f="$1" n="$2" s e
  s=$(now)
  ${prefix[@]+"${prefix[@]}"} "$BIN" derep --sample-name "$n" \
    -o "$OUT/derep/$n.json.gz" "$f" 2> "$OUT/logs/derep_$n.log"
  e=$(now)
  printf '%s\t%s\n' "$n" "$(awk -v a="$s" -v b="$e" 'BEGIN{printf "%.3f", b-a}')" \
    >> "$OUT/derep/times.tsv"
}
t0=$(now)
i=0
for f in "$@"; do
  derep_one "$f" "${names[$i]}" &
  i=$((i + 1))
  # Throttle to DEREP_JOBS; `wait -n` needs bash 4.3, so fall back to batches.
  if [ $((i % DEREP_JOBS)) -eq 0 ]; then wait; fi
done
wait
echo "$(now) $t0" | awk '{printf "%.3f\n", $1-$2}' > "$OUT/derep.wall"

jsons=()
for n in "${names[@]}"; do jsons+=("$OUT/derep/$n.json.gz"); done

# --- Arm B: derep JSON in ------------------------------------------------------
echo "==> arm B: dada-pooled on derep JSON"
t0=$(now)
${prefix[@]+"${prefix[@]}"} "$BIN" dada-pooled --error-model "$ERR" \
  --output-dir "$OUT/b_json" --threads "$THREADS" --sample-names "$sample_csv" \
  --metrics-json "$OUT/b_json.metrics.json" --verbose \
  ${dada_extra[@]+"${dada_extra[@]}"} "${jsons[@]}" 2> "$OUT/logs/b_json.log"
echo "$(now) $t0" | awk '{printf "%.3f\n", $1-$2}' > "$OUT/b_json.wall"

# --- Compare and summarise -----------------------------------------------------
python3 -I - "$OUT" <<'PY'
import json, sys
from pathlib import Path

out = Path(sys.argv[1])

def load(p):
    return json.loads(p.read_text())

def wall(name):
    return float((out / f"{name}.wall").read_text())

# Outputs: identical except the arm-specific input name and the version string.
IGNORE = {"input_file", "dada2_rs_version"}
a_files = sorted(p.name for p in (out / "a_fastq").glob("*.json"))
b_files = sorted(p.name for p in (out / "b_json").glob("*.json"))
diffs = []
if a_files != b_files:
    diffs.append(f"file sets differ: A-only {set(a_files)-set(b_files)}, "
                 f"B-only {set(b_files)-set(a_files)}")
for name in sorted(set(a_files) & set(b_files)):
    a = {k: v for k, v in load(out / "a_fastq" / name).items() if k not in IGNORE}
    b = {k: v for k, v in load(out / "b_json" / name).items() if k not in IGNORE}
    if a != b:
        keys = sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))
        diffs.append(f"{name}: fields differ: {keys}")

ma, mb = load(out / "a_fastq.metrics.json"), load(out / "b_json.metrics.json")
pa, pb = ma["pipeline"], mb["pipeline"]
da, db = pa["derep_detail"], pb["derep_detail"]

times = [float(l.split("\t")[1]) for l in (out / "derep/times.tsv").read_text().splitlines()]
json_bytes = sum(p.stat().st_size for p in (out / "derep").glob("*.json.gz"))

def peak(p):
    r = p.get("peak_rss_mb") or {}
    return max(r.values()) if r else float("nan")

rows = [
    ("metric", "A_fastq", "B_json"),
    ("wall_seconds (process)", ma["wall_seconds"], mb["wall_seconds"]),
    ("wall_seconds (outside)", wall("a_fastq"), wall("b_json")),
    ("load_front (pipeline.derep)", pa["derep"], pb["derep"]),
    ("  read (I/O + gunzip)", "", db.get("read")),
    ("  parse (serde_json)", "", db.get("parse")),
    ("  build", da["build"], db["build"]),
    ("  per_sample median", da["per_sample"]["median"], db["per_sample"]["median"]),
    ("  per_sample max", da["per_sample"]["max"], db["per_sample"]["max"]),
    ("  input MB (FASTQ on disk / JSON raw)", da["bytes"] / 2**20, db["bytes"] / 2**20),
    ("dada", pa["dada"], pb["dada"]),
    ("load_front share of wall", pa["derep"] / ma["wall_seconds"],
     pb["derep"] / mb["wall_seconds"]),
    ("peak_rss_mb", peak(pa), peak(pb)),
    ("derep step: sum of per-sample", "", sum(times)),
    ("derep step: wall", "", wall("derep")),
    ("derep step: max sample", "", max(times)),
    ("derep JSON.gz MB on disk", "", json_bytes / 2**20),
]

def fmt(v):
    if v == "" or v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4f}" if abs(v) < 1 else f"{v:.2f}"
    return str(v)

lines = ["\t".join(fmt(c) for c in r) for r in rows]
lines.append(f"outputs\t{len(a_files)} samples\t"
             + ("IDENTICAL" if not diffs else f"{len(diffs)} DIFFER"))
text = "\n".join(lines) + "\n"
(out / "summary.tsv").write_text(text)
print(text, end="")
for d in diffs[:20]:
    print("  " + d, file=sys.stderr)
sys.exit(1 if diffs else 0)
PY
