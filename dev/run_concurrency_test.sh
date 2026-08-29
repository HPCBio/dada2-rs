#!/usr/bin/env bash
# run_concurrency_test.sh — is it better to run one wide job or two narrow ones?
# ---------------------------------------------------------------------------
# WHY: the #152 thread sweep found that every arm measured buys wall time at
# roughly 4x the core-seconds — 16S R1 goes 1072.5s at 24 threads to 594.3s at
# 96, i.e. 45% parallel efficiency at the top. That says "pack jobs, don't scale
# threads". But nobody has measured two jobs sharing a node, and they would
# contend for exactly the memory bandwidth the k-mer screen is already limited
# by. The packing recommendation is currently an inference, not a result.
#
# THE NUMA QUESTION, which is the interesting half: dev/numa_pin.sh documents
# that binding to a single domain is the WRONG choice for one job — it forces
# every thread through one node's memory controllers and costs the parallel map
# 21%. With two concurrent jobs on a 2-domain node that reasoning may invert:
# each job gets a private set of controllers and neither disturbs the other.
# Interleaving both jobs, by contrast, has them share every controller. Which
# wins is not predictable from the single-job result, so both are measured.
#
# WHY THE SAME READ TWICE: pairing R1 with R2 (713.8s vs 767.8s at 48 threads
# on 16S) lets the shorter job finish first, leaving the longer one's tail
# running uncontended — which flatters the pair and blurs the contention it is
# meant to measure. Two copies of one read contend end to end.
#
# ARMS
#   solo   one job, N threads, --interleave=all          (the #152 baseline)
#   sbind  one job, N threads, bound to ONE NUMA domain  (locality, no sharing)
#   both   two jobs, N threads each, both interleaved    (naive packing)
#   split  two jobs, N threads each, each bound to its own NUMA domain
#
# Contention penalty = (arm wall) / (solo wall) - 1, per job. Ideal is 0%:
# two jobs finishing in the time one takes means the node was not the limit.
#
# Usage:
#   dev/run_concurrency_test.sh \
#       --bin ./target/release-native/dada2-rs \
#       --error-model err.json \
#       --out-root tmp/issue-152/concurrency \
#       --threads 48 \
#       derep/*.json
set -euo pipefail

BIN=./target/release-native/dada2-rs
ERRMODEL=
OUT_ROOT=concurrency-test
THREADS=48
ARMS="solo sbind both split"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)          BIN=$2; shift 2 ;;
    --error-model)  ERRMODEL=$2; shift 2 ;;
    --out-root)     OUT_ROOT=$2; shift 2 ;;
    --threads)      THREADS=$2; shift 2 ;;
    --arms)         ARMS=$2; shift 2 ;;
    --) shift; break ;;
    -*) echo "unknown option: $1" >&2; exit 2 ;;
    *) break ;;
  esac
done
INPUTS=("$@")

[[ -n "$ERRMODEL" ]] || { echo "--error-model is required" >&2; exit 2; }
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "no input derep files given" >&2; exit 2; }
command -v numactl >/dev/null || { echo "numactl not found" >&2; exit 2; }

NODES=$(numactl --hardware | awk '/^available:/ {print $2}')
echo "== node has ${NODES} NUMA domain(s); threads/job = ${THREADS}"
if [[ "$NODES" -lt 2 ]]; then
  echo "== only one NUMA domain: the 'split' arm is meaningless here, dropping it"
  ARMS=${ARMS//split/}
fi
mkdir -p "$OUT_ROOT"

# One job. $1 label, $2 numactl args. Writes <label>.log and <label>.wall.
run_job() {
  local label=$1 numa=$2 start end
  start=$(date +%s.%N)
  # shellcheck disable=SC2086
  DADA2RS_PROGRESS_SECS=15 numactl $numa "$BIN" dada-pooled \
      --verbose --threads "$THREADS" \
      --error-model "$ERRMODEL" \
      --output-dir "$OUT_ROOT/$label" \
      "${INPUTS[@]}" > "$OUT_ROOT/$label.log" 2>&1
  end=$(date +%s.%N)
  echo "$end - $start" | bc > "$OUT_ROOT/$label.wall"
}

pair_wall() {  # epoch seconds spanning both jobs, i.e. what the node was busy
  local a=$1 b=$2
  python3 - "$OUT_ROOT" "$a" "$b" <<'PY'
import sys, pathlib
root, a, b = sys.argv[1], sys.argv[2], sys.argv[3]
w = [float(pathlib.Path(f"{root}/{x}.wall").read_text()) for x in (a, b)]
print(f"{max(w):.1f}")
PY
}

for arm in $ARMS; do
  case "$arm" in
    solo)
      echo "== solo: one job, ${THREADS} threads, interleaved"
      run_job solo "--interleave=all"
      ;;
    sbind)
      # Separates LOCALITY from CONCURRENCY. On ITS2 R1 the split arm ran each
      # job 28% faster than solo -- but those jobs were still sharing the node,
      # so locality and packing are confounded there. If sbind also lands near
      # split, the whole gain is locality and packing is incidental; if it lands
      # near solo, the gain needs the second job.
      echo "== sbind: one job, ${THREADS} threads, bound to NUMA domain 0"
      run_job sbind "--cpunodebind=0 --membind=0"
      ;;
    both)
      echo "== both: two jobs, ${THREADS} threads each, both interleaved"
      run_job both_a "--interleave=all" &
      run_job both_b "--interleave=all" &
      wait
      ;;
    split)
      echo "== split: two jobs, ${THREADS} threads each, one NUMA domain apiece"
      # cpunodebind pins CPUs as well as memory, so THREADS must fit one domain.
      run_job split_a "--cpunodebind=0 --membind=0" &
      run_job split_b "--cpunodebind=1 --membind=1" &
      wait
      ;;
  esac
done

echo
echo "== results (wall seconds per job)"
solo_w=$(cat "$OUT_ROOT/solo.wall" 2>/dev/null || echo 0)
printf "%-10s %10s %10s %14s\n" arm job_a job_b "vs solo"
printf "%-10s %10.1f %10s %14s\n" solo "$solo_w" "-" "-"
for arm in sbind; do
  [[ -f "$OUT_ROOT/${arm}.wall" ]] || continue
  a=$(cat "$OUT_ROOT/${arm}.wall")
  pen=$(python3 -c "print(f'{($a/$solo_w-1)*100:+.1f}%')" 2>/dev/null || echo "n/a")
  printf "%-10s %10.1f %10s %14s\n" "$arm" "$a" "-" "$pen"
done
for arm in both split; do
  [[ -f "$OUT_ROOT/${arm}_a.wall" ]] || continue
  a=$(cat "$OUT_ROOT/${arm}_a.wall"); b=$(cat "$OUT_ROOT/${arm}_b.wall")
  pw=$(pair_wall "${arm}_a" "${arm}_b")
  pen=$(python3 -c "print(f'{(max($a,$b)/$solo_w-1)*100:+.1f}%')" 2>/dev/null || echo "n/a")
  printf "%-10s %10.1f %10.1f %14s   (pair wall %ss)\n" "$arm" "$a" "$b" "$pen" "$pw"
done
echo
echo "Throughput comparison: two jobs in 'pair wall' seconds against one job in"
echo "'solo' seconds. Packing wins if pair wall < 2x solo — and by how much is"
echo "the number the docs need. A penalty near 0% means the node was not the"
echo "limit and threads were being wasted on a single job."
echo
echo "NOTE: outputs across arms must be byte-identical. Concurrency and NUMA"
echo "policy change timing only; if they change results, something is wrong and"
echo "no timing here is worth reading."
