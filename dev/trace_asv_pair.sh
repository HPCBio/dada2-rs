#!/usr/bin/env bash
# trace_asv_pair.sh — for a pair of ASVs that swap names between two
# member-order arms (#157): are both present in the sample, do they share one
# cluster, and could the non-centre sibling ever bud?
#
# Usage:
#   trace_pair_157.sh <binary> <err.json> <sample.fastq.gz> <armA> <armB> <seq1> <seq2> [out-dir]
# e.g. armA=insertion armB=shuffle:2
#
# Prints (1) each sequence's derep count and rank, (2) per arm, the cluster(s)
# holding either sequence: centre, cluster reads, and for each of the pair its
# reads, lambda, e_reads, and the abundance p-value it would need to bud, as
# P(X >= reads | e_reads) x nraw against OMEGA_A.
set -euo pipefail
BIN=${1:?binary}; ERR=${2:?err.json}; FQ=${3:?fastq}; A=${4:?armA}; B=${5:?armB}
S1=${6:?seq1}; S2=${7:?seq2}; OUT=${8:-trace_pair_out}
mkdir -p "$OUT"
[ -f "$FQ" ] || { echo "ERROR: no such FASTQ: $FQ" >&2; exit 1; }

# Each step's stderr goes to a log, shown if the step fails: silencing it once
# turned a missing input into a run that printed one header and stopped.
step() {
  local log=$1; shift
  if ! "$@" 2> "$log"; then
    echo "ERROR: step failed: $*" >&2; tail -n 20 "$log" >&2; exit 1
  fi
}
step "$OUT/derep.log" env -u DADA2RS_MEMBER_ORDER "$BIN" derep "$FQ" -o "$OUT/derep.json"
for arm in "$A" "$B"; do
  tag=${arm/:/_}
  step "$OUT/$tag.log" env DADA2RS_MEMBER_ORDER="$arm" "$BIN" dada "$FQ" --error-model "$ERR" \
    --threads 1 --cluster-trace "$OUT/$tag.trace.json" -o "$OUT/$tag.json"
done

python3 - "$OUT" "$A" "$B" "$S1" "$S2" <<'EOF'
import json, math, sys
out, A, B, s1, s2 = sys.argv[1:6]
pair = {s1.upper(): "seq1", s2.upper(): "seq2"}

def upper_tail(n, lam):
    """P(X >= n) for Poisson(lam), summed directly (stable for tiny lam)."""
    if lam <= 0:
        return 0.0 if n > 0 else 1.0
    return sum(math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1)) for k in range(n, n + 200))

u = json.load(open(f"{out}/derep.json"))["uniques"]
print(f"derep: {len(u)} uniques, {sum(x['count'] for x in u)} reads")
for rank, x in enumerate(u):
    if x["sequence"].upper() in pair:
        print(f"  {pair[x['sequence'].upper()]}: {x['count']} reads (derep rank {rank})")
missing = set(pair.values()) - {pair[x["sequence"].upper()] for x in u if x["sequence"].upper() in pair}
for m in sorted(missing):
    print(f"  {m}: NOT PRESENT in this sample's derep")

for arm in (A, B):
    t = json.load(open(f"{out}/{arm.replace(':', '_')}.trace.json"))
    seqs = [s.upper() for s in t["sequences"]]
    print(f"\n== {arm}: {t['nclust']} clusters, nraw {t['nraw']}, OMEGA_A {t['omega_a']:.0e}")
    for c in t["clusters"]:
        members = c.get("members") or []
        hit = [m for m in members if seqs[m["raw_seq_id"]] in pair]
        centre = seqs[c["center_seq_id"]]
        if not hit and centre not in pair:
            continue
        print(f"  cluster {c['id']} ({c['birth_type']}): {c['abundance']} reads, "
              f"centre = {pair.get(centre, 'other')}")
        for m in hit:
            name = pair[seqs[m["raw_seq_id"]]]
            if seqs[m["raw_seq_id"]] == centre:
                print(f"    {name}: {m['abundance']} reads  [centre]")
                continue
            p = upper_tail(m["abundance"], m["e_reads"])
            print(f"    {name}: {m['abundance']} reads, lambda {m['lambda']:.2e}, e_reads {m['e_reads']:.2e}, "
                  f"pA x nraw {p * t['nraw']:.1e} -> "
                  f"{'could bud' if p * t['nraw'] < t['omega_a'] else 'cannot bud (not significant)'}")
EOF
