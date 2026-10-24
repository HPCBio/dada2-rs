#!/usr/bin/env python3
"""hp_variant_counts.py — homopolymer slippage rate per run, per parent ASV (issue #291).

For each deep parent ASV, counts the reads that match it exactly and the reads
that match it exactly except for one base added to or removed from one
homopolymer run. The ratio is that run's slippage rate, measured on reads
otherwise free of errors, so parents compare even when their depths differ.
Variant reads are counted whether or not dada split them: in a default run the
variants are absorbed and invisible in the table.

Use it to tell a sequence-context error (similar rates wherever the context
recurs) from something template-specific (high rates in some taxa, near zero
in others with the same context), e.g. the PacBio 95 hotspots: +G at
`GCAGC|GGGGG` (site A, V1) and `GGC|GGGG|TAAC` (site B, V2).

Reads must be prepped as for dada (primers trimmed, oriented), so an
error-free read equals its ASV. Reads truncated to one length (Illumina
--trunc-len) are detected: an insertion read then lacks its last base and a
deletion read carries one extra, so variants are matched at the read length.

Inputs:
  --seqtab     chimera-free table (run_indel_ab.sh <out>/default/nochim.json)
  --manifest   sample<TAB>fastq (run_indel_ab.sh's manifest)
  --taxonomy   optional assign-taxonomy JSON for the table's sequences

Output (--out TSV), one row per parent x run x edit:
  parent, parent_reads (table), parent_exact (exact reads), family, genus,
  run_pos, base, run_len, edit (+1/-1), variant_reads, frac
  (variant / (variant + parent_exact)), samples (with the variant), site (A, B
  or -), left, right (10 nt of context), to_end (bases after the run). Under
  fixed-length reads a run near the read end is unreliable: an insertion there
  is the tail shifted by one, indistinguishable from 3'-end errors.

Usage:
  hp_variant_counts.py --seqtab out/default/nochim.json --manifest out/manifest.tsv \\
      [--taxonomy tax.json] [--min-reads 3000] [--min-run 3] [--threads 8] --out hp_runs.tsv
"""

import argparse
import collections
import csv
import gzip
import json
import re
import sys
from multiprocessing import Pool

SITE_A = re.compile(r"AAGTCGA.{0,10}?GCAGC(G+)[ACT]")
SITE_B = re.compile(r"GGC(G+)TAACGGCCCACC")


def runs(seq, min_run):
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= min_run:
            yield i, j - i + 1
        i = j + 1


def site_of(seq, pos):
    for name, rx, lo, hi in (("A", SITE_A, 0, 150), ("B", SITE_B, 150, 400)):
        m = rx.search(seq, lo, hi)
        if m and m.start(1) == pos:
            return name
    return "-"


TARGETS = frozenset()
DEL_PREFIXES = frozenset()


def init(targets, del_prefixes):
    global TARGETS, DEL_PREFIXES
    TARGETS, DEL_PREFIXES = targets, del_prefixes


def count_file(path):
    hits = collections.Counter()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for k, line in enumerate(fh):
            if k % 4 == 1:
                s = line.rstrip("\n")
                if s in TARGETS:
                    hits[s] += 1
                elif DEL_PREFIXES and s[:-1] in DEL_PREFIXES:
                    hits[s[:-1]] += 1
    return hits


def fixed_length(parents, fastq, n=2000):
    """The common length when every parent and the first reads share one."""
    lens = {len(p) for p in parents}
    if len(lens) != 1:
        return None
    opener = gzip.open if fastq.endswith(".gz") else open
    with opener(fastq, "rt") as fh:
        seen = {len(line.rstrip("\n")) for k, line in enumerate(fh) if k % 4 == 1 and k < 4 * n}
    return lens.pop() if seen == lens else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seqtab", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--taxonomy")
    ap.add_argument("--min-reads", type=int, default=3000, help="parents with at least this many table reads")
    ap.add_argument("--min-run", type=int, default=3)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    t = json.load(open(a.seqtab))
    reads = {s: sum(row[i] for row in t["counts"]) for i, s in enumerate(t["sequences"])}
    parents = [s for s, n in sorted(reads.items(), key=lambda x: -x[1]) if n >= a.min_reads]
    tax = {}
    if a.taxonomy:
        tax = {x["sequence"]: x["taxonomy"] for x in json.load(open(a.taxonomy))["assignments"]}

    # Target string -> (parent index, run start, run length, edit). A string
    # reachable two ways (two parents, or two edits of one) is dropped: its
    # reads cannot be attributed.
    samples = [l.rstrip("\n").split("\t")[:2] for l in open(a.manifest) if l.strip()]
    fixed = fixed_length(parents, samples[0][1])
    print(f"read lengths: {'fixed at ' + str(fixed) if fixed else 'variable'}", file=sys.stderr)

    owner, clash, del_prefixes = {}, set(), set()
    def claim(s, key):
        if s in owner and owner[s] != key:
            clash.add(s)
        owner.setdefault(s, key)
    for pi, p in enumerate(parents):
        claim(p, (pi, None, None, 0))
        for pos, n in runs(p, a.min_run):
            ins, dele = p[:pos] + p[pos] + p[pos:], p[:pos] + p[pos + 1:]
            if fixed:
                # Truncated reads: the insertion loses its last base; the
                # deletion gains an unknown one, matched as a prefix.
                ins = ins[:fixed]
                del_prefixes.add(dele)
            claim(ins, (pi, pos, n, +1))
            claim(dele, (pi, pos, n, -1))
    for s in clash:
        del owner[s]
    print(f"{len(parents)} parents with >= {a.min_reads} reads; {len(owner)} target sequences "
          f"({len(clash)} ambiguous dropped)", file=sys.stderr)

    targets = frozenset(owner)
    prefixes = frozenset(d for d in del_prefixes if d in owner)
    total = collections.Counter()
    per_sample = collections.defaultdict(set)
    with Pool(a.threads, initializer=init, initargs=(targets, prefixes)) as pool:
        for (name, _), hits in zip(samples, pool.imap(count_file, [f for _, f in samples])):
            total.update(hits)
            for s in hits:
                per_sample[s].add(name)
            print(f"  {name}: {sum(hits.values())} target reads", file=sys.stderr)

    exact = {owner[s][0]: total[s] for s in owner if owner[s][1] is None}
    rows = []
    for s, (pi, pos, n, edit) in owner.items():
        if pos is None:
            continue
        p = parents[pi]
        v, e = total[s], exact.get(pi, 0)
        tx = tax.get(p, [""] * 6)
        rows.append(dict(parent=f"p{pi + 1}", parent_reads=reads[p], parent_exact=e,
                         family=tx[4] if len(tx) > 4 else "", genus=tx[5] if len(tx) > 5 else "",
                         run_pos=pos, base=p[pos], run_len=n, edit=f"{edit:+d}", variant_reads=v,
                         frac=f"{v / (v + e):.4g}" if v + e else "nan", samples=len(per_sample[s]),
                         site=site_of(p, pos), left=p[max(0, pos - 10):pos], right=p[pos + n:pos + n + 10],
                         to_end=len(p) - pos - n))
    rows.sort(key=lambda r: (int(r["parent"][1:]), r["run_pos"], r["edit"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print("\nsites A and B, +1 insertion, by parent:", file=sys.stderr)
    for r in rows:
        if r["site"] != "-" and r["edit"] == "+1":
            print(f"  {r['parent']:>4} {r['parent_reads']:>7} reads ({r['parent_exact']:>6} exact) "
                  f"{(r['family'] or '?')[:18]:18s} {(r['genus'] or '')[:16]:16s} site {r['site']} "
                  f"{r['base']}{r['run_len']}: {r['variant_reads']:>5} reads, frac {r['frac']}", file=sys.stderr)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
