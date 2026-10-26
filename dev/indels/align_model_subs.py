#!/usr/bin/env python3
"""align_model_subs.py — where does a gap-opening cost move substitutions? (#309)

Under affine gaps a gap-opening cost changes which alignment wins, not just
its score: a nearby insertion/deletion pair can become substitutions, and an
end gap can become terminal mismatches. This realigns a fixed partition (each
unique against its final centre) under two models and classifies every
substitution one model calls and the other does not, so the shift can be
mapped along the amplicon and judged against the substitutions both models
agree on.

  A  ends-free, linear gaps (dada2-rs NW default): end gaps free
  B  global, affine gaps: a gap of length L costs gap*L + open

Scoring as dada2-rs: match 5, mismatch -4, gap -8. Fixing the partition
separates the alignment effect from re-partitioning across iterations, which
the learn-errors totals mix together.

Inputs: a final-iteration cluster trace (learn-errors --cluster-trace-dir)
and the derep JSON for the same reads (per-position quality sums).

Check first (--check): the read-weighted substitution total under model A
should match the learned model's transition counts (err.json `trans`);
alignment ties may move a few.

Classes for a substitution only one model calls:
  pair     the other model has an insertion and a deletion within --window
           bases of it (an indel pair re-described as substitutions)
  end      within --window bases of an end where the other model has an end gap
  single   the other model has one indel within --window (a gap moved)
  other    none of the above

Requires Biopython (PairwiseAligner).

Usage:
  align_model_subs.py TRACE.json DEREP.json --open 8 [--check err.json] \\
      [--out events.tsv] [--threads 8]
"""

import argparse
import collections
import csv
import json
import multiprocessing as mp
import sys

from Bio.Align import PairwiseAligner

MATCH, MISMATCH, GAP = 5, -4, -8
_aligners = {}


def aligner(model, gap_open):
    key = (model, gap_open)
    if key not in _aligners:
        a = PairwiseAligner(mode="global", match_score=MATCH, mismatch_score=MISMATCH)
        if model == "A":
            a.open_gap_score = a.extend_gap_score = GAP
            a.end_gap_score = 0
        else:
            a.open_gap_score = GAP - gap_open
            a.extend_gap_score = GAP
        _aligners[key] = a
    return _aligners[key]


def events(center, read, model, gap_open):
    """Substitutions {centre pos: read pos}, gap positions (centre coords), end gaps."""
    aln = aligner(model, gap_open).align(center, read)[0]
    c_al, r_al = aln[0], aln[1]
    subs, ins, dels, gaps = {}, [], [], []
    ci = ri = 0
    cols = []  # alignment columns holding a base in both sequences
    for k, (x, y) in enumerate(zip(c_al, r_al)):
        if x != "-" and y != "-":
            cols.append(k)
            if x != y:
                subs[ci] = ri
            ci += 1
            ri += 1
        else:
            gaps.append(k)
            if x == "-":
                ins.append((k, ci))  # inserted before centre position ci
                ri += 1
            else:
                dels.append((k, ci))
                ci += 1
    first, last = (cols[0], cols[-1]) if cols else (len(c_al), -1)
    end5 = any(k < first for k in gaps)
    end3 = any(k > last for k in gaps)
    internal_ins = [p for k, p in ins if first < k < last]
    internal_del = [p for k, p in dels if first < k < last]
    return subs, internal_ins, internal_del, end5, end3


def classify(pos, other, n, window):
    _, oins, odel, end5, end3 = other
    near_i = any(abs(p - pos) <= window for p in oins)
    near_d = any(abs(p - pos) <= window for p in odel)
    if near_i and near_d:
        return "pair"
    if (end5 and pos < window) or (end3 and pos >= n - window):
        return "end"
    if near_i or near_d:
        return "single"
    return "other"


def run_length(seq, i):
    if i < 0 or i >= len(seq):
        return 0
    b, lo, hi = seq[i], i, i
    while lo > 0 and seq[lo - 1] == b:
        lo -= 1
    while hi + 1 < len(seq) and seq[hi + 1] == b:
        hi += 1
    return hi - lo + 1


def work(job):
    center, read, reads, qual, gap_open, window = job
    if center == read:
        return reads, 0, 0, []
    a = events(center, read, "A", gap_open)
    b = events(center, read, "B", gap_open)
    n = len(center)
    rows = []
    for tag, mine, other in (("B_only", b, a), ("A_only", a, b), ("both", b, None)):
        for pos, rpos in mine[0].items():
            if tag == "both":
                if pos not in a[0]:
                    continue
                cls = "both"
            elif pos in other[0]:
                continue
            else:
                cls = classify(pos, other, n, window)
            rows.append(dict(set=tag, cls=cls, pos=pos, to_end=min(pos, n - 1 - pos), len=n,
                             cbase=center[pos], rbase=read[rpos], qual=qual[rpos] if qual else "",
                             hp_left=run_length(center, pos - 1) if pos > 0 and center[pos - 1] == read[rpos] else 0,
                             hp_right=run_length(center, pos + 1) if pos + 1 < n and center[pos + 1] == read[rpos] else 0,
                             reads=reads))
    return reads, len(a[0]), len(b[0]), rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("trace", help="final-iteration cluster trace JSON")
    ap.add_argument("derep", help="derep JSON for the same reads (quality sums)")
    ap.add_argument("--open", type=float, required=True, help="model B gap-opening cost")
    ap.add_argument("--window", type=int, default=10, help="bases for pair/end/single classes")
    ap.add_argument("--check", help="err.json whose trans substitution total model A should match")
    ap.add_argument("--out", help="TSV of every differing (and shared) substitution")
    ap.add_argument("--bin", type=int, default=100, help="position bin for the map")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    t = json.load(open(a.trace))
    d = json.load(open(a.derep))
    qual = {}
    for u in d["uniques"]:
        qs = json.loads(u["qual_sum"]) if isinstance(u["qual_sum"], str) else u["qual_sum"]
        qual[u["sequence"]] = [q / int(u["count"]) for q in qs]
    seqs = t["sequences"]
    jobs = []
    for c in t["clusters"]:
        center = seqs[c["center_seq_id"]]
        for m in c["members"]:
            r = seqs[m["raw_seq_id"]]
            jobs.append((center, r, m["abundance"], qual.get(r), a.open, a.window))

    tot_a = tot_b = 0
    allrows = []
    with mp.Pool(a.threads) as pool:
        for reads, na, nb, rows in pool.imap_unordered(work, jobs, chunksize=16):
            tot_a += reads * na
            tot_b += reads * nb
            allrows.extend(rows)

    print(f"{len(jobs)} uniques; read-weighted substitutions: A (ends-free linear) {tot_a}, "
          f"B (global, open {a.open:g}) {tot_b} ({(tot_b - tot_a) / max(tot_a, 1) * 100:+.1f}%)")
    if a.check:
        tr = json.load(open(a.check))["trans"]
        ref = sum(sum(tr[r]) for r in range(16) if r // 4 != r % 4)
        print(f"check: err.json trans substitutions {ref}; model A {tot_a} ({(tot_a - ref) / max(ref, 1) * 100:+.2f}%)")

    w = collections.Counter()
    for r in allrows:
        w[(r["set"], r["cls"])] += r["reads"]
    print("\nread-weighted substitutions by set and class (class = what the other model did nearby):")
    for s in ("B_only", "A_only", "both"):
        tot = sum(v for (ss, _), v in w.items() if ss == s)
        parts = ", ".join(f"{c} {w[(s, c)]}" for c in ("pair", "single", "end", "other", "both") if w[(s, c)])
        print(f"  {s:7s} {tot:7d}  ({parts})")

    def qstats(rows):
        qs = sorted(r["qual"] for r in rows for _ in range(1) if r["qual"] != "")
        if not qs:
            return "-"
        wq = [(r["qual"], r["reads"]) for r in rows if r["qual"] != ""]
        tot = sum(n for _, n in wq)
        mean = sum(q * n for q, n in wq) / tot
        return f"mean Q {mean:5.1f}"

    def shift_frac(rows):
        tot = sum(r["reads"] for r in rows)
        sh = sum(r["reads"] for r in rows if r["hp_left"] or r["hp_right"])
        return f"read base extends a neighbouring run {sh / max(tot, 1) * 100:5.1f}%"

    print("\nfeatures (read-weighted):")
    groups = [("both", lambda r: r["set"] == "both")]
    for c in ("pair", "single", "end", "other"):
        groups.append((f"B_only/{c}", lambda r, c=c: r["set"] == "B_only" and r["cls"] == c))
    groups.append(("A_only", lambda r: r["set"] == "A_only"))
    for name, f in groups:
        rows = [r for r in allrows if f(r)]
        if rows:
            print(f"  {name:14s} {sum(r['reads'] for r in rows):7d}  {qstats(rows)}  {shift_frac(rows)}")

    print(f"\nmap: read-weighted B_only substitutions per {a.bin}-nt bin of centre position (pair/single/end/other)")
    bins = collections.defaultdict(collections.Counter)
    shared = collections.Counter()
    for r in allrows:
        if r["set"] == "B_only":
            bins[r["pos"] // a.bin][r["cls"]] += r["reads"]
        elif r["set"] == "both":
            shared[r["pos"] // a.bin] += r["reads"]
    for b in sorted(set(bins) | set(shared)):
        c = bins[b]
        print(f"  {b * a.bin:5d}-{(b + 1) * a.bin - 1:<5d} B_only {sum(c.values()):6d} "
              f"({c['pair']}/{c['single']}/{c['end']}/{c['other']})  shared {shared[b]:6d}")

    if a.out:
        with open(a.out, "w", newline="") as fh:
            cols = ["set", "cls", "pos", "to_end", "len", "cbase", "rbase", "qual", "hp_left", "hp_right", "reads"]
            wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            wr.writeheader()
            for r in allrows:
                r = dict(r)
                if r["qual"] != "":
                    r["qual"] = f"{r['qual']:.1f}"
                wr.writerow(r)


if __name__ == "__main__":
    sys.exit(main())
