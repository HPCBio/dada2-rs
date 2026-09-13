#!/usr/bin/env python3
"""compare_counts_correlation.py — is the count matrix distorted, and how?

`compare_seqtab_matrix.py` answers "do the counts differ, and by how much in
total". That total is a **global sum**, and three very different situations
produce the same number:

  * a small distortion spread evenly over every cell,
  * a few samples badly mangled while the rest are exact,
  * a systematic bias in one direction, which shifts every downstream
    composition even when the magnitude looks small.

This reports the shape of the disagreement rather than its size. On ASVs present
in BOTH tables it gives:

  * **log10 Pearson / Spearman** -- raw-count Pearson is useless here (abundant
    ASVs pin it near 1.0000 regardless), so counts are compared on a log scale
    and by rank, where a rare ASV counts as much as a dominant one.
  * **Bray-Curtis per sample** -- the dissimilarity ecologists actually act on.
    Reported as a distribution, so one wrecked sample cannot hide behind 361
    good ones.
  * **signed bias** -- whether the comparison arm systematically gains or loses
    reads. Random reassignment sums to ~0; a real bias does not.
  * **where the error lives** -- deviation split by ASV abundance decile, which
    separates "rare ASVs jitter" from "abundant ASVs moved".

Usage:
  compare_counts_correlation.py <a.json> <b.json> [--label-a A] [--label-b B]
"""

import argparse
import json
import math
import statistics
import sys


def load(path):
    j = json.load(open(path))
    seqs, samples, counts = j["sequences"], j["samples"], j["counts"]
    cell, per_asv = {}, {}
    for smp, row in zip(samples, counts):
        for q, c in zip(seqs, row):
            if c:
                cell[(smp, q)] = c
                per_asv[q] = per_asv.get(q, 0) + c
    return set(seqs), cell, per_asv, samples


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    args = ap.parse_args()

    sa, ca, aa, samples_a = load(args.a)
    sb, cb, ab, _ = load(args.b)
    la, lb = args.label_a, args.label_b
    shared = sa & sb

    print(f"{la}: {len(sa)} ASVs, {sum(ca.values()):,} reads")
    print(f"{lb}: {len(sb)} ASVs, {sum(cb.values()):,} reads")
    print(f"shared ASVs: {len(shared)}  (only_{la}={len(sa-sb)}, only_{lb}={len(sb-sa)})")

    # --- cell-level, restricted to shared ASVs ---
    xs, ys = [], []
    for (smp, q), v in ca.items():
        if q in shared:
            xs.append(v)
            ys.append(cb.get((smp, q), 0))
    for (smp, q), v in cb.items():
        if q in shared and (smp, q) not in ca:
            xs.append(0)
            ys.append(v)

    lx = [math.log10(v + 1) for v in xs]
    ly = [math.log10(v + 1) for v in ys]
    print(f"\ncells compared (shared ASVs): {len(xs):,}")
    print(f"  Pearson  r  on RAW counts : {pearson(xs, ys):.6f}   <- near 1 regardless; not diagnostic")
    print(f"  Pearson  r  on log10 counts: {pearson(lx, ly):.6f}")
    print(f"  Spearman rho (rank)        : {spearman(xs, ys):.6f}")

    exact = sum(1 for x, y in zip(xs, ys) if x == y)
    within10 = sum(1 for x, y in zip(xs, ys) if max(x, y) and abs(x - y) <= 0.10 * max(x, y))
    within2x = sum(1 for x, y in zip(xs, ys) if max(x, y) and max(x, y) <= 2 * max(min(x, y), 1))
    print(f"  cells exactly equal        : {100*exact/len(xs):.3f}%")
    print(f"  cells within 10%           : {100*within10/len(xs):.3f}%")
    print(f"  cells within 2x            : {100*within2x/len(xs):.3f}%")

    # --- signed bias: random reassignment cancels, a real bias does not ---
    net = sum(y - x for x, y in zip(xs, ys))
    gross = sum(abs(y - x) for x, y in zip(xs, ys))
    tot = sum(xs)
    print(f"\nsigned bias on shared ASVs:")
    print(f"  net   change: {net:+,} reads ({100*net/max(tot,1):+.4f}% of {la})")
    print(f"  gross change: {gross:,} reads ({100*gross/max(tot,1):.4f}% of {la})")
    if gross:
        print(f"  net/gross   : {abs(net)/gross:.3f}   "
              f"({'SYSTEMATIC -- one direction' if abs(net)/gross > 0.3 else 'cancels out; reassignment, not bias'})")

    # --- per-sample Bray-Curtis: one wrecked sample cannot hide in an average ---
    bc = []
    for smp in samples_a:
        num = den = 0
        for q in shared:
            x, y = ca.get((smp, q), 0), cb.get((smp, q), 0)
            num += abs(x - y)
            den += x + y
        if den:
            bc.append((num / den, smp))
    if bc:
        vals = sorted(v for v, _ in bc)
        bc.sort(reverse=True)
        print(f"\nper-sample Bray-Curtis dissimilarity ({len(bc)} samples):")
        print(f"  median {statistics.median(vals):.6f}   p90 {vals[int(0.9*len(vals))]:.6f}   max {vals[-1]:.6f}")
        print(f"  samples with BC > 0.01 : {sum(1 for v in vals if v > 0.01)}")
        print(f"  worst 3: " + ", ".join(f"{s} ({v:.4f})" for v, s in bc[:3]))

    # --- where the error lives, by ASV abundance ---
    ranked = sorted(shared, key=lambda q: -aa.get(q, 0))
    dec = max(1, len(ranked) // 10)
    print(f"\ndeviation by ASV abundance decile (1 = most abundant):")
    print(f"  {'decile':>7s} {'ASVs':>5s} {'reads in ' + la:>14s} {'gross change':>13s} {'% of decile':>12s}")
    for d in range(10):
        grp = ranked[d*dec:(d+1)*dec] if d < 9 else ranked[9*dec:]
        if not grp:
            continue
        g = set(grp)
        tot_d = sum(v for (s, q), v in ca.items() if q in g)
        ch = 0
        for (s, q), v in ca.items():
            if q in g:
                ch += abs(cb.get((s, q), 0) - v)
        for (s, q), v in cb.items():
            if q in g and (s, q) not in ca:
                ch += v
        print(f"  {d+1:7d} {len(grp):5d} {tot_d:14,d} {ch:13,d} {100*ch/max(tot_d,1):11.4f}%")


if __name__ == "__main__":
    main()
