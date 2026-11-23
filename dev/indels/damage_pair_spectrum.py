#!/usr/bin/env python3
"""damage_pair_spectrum.py — molecular vs sequencing substitution rates from read pairs (issue #291).

Usage:
  damage_pair_spectrum.py EVENTS.tsv [--overlap-start N] [--max-site-frac 0.01]

Reads the output of `damage_pair_events.py`. Within each pair's overlap (the
reference positions both reads cover), a substitution in both reads is
molecular: template damage or a PCR error, copied into both strands. One in a
single read is a sequencing error. Rates are events per reference base of the
substituted type in the overlap, x 1e-5.

--overlap-start clips the overlap's start, to keep out the low-quality start of
R2 or a primer-adjacent region. Sites (ref, position, new base) carried by more
than --max-site-frac of a reference's pairs are dropped as unlisted alleles or
position-specific errors, from both the events and that type's opportunities.
Class rates pool their types: events over bases of either source type.
"""
import argparse, collections

NT = "ACGT"
GROUP = {"C>T": "deamination", "G>A": "deamination", "G>T": "oxidation", "C>A": "oxidation",
         "A>G": "other transition", "T>C": "other transition"}
GROUPS = ("deamination", "oxidation", "other transition", "other transversion")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("events")
    ap.add_argument("--overlap-start", type=int, default=0)
    ap.add_argument("--max-site-frac", type=float, default=0.01)
    a = ap.parse_args()

    refs = {}
    for l in open(a.events + ".refs"):
        j, _, s = l.rstrip("\n").split("\t")
        refs[int(j)] = s
    # cumulative base counts per reference, for opportunities over any span
    cum = {}
    for j, s in refs.items():
        c = {b: [0] for b in NT}
        for x in s:
            for b in NT:
                c[b].append(c[b][-1] + (x == b))
        cum[j] = c

    span, npairs = {}, collections.Counter()
    with open(a.events + ".pairs") as fh:
        next(fh)
        for l in fh:
            pid, j, s1, e1, s2, e2 = map(int, l.split())
            lo, hi = max(s1, s2, a.overlap_start), min(e1, e2)
            if hi > lo:
                span[pid] = (j, lo, hi)
            npairs[j] += 1
    opp = collections.Counter()
    for j, lo, hi in span.values():
        for b in NT:
            opp[b] += cum[j][b][hi] - cum[j][b][lo]

    sides = collections.defaultdict(set)
    with open(a.events) as fh:
        next(fh)
        for l in fh:
            pid, j, pos, side, x, y, _ = l.rstrip("\n").split("\t")
            pid, j, pos = int(pid), int(j), int(pos)
            sp = span.get(pid)
            if sp and sp[1] <= pos < sp[2]:
                sides[(pid, j, pos, x, y)].add(side)
    carriers = collections.Counter((j, pos, y) for (_, j, pos, _, y) in sides)
    drop = {k for k, n in carriers.items() if n > a.max_site_frac * npairs[k[0]]}

    # per-type opportunities, less the pairs covering each dropped site
    optype = collections.Counter({f"{x}>{y}": opp[x] for x in NT for y in NT if x != y})
    for j, pos, y in drop:
        cover = sum(1 for jj, lo, hi in span.values() if jj == j and lo <= pos < hi)
        optype[f"{refs[j][pos]}>{y}"] -= cover

    mol, r1, r2 = collections.Counter(), collections.Counter(), collections.Counter()
    for (pid, j, pos, x, y), s in sides.items():
        if (j, pos, y) in drop:
            continue
        t = f"{x}>{y}"
        if len(s) == 2:
            mol[t] += 1
        elif "R1" in s:
            r1[t] += 1
        else:
            r2[t] += 1

    print(f"{len(span)} pairs with an overlap; {len(drop)} sites dropped above {a.max_site_frac:g} of pairs")
    print(f"{'type':18s} {'opp':>11s} {'molecular':>10s} {'R1 only':>9s} {'R2 only':>9s}  (x 1e-5)")
    rate = lambda c, ts: 1e5 * sum(c[t] for t in ts) / max(1, sum(optype[t] for t in ts))
    types = [f"{x}>{y}" for x in NT for y in NT if x != y]
    for t in types:
        print(f"{t:18s} {optype[t]:11d} {rate(mol, [t]):10.2f} {rate(r1, [t]):9.2f} {rate(r2, [t]):9.2f}")
    print()
    for g in GROUPS:
        ts = [t for t in types if GROUP.get(t, "other transversion") == g]
        print(f"{g:18s} {sum(optype[t] for t in ts):11d} {rate(mol, ts):10.2f} {rate(r1, ts):9.2f} {rate(r2, ts):9.2f}")

if __name__ == "__main__":
    main()
