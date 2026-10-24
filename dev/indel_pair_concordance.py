#!/usr/bin/env python3
"""indel_pair_concordance.py — molecular vs sequencing indels from read pairs (issue #291).

A sequencing error happens independently in the two reads of a pair; an error
in the molecule (PCR slippage, a template variant) appears in both. For each
pair whose forward and reverse reads denoise to the same reference sequence,
every indel inside the reads' overlap is classified as concordant (in both
reads) or one-read-only. Homopolymer events are keyed by the run they change,
so where the aligner places a gap inside the run, which can differ between
the strands, does not matter.

Inputs per direction: the `--indel-events` table, the dada JSON (cluster map
and centres) and the derep JSON with `--show-map` (read -> unique), for the
same reads in the same order. Pass a file for one sample, or a directory of
`<sample>.json` for several (the events table's `sample` column selects).
The reference is a FASTA of amplicons with primers removed: a mock's truth
set, or a run's own non-chimeric merged ASVs. A forward centre matches a
reference that starts with it, a reverse centre one that ends with its
reverse complement; reads truncated to --trunc-fwd / --trunc-rev overlap over
reference positions [len - trunc_rev, trunc_fwd).

Usage:
  indel_pair_concordance.py --truth truth.fa --trunc-fwd 240 --trunc-rev 200 \\
      --fwd-events indel_fwd.events.tsv --rev-events indel_rev.events.tsv \\
      --fwd-dada dada_fwd/ --rev-dada dada_rev/ --fwd-derep derepmap_fwd/ --rev-derep derepmap_rev/
"""

import argparse
import collections
import csv
import glob
import json
import os

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def by_sample(path):
    """{sample: json} from one JSON file or a directory of them."""
    paths = [path] if os.path.isfile(path) else sorted(glob.glob(os.path.join(path, "*.json")))
    out = {}
    for p in paths:
        d = json.load(open(p))
        out[d["sample"]] = d
    return out


def run_of(seq, i):
    b, s, e = seq[i], i, i
    while s > 0 and seq[s - 1] == b:
        s -= 1
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return s, e


def event_key(e, d, ref):
    """(type, start, end, bases, len, run length) in reference coordinates."""
    L, k, p = len(ref), int(e["len"]), int(e["center_pos"])
    bases = e["bases"] if d == "fwd" else rc(e["bases"])
    if e["type"] == "D":
        x = p if d == "fwd" else L - p - k
        if len(set(bases)) == 1:
            s, t = run_of(ref, x)
            return ("D", s, t, bases[0], k, t - s + 1)
        return ("D", x, x + k - 1, bases, k, 0)
    x = p if d == "fwd" else L - 2 - p          # inserted between x and x + 1
    if len(set(bases)) == 1:
        for y in (x, x + 1):
            if 0 <= y < L and ref[y] == bases[0]:
                s, t = run_of(ref, y)
                return ("I", s, t, bases[0], k, t - s + 1)
    return ("I", x, x, bases, k, 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--trunc-fwd", type=int, required=True)
    ap.add_argument("--trunc-rev", type=int, required=True)
    for d in ("fwd", "rev"):
        ap.add_argument(f"--{d}-events", required=True)
        ap.add_argument(f"--{d}-dada", required=True)
        ap.add_argument(f"--{d}-derep", required=True)
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()

    refs, name = {}, None
    for line in open(a.truth):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split(";")[0]
        elif name:
            refs[name] = refs.get(name, "") + line
    data = {}
    for d in ("fwd", "rev"):
        events = collections.defaultdict(list)
        for e in read_tsv(getattr(a, f"{d}_events")):
            if e["terminal"] == "0":
                events[(e["sample"], int(e["unique"]))].append(e)
        data[d] = dict(events=events, dada=by_sample(getattr(a, f"{d}_dada")),
                       derep=by_sample(getattr(a, f"{d}_derep")))

    stats, site, pairs_by_ref = collections.Counter(), collections.Counter(), collections.Counter()
    overlap_bases = 0
    for sample in sorted(set(data["fwd"]["dada"]) & set(data["rev"]["dada"])):
        D = {}
        for d in ("fwd", "rev"):
            dj = data[d]["dada"][sample]
            centres = [{n for n, t in refs.items()
                        if (t.startswith(c["sequence"]) if d == "fwd" else t.endswith(rc(c["sequence"])))}
                       for c in dj["asvs"]]
            D[d] = (data[d]["derep"][sample]["map"], dj["map"], centres)
        (fmap, fum, fc), (rmap, rum, rcs) = D["fwd"], D["rev"]
        for i in range(len(fmap)):
            uf, ur = fmap[i], rmap[i]
            cf, cr = fum[uf], rum[ur]
            if cf is None or cr is None:
                continue
            common = fc[cf] & rcs[cr]
            if len(common) != 1:
                continue
            ref_name = next(iter(common))
            ref = refs[ref_name]
            lo, hi = len(ref) - a.trunc_rev, a.trunc_fwd - 1
            if hi - lo < 10:
                continue
            pairs_by_ref[ref_name] += 1
            overlap_bases += hi - lo + 1
            keys = {}
            for d, u in (("fwd", uf), ("rev", ur)):
                keys[d] = {k for k in (event_key(e, d, ref) for e in data[d]["events"].get((sample, u), []))
                           if lo + 2 <= k[1] and k[2] <= hi - 2}
            for k in keys["fwd"] | keys["rev"]:
                cls = ("concordant" if k in keys["fwd"] and k in keys["rev"]
                       else "fwd only" if k in keys["fwd"] else "rev only")
                run = f"run {min(k[5], 5)}{'+' if k[5] >= 5 else ''}" if k[5] else "non-hp"
                stats[(f"{k[0]}{k[4] if k[4] < 3 else '3+'} {run}", cls)] += 1
                if cls == "concordant":
                    site[(ref_name, k)] += 1

    print(f"pairs with both reads on one reference: {sum(pairs_by_ref.values()):,}; overlap bases: {overlap_bases:,}")
    print(f"   {'event':14s} {'concordant':>10s} {'fwd only':>9s} {'rev only':>9s}  {'% concordant':>12s}")
    for kind in sorted({k for k, _ in stats}):
        c, f, r = (stats[(kind, x)] for x in ("concordant", "fwd only", "rev only"))
        print(f"   {kind:14s} {c:10d} {f:9d} {r:9d}  {100 * c / max(c + f + r, 1):11.1f}%")
    print("top concordant (molecular) sites, as % of the pairs from that reference:")
    for (ref_name, k), n in site.most_common(a.top):
        print(f"   {ref_name:36s} {k[0]}{k[4]} {k[3]} run {k[5]} @ {k[1]}-{k[2]}: "
              f"{n} pairs = {100 * n / pairs_by_ref[ref_name]:.2f}% of {pairs_by_ref[ref_name]:,}")


if __name__ == "__main__":
    main()
