#!/usr/bin/env python3
"""damage_pair_events.py — substitutions per read of each pair against a reference set (issue #291).

Usage:
  damage_pair_events.py REFS.fa OUT.tsv R1.fastq.gz R2.fastq.gz [R1 R2 ...]

Each pair is assigned to its closest reference by shared 12-mers; R1 and the
reverse-complemented R2 are aligned to it with free reference ends. Pairs with
more than 3 edits in either read (contaminants, chimeras, unlisted alleles) are
skipped. Everything is in forward reference coordinates, so template damage
copied by PCR shows the same substitution in R1 and R2, and a sequencing error
does not. The references are a mock's truth set cut to the amplicon, or a
sample's own abundant ASVs (on the PRJNA674596 Zymo mocks the two give the
same molecular rate).

Writes OUT.tsv (one row per substitution: pair, ref, pos, side, from, to, q),
OUT.tsv.pairs (each kept pair's reference spans per read), OUT.tsv.opp (bases
observed per ref, side, base and Q) and OUT.tsv.refs (reference index, names,
sequence). `damage_pair_spectrum.py` reads them. Needs Biopython.
"""
import collections, gzip, sys
from multiprocessing import Pool
from Bio.Align import PairwiseAligner

K = 12
RC = str.maketrans("ACGT", "TGCA")
_REFS = _IDX = _AL = None


def rc(s):
    return s.translate(RC)[::-1]


def init(refs):
    global _REFS, _IDX, _AL
    _REFS = refs
    _IDX = [{s[i:i + K] for i in range(len(s) - K + 1)} for _, s in refs]
    _AL = PairwiseAligner(mode="global", match_score=5, mismatch_score=-4,
                          open_gap_score=-9, extend_gap_score=-8,
                          end_deletion_score=0, end_insertion_score=-9)


def align(ref, seq, qual):
    """(events, opportunities, edits, span); events (pos, ref, alt, q), opportunities (ref base, q) counts."""
    a = _AL.align(ref, seq)[0]
    R, Q = a[0], a[1]
    # trim the free reference overhang
    lo = next(i for i, c in enumerate(Q) if c != "-")
    hi = len(Q) - next(i for i, c in enumerate(reversed(Q)) if c != "-")
    ri = sum(1 for c in R[:lo] if c != "-")
    start = ri
    qi, ev, opp, edits = 0, [], collections.Counter(), 0
    for x, y in zip(R[lo:hi], Q[lo:hi]):
        if x == "-":
            qi += 1; edits += 1; continue
        if y == "-":
            ri += 1; edits += 1; continue
        q = qual[qi]
        opp[(x, q)] += 1
        if x != y:
            ev.append((ri, x, y, q)); edits += 1
        ri += 1; qi += 1
    return ev, opp, edits, (start, ri)


def one(item):
    pid, s1, q1, s2, q2 = item
    ks = {s1[i:i + K] for i in range(0, len(s1) - K + 1, 4)}
    j = max(range(len(_REFS)), key=lambda i: len(ks & _IDX[i]))
    ref = _REFS[j][1]
    e1, o1, d1, sp1 = align(ref, s1, [ord(c) - 33 for c in q1])
    e2, o2, d2, sp2 = align(ref, rc(s2), [ord(c) - 33 for c in q2][::-1])
    if d1 > 3 or d2 > 3:
        return None
    return pid, j, e1, o1, e2, o2, (sp1, sp2)


def fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
            yield s, q


def main():
    truth, out, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    if len(files) < 2 or len(files) % 2:
        sys.exit(__doc__)
    seqs, name = {}, None
    for l in open(truth):
        l = l.strip()
        if l.startswith(">"):
            name = l[1:].split(";")[0]
        else:
            seqs[name] = seqs.get(name, "") + l
    uniq = collections.defaultdict(list)
    for n, s in seqs.items():
        uniq[s].append(n)
    refs = [("+".join(ns), s) for s, ns in uniq.items()]
    print(f"{len(seqs)} alleles, {len(refs)} distinct in the region", file=sys.stderr)

    pairs = []
    for f1, f2 in zip(files[::2], files[1::2]):
        pairs += [(len(pairs) + k, a[0], a[1], b[0], b[1]) for k, (a, b) in enumerate(zip(fastq(f1), fastq(f2)))]
    opp = collections.Counter()
    used = 0
    with Pool(8, initializer=init, initargs=(refs,)) as p, open(out, "w") as fh, open(out + ".pairs", "w") as ph:
        fh.write("pair\tref\tpos\tside\tfrom\tto\tq\n")
        ph.write("pair\tref\tr1_start\tr1_end\tr2_start\tr2_end\n")
        for r in p.imap_unordered(one, pairs, chunksize=200):
            if r is None:
                continue
            used += 1
            pid, j, e1, o1, e2, o2, (sp1, sp2) = r
            ph.write(f"{pid}\t{j}\t{sp1[0]}\t{sp1[1]}\t{sp2[0]}\t{sp2[1]}\n")
            for side, ev, o in (("R1", e1, o1), ("R2", e2, o2)):
                for (b, q), c in o.items():
                    opp[(j, side, b, q)] += c
                for pos, x, y, q in ev:
                    fh.write(f"{pid}\t{j}\t{pos}\t{side}\t{x}\t{y}\t{q}\n")
    with open(out + ".opp", "w") as fh:
        fh.write("ref\tside\tbase\tq\tcount\n")
        for (j, side, b, q), c in sorted(opp.items()):
            fh.write(f"{j}\t{side}\t{b}\t{q}\t{c}\n")
    with open(out + ".refs", "w") as fh:
        for j, (n, s) in enumerate(refs):
            fh.write(f"{j}\t{n}\t{s}\n")
    print(f"{used} of {len(pairs)} pairs within 3 edits per read", file=sys.stderr)


if __name__ == "__main__":
    main()
