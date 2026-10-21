#!/usr/bin/env python3
"""insilico_pcr.py — cut amplicon truth sets from genomes by in-silico PCR.

Builds a mock community's truth set for a given primer pair (issue #166: a
truth set is specific to the mock, the region AND the primers). Each genome
FASTA is one organism. Both strands are scanned for the forward primer and, in
amplicon orientation, the reverse-complemented reverse primer, each within
--max-mismatch (IUPAC-aware); a forward and reverse hit on the same strand
within --min-len..--max-len make one copy; each reverse site keeps its
best-matching forward site, so a nested off-target site does not add a second,
longer product. Primers are excluded from the
amplicon, as they are trimmed from reads.

Primer hits are found from an exact match of each primer's 3' --seed bases,
the end that must anneal for extension, then scored over the full primer.

Outputs:
  --fasta   alleles: identical copies within an organism merged, `copies=N`;
            a sequence shared between organisms is reported once per organism
  --tsv     every copy: organism, contig, strand, position, primer mismatches,
            length, allele

Usage:
  insilico_pcr.py --fwd CCTACGGGNGGCWGCAG --rev GACTACHVGGGTATCTAATCC \\
      --min-len 300 --max-len 600 --fasta truth.fa --tsv copies.tsv genomes/*.fna
"""

import argparse
import collections
import gzip
import os
import re
import sys

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "M": "AC", "K": "GT",
    "S": "CG", "W": "AT", "N": "ACGT", "H": "ACT", "V": "ACG", "B": "CGT", "D": "AGT",
}
COMP = str.maketrans("ACGTRYMKSWNHVBD", "TGCAYRKMSWNDBVH")


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name:
        yield name, "".join(buf).upper()


def mismatches(primer, site):
    return sum(c not in IUPAC.get(p, p) for p, c in zip(primer, site))


def hits(primer, seq, seed, max_mm, seed_at_end):
    """Start positions where `primer` matches `seq` within max_mm, seeded by an
    exact match of its 3' `seed` bases (`seed_at_end`) or of its first `seed`
    bases (for a reverse-complemented reverse primer, whose 3' end is first)."""
    k = len(primer)
    core = primer[-seed:] if seed_at_end else primer[:seed]
    rx = re.compile("(?=" + "".join(f"[{IUPAC[c]}]" for c in core) + ")")
    out = []
    for m in rx.finditer(seq):
        start = m.start() - (k - seed) if seed_at_end else m.start()
        if start < 0 or start + k > len(seq):
            continue
        mm = mismatches(primer, seq[start:start + k])
        if mm <= max_mm:
            out.append((start, mm))
    return out


def organism(path, first_header):
    words = first_header.split(None, 1)[1].split(",")[0].split() if " " in first_header else [os.path.basename(path)]
    return "_".join(w for w in words[:2]) if len(words) >= 2 else words[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("genomes", nargs="+")
    ap.add_argument("--fwd", required=True, help="forward primer, 5'->3'")
    ap.add_argument("--rev", required=True, help="reverse primer, 5'->3' (catalog orientation)")
    ap.add_argument("--max-mismatch", type=int, default=3)
    ap.add_argument("--seed", type=int, default=5, help="exact 3' bases required (default 5)")
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--max-len", type=int, default=2000)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tsv", required=True)
    a = ap.parse_args()
    fwd, rrc = a.fwd.upper(), rc(a.rev.upper())

    copies = []
    for path in a.genomes:
        org = None
        for name, g in fasta(path):
            org = org or organism(path, name)
            contig = name.split()[0]
            for strand, s in (("+", g), ("-", rc(g))):
                fh = hits(fwd, s, a.seed, a.max_mismatch, True)
                rh = hits(rrc, s, a.seed, a.max_mismatch, False)
                # One product per reverse site: its best-matching forward site,
                # the shortest on a tie (a nested off-target site primes poorly).
                best = {}
                for rs, rmm in rh:
                    for fs, fmm in fh:
                        size = rs - (fs + len(fwd))
                        if a.min_len <= size <= a.max_len:
                            key = (fmm, size)
                            if rs not in best or key < best[rs][0]:
                                best[rs] = (key, fs, fmm, rmm)
                for rs, (_, fs, fmm, rmm) in sorted(best.items()):
                    pos = fs if strand == "+" else len(g) - (rs + len(rrc))
                    copies.append((org, contig, strand, pos, fmm, rmm, s[fs + len(fwd):rs]))

    alleles = collections.OrderedDict()
    for org, *_, amp in copies:
        alleles.setdefault((org, amp), 0)
        alleles[(org, amp)] += 1
    names, per = {}, collections.Counter()
    with open(a.fasta, "w") as fh:
        for (org, amp), n in alleles.items():
            per[org] += 1
            names[(org, amp)] = f"{org}_a{per[org]}"
            fh.write(f">{org}_a{per[org]};copies={n};len={len(amp)}\n{amp}\n")
    with open(a.tsv, "w") as fh:
        fh.write("organism\tcontig\tstrand\tposition\tfwd_mismatches\trev_mismatches\tlength\tallele\n")
        for org, contig, strand, pos, fmm, rmm, amp in copies:
            fh.write(f"{org}\t{contig}\t{strand}\t{pos}\t{fmm}\t{rmm}\t{len(amp)}\t{names[(org, amp)]}\n")

    ncopy = collections.Counter(c[0] for c in copies)
    seen = {organism(p, next(fasta(p))[0]) for p in a.genomes}
    for org in sorted(seen):
        mm = sorted({(c[4], c[5]) for c in copies if c[0] == org})
        print(f"{org:34s} copies={ncopy[org]:3d} alleles={per[org]:2d} primer mismatches F/R={mm}", file=sys.stderr)
    shared = collections.Counter(amp for (org, amp) in alleles)
    print(f"{len(copies)} copies, {len(alleles)} organism alleles, {len(shared)} distinct sequences "
          f"({sum(1 for v in shared.values() if v > 1)} shared between organisms)", file=sys.stderr)


if __name__ == "__main__":
    main()
