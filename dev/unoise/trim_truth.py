#!/usr/bin/env python3
"""Trim a reference FASTA to the region between two primers (#323).

In-silico PCR for a mock-community truth set: for each reference, find the
forward primer and, after it, the reverse complement of the reverse primer,
and keep the sequence between them, primers excluded, as in primer-trimmed
reads. Identical regions collapse into one record whose header joins the
source labels with `|` (e.g. S. aureus and S. epidermidis share a V4).

Primers are matched IUPAC-aware with up to --max-mismatch mismatches: in a
genome the primer site is genomic, not the synthetic primer that PCR writes
into a read, so it can differ (HMP mock P. acnes has one per site at 515F/806R).

Defaults are 515F/806R (V4), the MiSeq SOP primers:
  dev/unoise/trim_truth.py HMP_MOCK.v35.fasta > hmp_mock_v4.fasta

The report (stderr) lists references without both primers, the region
lengths, and which references collapsed.
"""
import argparse
import sys
from collections import Counter

AMB = {"M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
       "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "N": "ACGT"}
COMP = str.maketrans("ACGTMRWSYKVHDBN", "TGCAKYWSRMBDHVN")


def rc(s):
    return s[::-1].translate(COMP)


def find(seq, primer, start, max_mm):
    """Leftmost best-scoring match at or after start, or -1."""
    best, pos = max_mm + 1, -1
    for i in range(start, len(seq) - len(primer) + 1):
        mm = 0
        for a, b in zip(seq[i:i + len(primer)], primer):
            if a != b and a not in AMB.get(b, ""):
                mm += 1
                if mm >= best:
                    break
        if mm < best:
            best, pos = mm, i
    return pos


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("reference")
    ap.add_argument("--fwd", default="GTGCCAGCMGCCGCGGTAA", help="forward primer (515F)")
    ap.add_argument("--rev", default="GGACTACHVGGGTWTCTAAT", help="reverse primer (806R)")
    ap.add_argument("--max-mismatch", type=int, default=2)
    a = ap.parse_args()
    fwd, rev = a.fwd.upper(), rc(a.rev.upper())

    recs = []
    with open(a.reference) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                recs.append([line[1:].split()[0], []])
            elif line:
                recs[-1][1].append(line.upper())
    regions, missing = {}, []
    for label, parts in recs:
        s = "".join(parts)
        i = find(s, fwd, 0, a.max_mismatch)
        j = find(s, rev, i + len(fwd), a.max_mismatch) if i >= 0 else -1
        if j < 0:
            missing.append(label)
            continue
        regions.setdefault(s[i + len(fwd):j], []).append(label)

    print(f"{len(recs)} references, {len(regions)} distinct regions", file=sys.stderr)
    if missing:
        print(f"without both primers: {', '.join(missing)}", file=sys.stderr)
    print("lengths:", dict(sorted(Counter(map(len, regions)).items())), file=sys.stderr)
    for seq, labels in regions.items():
        if len(labels) > 1:
            print(f"collapsed: {', '.join(labels)}", file=sys.stderr)
        print(f">{'|'.join(labels)}\n{seq}")


if __name__ == "__main__":
    main()
