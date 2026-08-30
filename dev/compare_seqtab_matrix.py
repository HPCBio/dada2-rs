#!/usr/bin/env python3
"""compare_seqtab_matrix.py — exact sample x ASV count-matrix equality.

`compare_asvs.py` answers the *set* question (which ASVs churned) and summarises
pooled abundance. This answers the stricter one the concordance standard actually
implies: is the full **sample x ASV count matrix** identical between two runs?

Order-insensitive on both axes, because ASV ordering in `seqtab.nochim.json` is
NOT stable run-to-run -- two runs of the same binary emit the same ASVs and the
same counts in a different order (observed on the default k-mer backend). A
byte-diff or a naive list comparison therefore reports spurious differences; the
matrix is the invariant, not its layout.

Usage:
  compare_seqtab_matrix.py <a.json> <b.json> [--label-a A] [--label-b B]

Exit status is 0 when the matrices are identical, 1 otherwise, so it can gate CI.
"""

import argparse
import json
import sys


def load(path):
    """Return {sample: {sequence: count}} from a seqtab JSON."""
    j = json.load(open(path))
    samples, seqs, counts = j["samples"], j["sequences"], j["counts"]
    if len(counts) != len(samples):
        sys.exit(f"{path}: {len(counts)} count rows for {len(samples)} samples")
    out = {}
    for s, row in zip(samples, counts):
        if len(row) != len(seqs):
            sys.exit(f"{path}: sample {s} has {len(row)} counts for {len(seqs)} ASVs")
        # Drop zeros so a run that merely *carries* an ASV column at zero is not
        # counted as disagreeing with one that never emitted the column.
        out[s] = {q: c for q, c in zip(seqs, row) if c}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    args = ap.parse_args()

    A, B = load(args.a), load(args.b)
    la, lb = args.label_a, args.label_b

    ok = True
    if set(A) != set(B):
        ok = False
        print(f"SAMPLE SETS DIFFER: only {la}={sorted(set(A)-set(B))} "
              f"only {lb}={sorted(set(B)-set(A))}")

    shared_samples = sorted(set(A) & set(B))
    asvs_a = {q for s in A.values() for q in s}
    asvs_b = {q for s in B.values() for q in s}
    print(f"{la}: {len(A)} samples, {len(asvs_a)} ASVs, {sum(sum(s.values()) for s in A.values())} reads")
    print(f"{lb}: {len(B)} samples, {len(asvs_b)} ASVs, {sum(sum(s.values()) for s in B.values())} reads")
    print(f"ASV set identical: {asvs_a == asvs_b}"
          + ("" if asvs_a == asvs_b
             else f"  (only_{la}={len(asvs_a-asvs_b)} only_{lb}={len(asvs_b-asvs_a)})"))
    if asvs_a != asvs_b:
        ok = False

    # The count axis, cell by cell.
    diff_cells = 0
    diff_abs = 0
    total_cells = 0
    worst = []
    for s in shared_samples:
        for q in set(A[s]) | set(B[s]):
            total_cells += 1
            ca, cb = A[s].get(q, 0), B[s].get(q, 0)
            if ca != cb:
                diff_cells += 1
                diff_abs += abs(ca - cb)
                worst.append((abs(ca - cb), s, ca, cb))
    reads = sum(sum(s.values()) for s in A.values())
    print(f"count matrix: {diff_cells} of {total_cells} non-zero cells differ; "
          f"L1 = {diff_abs} reads ({100.0 * diff_abs / max(reads, 1):.4f}% of {la} total)")
    if diff_cells:
        ok = False
        worst.sort(reverse=True)
        print(f"  largest {min(5, len(worst))} disagreements (|d|, sample, {la}, {lb}):")
        for d, s, ca, cb in worst[:5]:
            print(f"    {d:>8}  {s}  {ca} vs {cb}")

    print("==> IDENTICAL at the ASV + count level" if ok else "==> DIFFERS")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
