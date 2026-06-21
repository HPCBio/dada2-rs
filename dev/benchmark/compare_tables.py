#!/usr/bin/env python3
"""Compare final sequence tables (seqtab / seqtab_nochim JSON) across runs.

A standalone version of the ASV-table equivalence that bench_pooled.py reports
in its WFA sweep CSV (`jaccard_vs_nw`, `count_frac_vs_nw`). Use it to compare any
two+ runs' tables directly — e.g. nw vs wfa arms, k5 vs k7, two pool modes —
without re-running a sweep.

The `load_table` / `table_equiv` definitions are kept byte-for-byte identical to
bench_pooled.py so the numbers here match that CSV exactly. (compare_asvs.py is
the richer tool — abundance-stratified churn + nearest-neighbour Hamming — when
you need to characterize *what* diverged; this one is the quick jaccard/count
equivalence check.)

Metrics (each comparison vs the baseline):
  jaccard     — |shared ASVs| / |union ASVs| (sequence sets; 1.0 = same ASV set)
  count_frac  — of the SHARED ASVs, fraction whose full per-sample count vector
                is identical (the CSV's count_frac_vs_nw)
  identical   — True iff same ASV set AND every per-sample count matches
  +only / -only — ASVs the comparison has that the baseline lacks / vice-versa

Usage:
  compare_tables.py BASELINE=path/seqtab_nochim.json LABEL=path ... [--abundance N]

  --abundance N  also report jaccard restricted to ASVs with total count >= N
                 (mirrors the MIN_ABUNDANCE filter; default 1 = no filter)

Example:
  compare_tables.py \\
      nw=tmp/pb_k7/nw/seqtab_nochim.json \\
      e30=tmp/pb_k7/e30/seqtab_nochim.json \\
      e0=tmp/pb_k7/e0/seqtab_nochim.json
"""
import json
import sys


def load_table(nochim_json):
    """Load a make-sequence-table / remove-bimera-denovo JSON into a dict
    {sequence -> {sample -> count}} (per-sample, order-independent). Returns None
    if the file is missing/unreadable."""
    try:
        with open(nochim_json) as fh:
            d = json.load(fh)
        samples, seqs, counts = d["samples"], d["sequences"], d["counts"]
    except (OSError, ValueError, KeyError):
        return None
    # counts[i][j] = sample i, sequence j
    table = {}
    for j, seq in enumerate(seqs):
        table[seq] = {samples[i]: counts[i][j] for i in range(len(samples))}
    return table


def table_equiv(a, b):
    """Compare two {seq -> {sample -> count}} tables. Returns
    (jaccard, exact_frac, identical) where jaccard is over the ASV sequence sets,
    exact_frac is the fraction of shared ASVs whose full per-sample count vector
    matches, and identical is True iff the tables are byte-equal (same ASV set
    AND every count). None-safe: returns (None, None, None) if either is None."""
    if a is None or b is None:
        return (None, None, None)
    sa, sb = set(a), set(b)
    union = sa | sb
    shared = sa & sb
    jaccard = 1.0 if not union else len(shared) / len(union)
    exact = sum(1 for s in shared if a[s] == b[s])
    exact_frac = 1.0 if not shared else exact / len(shared)
    identical = (sa == sb) and all(a[s] == b[s] for s in shared)
    return (jaccard, exact_frac, identical)


def total_count(table):
    """ASV -> summed-over-samples count, for abundance filtering."""
    return {seq: sum(per.values()) for seq, per in table.items()}


def filter_abundance(table, min_abund):
    if min_abund <= 1:
        return table
    tot = total_count(table)
    return {s: table[s] for s in table if tot[s] >= min_abund}


def main(argv):
    # Parse --abundance[=N | N] first, then treat the rest as LABEL=path specs.
    args, min_abund = [], 1
    it = iter(argv[1:])
    for a in it:
        if a.startswith("--abundance"):
            min_abund = int(a.split("=", 1)[1]) if "=" in a else int(next(it))
        else:
            args.append(a)
    specs = [a for a in args if "=" in a]
    if len(specs) < 2:
        sys.exit(f"usage: {argv[0]} BASELINE=path LABEL=path [...] [--abundance N]")

    base_label, base_path = specs[0].split("=", 1)
    base = load_table(base_path)
    if base is None:
        sys.exit(f"baseline {base_label} ({base_path}) unreadable")
    base_f = filter_abundance(base, min_abund)

    width = max(len(s.split("=", 1)[0]) for s in specs)
    af = f" (abundance>={min_abund})" if min_abund > 1 else ""
    print(f"baseline: {base_label} = {base_path}  ({len(base_f)} ASVs{af})\n")
    hdr = (f"{'label':>{width}}  {'ASVs':>6} {'jaccard':>8} {'count_frac':>11} "
           f"{'+only':>6} {'-only':>6} {'identical':>9}")
    print(hdr)
    print("-" * len(hdr))
    for spec in specs[1:]:
        label, path = spec.split("=", 1)
        t = load_table(path)
        if t is None:
            print(f"{label:>{width}}  <unreadable: {path}>")
            continue
        tf = filter_abundance(t, min_abund)
        jac, cnt, ident = table_equiv(base_f, tf)
        plus = len(set(tf) - set(base_f))   # ASVs in comparison not in baseline
        minus = len(set(base_f) - set(tf))  # ASVs in baseline not in comparison
        print(f"{label:>{width}}  {len(tf):>6} {jac:>8.4f} {cnt:>11.4f} "
              f"{plus:>6} {minus:>6} {str(ident):>9}")


if __name__ == "__main__":
    main(sys.argv)
