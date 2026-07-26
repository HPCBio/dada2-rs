#!/usr/bin/env python3
"""Join the Rust arms against the R ground truth and report divergence.

The metric is |log10(rate_arm / rate_R)| on the CLAMPED ERROR RATE, not on the
raw log10 fit. That is the number that reaches inference: 0.30 means the arm
disagrees with R by a factor of two on an error rate, whereas a large
disagreement in raw fit space may vanish once `[1e-7, 0.25]` clamping is
applied (or vice versa).

Each arm is compared against the R surface it is configured to reproduce:
`*_interp` arms against R `surface="interpolate"`, everything else against
`surface="direct"`.

  --by-q       divergence profile per quality column, to localise it (boundary
               vertices vs interior vs high-Q tail)
  --signed     mean signed log10 ratio per q: does the arm push rates up or down?
  --floor      count cells pinned at the 1e-7 floor, i.e. fits that failed

Usage: compare.py <r_pred.tsv> <rust_pred.tsv> [--by-q] [--signed] [--floor]
"""

import argparse
import csv
import math
import statistics as st

MIN_RATE = 1e-7

# Which R surface each arm is trying to reproduce.
REF = {
    "ours_direct": "direct",
    "ours_interp": "interpolate",
    "lrs_direct": "direct",
    "lrs_interp": "interpolate",
    "lrs_robust_mad": "direct",
    "lrs_reflect": "direct",
    "lrs_extend": "direct",
}


def load(path, keycols):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                val = float(r["rate"])
            except (ValueError, KeyError):
                val = float("nan")
            out[tuple(r[c] for c in keycols)] = val
    return out


def ratios(R, S, arm, surf, pairs, qs):
    """|log10(arm/R)| over all (pair, q), skipping cells either side lost."""
    out = []
    for p in pairs:
        for q in qs:
            sv = S.get((p, q, arm))
            rv = R.get((p, q, surf))
            if sv is None or rv is None:
                continue
            if not (math.isfinite(sv) and math.isfinite(rv) and sv > 0 and rv > 0):
                continue
            out.append((p, q, math.log10(sv / rv)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r_pred")
    ap.add_argument("rust_pred")
    ap.add_argument("--by-q", action="store_true")
    ap.add_argument("--signed", action="store_true")
    ap.add_argument("--floor", action="store_true")
    args = ap.parse_args()

    R = load(args.r_pred, ["pair", "q", "surface"])
    S = load(args.rust_pred, ["pair", "q", "arm"])

    pairs = sorted({k[0] for k in S})
    qs = sorted({k[1] for k in S}, key=float)
    arms = [a for a in REF if any(k[2] == a for k in S)]

    print(f"{len(pairs)} transitions x {len(qs)} quality columns\n")
    hdr = f"{'arm':<16}{'vs R':<13}{'median':>11}{'p95':>11}{'max':>11}{'n>2x':>7}{'cells':>8}"
    print(hdr)
    print("-" * len(hdr))
    for arm in arms:
        surf = REF[arm]
        d = [abs(v) for _, _, v in ratios(R, S, arm, surf, pairs, qs)]
        if not d:
            print(f"{arm:<16}{surf:<13}{'no comparable cells':>48}")
            continue
        d.sort()
        big = sum(1 for v in d if v > math.log10(2))
        print(
            f"{arm:<16}{surf:<13}{st.median(d):>11.3e}"
            f"{d[int(0.95 * (len(d) - 1))]:>11.3e}{d[-1]:>11.3e}{big:>7}{len(d):>8}"
        )

    if args.floor:
        print("\ncells pinned at the 1e-7 floor (a failed fit):")
        for arm in arms:
            n = sum(1 for k in S if k[2] == arm and S[k] == MIN_RATE)
            tot = sum(1 for k in S if k[2] == arm)
            flag = "  <-- ALL FLOORED" if n == tot and tot else ""
            print(f"  {arm:<16}{n:>6}/{tot}{flag}")

    if args.by_q:
        print("\nmax |log10 ratio| by quality column:")
        print(f"{'q':>5}" + "".join(f"{a:>16}" for a in arms))
        for q in qs:
            row = f"{q:>5}"
            for arm in arms:
                d = [abs(v) for _, qq, v in ratios(R, S, arm, REF[arm], pairs, [q])]
                row += f"{max(d):>16.3e}" if d else f"{'-':>16}"
            print(row)

    if args.signed:
        print("\nmean SIGNED log10(arm/R) by quality column (+ = arm rate higher):")
        print(f"{'q':>5}" + "".join(f"{a:>16}" for a in arms))
        for q in qs:
            row = f"{q:>5}"
            for arm in arms:
                d = [v for _, qq, v in ratios(R, S, arm, REF[arm], pairs, [q])]
                row += f"{st.mean(d):>+16.3f}" if d else f"{'-':>16}"
            print(row)


if __name__ == "__main__":
    main()
