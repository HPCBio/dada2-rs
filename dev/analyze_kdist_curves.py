#!/usr/bin/env python3
"""analyze_kdist_curves.py — derive a screen's cutoff from its calibration curve.

Takes `kdist-calibrate` CSVs (one per screen backend, produced on the SAME input
so the sampled pairs and their true divergences match) and answers two questions
that a raw curve does not:

  1. RECALL: at a given cutoff, what fraction of genuine near-neighbour pairs
     does the screen still let through? This is the quantity that decides
     correctness. A screen that drops real neighbours starves raws of a cluster
     match and they birth spurious low-abundance ASVs.

  2. THE EQUIVALENT CUTOFF: DADA2's KDIST_CUTOFF = 0.42 is documented as
     "~10% nucleotide divergence". Empirically that is exactly what it is -- the
     largest k-mer distance among pairs at <=10% true divergence lands just under
     0.42. Applying the SAME RULE to another screen gives its equivalent cutoff.
     This is the only defensible way to port a cutoff: a cutoff is a property of
     the distance DISTRIBUTION a metric induces, not of its formula, so two
     screens sharing an algebraic form still need different numbers.

MEMORY: these CSVs are one row per sampled pair and reach tens of GB on a large
run. Everything here is computed from running tallies in a SINGLE STREAMING PASS
per file, with O(#cutoffs) memory -- the rows are never retained. Files may be
plain or gzipped.

Usage:
  analyze_kdist_curves.py LABEL=curve.csv[.gz] [LABEL=curve.csv[.gz] ...]

Example:
  dada2-rs kdist-calibrate derep/*.json.gz --k 7 --per-sample \\
      --max-pairs 500000 --threads 32 -o kmer.csv
  dada2-rs kdist-calibrate derep/*.json.gz --k 7 --per-sample \\
      --screen-backend minimizer --minimizer-k 8 \\
      --max-pairs 500000 --threads 32 -o mini.csv
  analyze_kdist_curves.py kmer=kmer.csv minimizer=mini.csv
"""

import argparse
import bisect
import csv
import gzip
import sys

# Cutoffs to report retention at. Spans the DADA2 default through the range a
# sketch-based screen needs on either platform.
CUTOFFS = [0.30, 0.40, 0.42, 0.45, 0.48, 0.50, 0.55, 0.60, 0.65, 0.70, 0.72, 0.75, 0.80]
DIV_BANDS = [1, 2, 3, 5, 10]


def scan(path, target_div):
    """One streaming pass. Returns tallies only -- never the rows.

    Both `DIV_BANDS` and `CUTOFFS` are sorted, so a row's contribution is a
    SUFFIX of each list ("passes every cutoff at or above its distance"). Rather
    than testing all 18 thresholds per row -- ~1.7e9 comparisons on a 94M-row
    file -- each row does two bisects into per-bucket tallies, and the suffix
    sums are accumulated once at the end.
    """
    opener = gzip.open if path.endswith(".gz") else open
    nb, nc = len(DIV_BANDS), len(CUTOFFS)
    # bucket i = rows whose pd first qualifies at DIV_BANDS[i]
    band_bucket_n = [0] * (nb + 1)
    band_bucket_max = [None] * (nb + 1)
    # bucket i = rows whose d first passes at CUTOFFS[i]
    cut_bucket_n = [0] * (nc + 1)
    cut_bucket_close = [0] * (nc + 1)
    n = close_n = 0
    close_max = None

    bl, cl = DIV_BANDS, CUTOFFS
    with opener(path, "rt", newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr, None)
        if not header:
            sys.exit(f"{path}: empty")
        try:
            ki = header.index("kdist")
            pi = header.index("pct_div")
        except ValueError:
            sys.exit(f"{path}: expected kdist,pct_div columns; got {header}")
        for row in rdr:
            try:
                d = float(row[ki]); pd = float(row[pi])
            except (IndexError, ValueError):
                continue
            n += 1
            bi = bisect.bisect_left(bl, pd)
            band_bucket_n[bi] += 1
            m = band_bucket_max[bi]
            if m is None or d > m:
                band_bucket_max[bi] = d
            ci = bisect.bisect_left(cl, d)
            cut_bucket_n[ci] += 1
            if pd <= target_div:
                close_n += 1
                cut_bucket_close[ci] += 1
                if close_max is None or d > close_max:
                    close_max = d
    if n == 0:
        sys.exit(f"{path}: no usable rows")

    # Prefix-accumulate: band b covers every bucket at or below it.
    band_max, band_n, run_n, run_max = {}, {}, 0, None
    for i, b in enumerate(bl):
        run_n += band_bucket_n[i]
        m = band_bucket_max[i]
        if m is not None and (run_max is None or m > run_max):
            run_max = m
        band_n[b], band_max[b] = run_n, run_max
    # Same for cutoffs: cutoff c passes every bucket at or below it.
    pass_n, close_pass, rn, rc = {}, {}, 0, 0
    for i, c in enumerate(cl):
        rn += cut_bucket_n[i]
        rc += cut_bucket_close[i]
        pass_n[c], close_pass[c] = rn, rc

    return dict(n=n, band_max=band_max, band_n=band_n, pass_n=pass_n,
                close_pass=close_pass, close_n=close_n, close_max=close_max)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("curves", nargs="+", metavar="LABEL=CSV")
    ap.add_argument(
        "--target-div",
        type=float,
        default=10.0,
        help="divergence (%%) defining a 'genuine near neighbour' (default 10, "
        "which is what DADA2's 0.42 was calibrated to)",
    )
    args = ap.parse_args()

    data = {}
    for spec in args.curves:
        label, _, path = spec.partition("=")
        print(f"scanning {label} ({path}) ...", file=sys.stderr)
        data[label] = scan(path, args.target_div)

    counts = {k: v["n"] for k, v in data.items()}
    if len(set(counts.values())) > 1:
        print(f"WARNING: curves have different pair counts {counts}.", file=sys.stderr)
        print(
            "  They must come from the SAME input and seed, or the comparison is\n"
            "  between different pair samples rather than between screens.",
            file=sys.stderr,
        )

    tgt = args.target_div
    print(f"pairs per curve: {counts}")
    print()
    print("Max screen distance among pairs at or below a given TRUE divergence.")
    print(f"DADA2's 0.42 sits just above the k-mer row at {tgt:g}% -- that IS the")
    print(f"documented '0.42 ~ {tgt:g}% divergence' calibration. Read each other")
    print("screen's equivalent cutoff off the same row.")
    print()
    print(f"{'true div <=':>12s} {'n':>10s}" + "".join(f"{l:>16s}" for l in data))
    for band in DIV_BANDS:
        first = next(iter(data.values()))
        cells = ""
        for v in data.values():
            m = v["band_max"][band]
            cells += f"{m:16.4f}" if m is not None else f"{'-':>16s}"
        print(f"{band:11d}% {first['band_n'][band]:10d}{cells}")

    print()
    print(f"RECALL: fraction of pairs at <= {tgt:g}% true divergence that the screen")
    print("still passes, and the share of ALL pairs it passes (the cost side).")
    print()
    print(f"{'cutoff':>8s}" + "".join(f"{l + ' recall':>18s}{l + ' pass':>16s}" for l in data))
    for c in CUTOFFS:
        line = f"{c:8.2f}"
        for v in data.values():
            r = 100 * v["close_pass"][c] / max(v["close_n"], 1)
            p = 100 * v["pass_n"][c] / v["n"]
            line += f"{r:17.2f}%{p:15.2f}%"
        print(line)

    print()
    print("Recommended cutoff per screen = the smallest listed cutoff reaching")
    print(f"100% recall at {tgt:g}% divergence (i.e. lossless on this data):")
    for label, v in data.items():
        need = v["close_max"]
        rec = next((c for c in CUTOFFS if need is not None and c >= need), None)
        p = 100 * v["pass_n"][rec] / v["n"] if rec else float("nan")
        print(f"  {label:>12s}: >= {need:.4f}  -> use {rec}  (passes {p:.1f}% of all pairs)")
    print()
    print("DO NOT ADOPT THAT NUMBER WITHOUT CHECKING ALIGNMENT COUNT.")
    print("b_compare is 93-97% align-dominated, so the cutoff sets the cost, and")
    print("100% recall is the expensive end of the knob. On PacBio it picked 0.50,")
    print("which aligns ~10% MORE pairs than the k-mer screen and runs SLOWER,")
    print("while 0.42-0.45 aligns 7-11% FEWER for a bit-identical ASV set -- the")
    print("last 1.5% of recall cost ~20% of the alignment work and bought nothing.")
    print("Use this to BRACKET the region, then let a sweep choose within it.")


if __name__ == "__main__":
    main()
