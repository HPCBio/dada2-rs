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

Usage:
  analyze_kdist_curves.py LABEL=curve.csv [LABEL=curve.csv ...] [--target-div 10]

Example:
  dada2-rs kdist-calibrate derep/*.json.gz --k 7 --per-sample \\
      --max-pairs 500000 --threads 32 -o kmer.csv
  dada2-rs kdist-calibrate derep/*.json.gz --k 7 --per-sample \\
      --screen-backend minimizer --minimizer-k 8 \\
      --max-pairs 500000 --threads 32 -o mini.csv
  analyze_kdist_curves.py kmer=kmer.csv minimizer=mini.csv
"""

import argparse
import csv
import sys

# Cutoffs to report retention at. Spans the DADA2 default through the range a
# sketch-based screen needs.
CUTOFFS = [0.30, 0.42, 0.50, 0.55, 0.60, 0.65, 0.70, 0.72, 0.75, 0.80]
DIV_BANDS = [1, 2, 3, 5, 10]


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                rows.append((float(r["kdist"]), float(r["pct_div"])))
            except (KeyError, ValueError):
                continue
    if not rows:
        sys.exit(f"{path}: no usable rows (expected kdist,pct_div columns)")
    return rows


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
        data[label] = load(path)

    n = {k: len(v) for k, v in data.items()}
    if len(set(n.values())) > 1:
        print(f"WARNING: curves have different pair counts {n}.", file=sys.stderr)
        print(
            "  They must come from the SAME input and seed, or the comparison is\n"
            "  between different pair samples rather than between screens.",
            file=sys.stderr,
        )

    tgt = args.target_div
    print(f"pairs per curve: {n}")
    print()
    print(f"Max screen distance among pairs at or below a given TRUE divergence.")
    print(f"DADA2's 0.42 sits just above the k-mer row at {tgt:g}% -- that IS the")
    print(f"documented '0.42 ~ {tgt:g}% divergence' calibration. Read each other")
    print("screen's equivalent cutoff off the same row.")
    print()
    hdr = f"{'true div <=':>12s} {'n':>8s}" + "".join(f"{l:>16s}" for l in data)
    print(hdr)
    for band in DIV_BANDS:
        counts = [sum(1 for d, pd in v if pd <= band) for v in data.values()]
        cells = ""
        for v in data.values():
            sel = [d for d, pd in v if pd <= band]
            cells += f"{max(sel):16.4f}" if sel else f"{'-':>16s}"
        print(f"{band:11d}% {counts[0]:8d}{cells}")

    print()
    print(f"RECALL: fraction of pairs at <= {tgt:g}% true divergence that the screen")
    print("still passes, and the share of ALL pairs it passes (the cost side).")
    print()
    print(f"{'cutoff':>8s}" + "".join(f"{l + ' recall':>18s}{l + ' pass':>16s}" for l in data))
    for c in CUTOFFS:
        line = f"{c:8.2f}"
        for v in data.values():
            close = [d for d, pd in v if pd <= tgt]
            kept = sum(1 for d in close if d <= c)
            passed = sum(1 for d, _ in v if d <= c)
            line += f"{100 * kept / max(len(close), 1):17.2f}%{100 * passed / len(v):15.2f}%"
        print(line)

    print()
    print("Recommended cutoff per screen = the smallest listed cutoff reaching")
    print(f"100% recall at {tgt:g}% divergence (i.e. lossless on this data):")
    for label, v in data.items():
        close = [d for d, pd in v if pd <= tgt]
        need = max(close) if close else float("nan")
        rec = next((c for c in CUTOFFS if c >= need), None)
        passed = sum(1 for d, _ in v if d <= rec) / len(v) * 100 if rec else float("nan")
        print(
            f"  {label:>12s}: >= {need:.4f}  -> use {rec}  (passes {passed:.1f}% of all pairs)"
        )
    print()
    print("Cross-check the choice against ALIGNMENT COUNT before adopting it:")
    print("b_compare is align-dominated, so a cutoff that passes materially more")
    print("pairs than the k-mer screen makes the run slower even when it is more")
    print("accurate. Buying the last fraction of a percent of recall is expensive.")


if __name__ == "__main__":
    main()
