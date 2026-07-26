#!/usr/bin/env python3
"""Extract the exact (q, errs, tot, rlogp) triples that loess_errfun feeds its
smoother, from a learn-errors JSON, into a TSV the other arms can all read.

This is deliberately a re-derivation of loess_errfun's pre-smoothing arithmetic
(src/error_models.rs) rather than a dump of its output: the whole point is to
hand every implementation byte-identical inputs.

    rlogp = log10((errs + 1) / tot)      # NaN where tot == 0
    tot   = colSums over the 4 destinations from this source nucleotide

Optional --bins re-derives the triples as if the data were binned, by
collapsing every quality column onto its nearest bin value. Two shapes:

  --bins 12,24,38            full 40-column grid, counts only at those columns
                             (the production shape: what learn-errors actually
                             sees on binned input)
  --bins 12,24,38 --anchors-only
                             just the bin columns, nothing else
                             (isolates n_valid, the degenerate-neighborhood axis)

Usage:
  extract_triples.py <learn_errors.json> <out.tsv> [--bins 12,24,38] [--anchors-only]
  extract_triples.py <learn_errors.json> <out.tsv> --nbins 5
"""

import argparse
import gzip
import json
import math
import sys

NT = "ACGT"


def load_trans(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        d = json.load(fh)
    if "trans" not in d:
        sys.exit(
            f"{path} has no `trans` block. learn-errors JSONs written by "
            "learnerrors_to_dada2rs.R do not carry one; use a native "
            "dada2-rs learn-errors output."
        )
    trans, nq = d["trans"], d["nq"]
    if len(trans) != 16:
        sys.exit(f"expected 16 transition rows, got {len(trans)}")
    return {NT[i] + "2" + NT[j]: trans[i * 4 + j] for i in range(4) for j in range(4)}, nq


def rebin(rows, nq, bins, anchors_only):
    """Collapse counts onto the nearest bin value.

    Empty columns are skipped rather than mapped, so a column with no
    observations cannot invent counts at a bin.
    """
    grid = bins if anchors_only else list(range(nq))
    out = {k: {q: 0 for q in grid} for k in rows}
    for k, counts in rows.items():
        for q, c in enumerate(counts):
            if c:
                out[k][min(bins, key=lambda b: abs(b - q))] += c
    return out, grid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json")
    ap.add_argument("out")
    ap.add_argument("--bins", help="comma-separated bin quality values, e.g. 12,24,38")
    ap.add_argument(
        "--nbins",
        type=int,
        help="synthesise N evenly-spaced bins over q=12..39 (sparsity sweep)",
    )
    ap.add_argument(
        "--anchors-only",
        action="store_true",
        help="emit only the bin columns, isolating n_valid",
    )
    args = ap.parse_args()

    rows, nq = load_trans(args.json)

    if args.nbins:
        bins = [int(round(12 + (39 - 12) * i / (args.nbins - 1))) for i in range(args.nbins)]
    elif args.bins:
        bins = [int(b) for b in args.bins.split(",")]
    else:
        bins = None

    if bins:
        bad = [b for b in bins if not 0 <= b < nq]
        if bad:
            sys.exit(f"bin values {bad} outside the observed q range 0..{nq - 1}")
        counts, grid = rebin(rows, nq, bins, args.anchors_only)
    else:
        counts = {k: {q: v[q] for q in range(nq)} for k, v in rows.items()}
        grid = list(range(nq))

    with open(args.out, "w") as fh:
        fh.write("pair\tq\terrs\ttot\trlogp\n")
        for nti in NT:
            tot = {q: sum(counts[nti + "2" + y][q] for y in NT) for q in grid}
            for ntj in NT:
                if nti == ntj:
                    continue
                for q in grid:
                    e, t = counts[nti + "2" + ntj][q], tot[q]
                    r = math.log10((e + 1) / t) if t > 0 else float("nan")
                    fh.write(f"{nti}2{ntj}\t{q}\t{e}\t{t}\t{r}\n")

    populated = sum(1 for q in grid if any(counts[k][q] for k in counts))
    print(
        f"wrote {args.out}: {len(grid)} q columns, {populated} populated"
        + (f", bins={bins}" if bins else "")
    )


if __name__ == "__main__":
    main()
