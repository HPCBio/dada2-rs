#!/usr/bin/env python3
"""Requantize FASTQ per-base quality into a fixed set of bins.

Simulates Illumina binned-quality output (NovaSeq / NextSeq / i100) from
full-resolution reads, so the error-model's response to quality binning can be
tested in isolation — same reads, same organism, only the quality resolution
changes. Each base's Phred score is mapped to the NEAREST bin value.

Read-level (not derep-level) on purpose: binning is non-linear, so binning each
read then dereplicating reproduces the real binned distribution, whereas binning
a per-unique mean quality would not.

Usage:
  # exact scheme (RECOMMENDED — matches production range->representative maps)
  rebin_quality.py --preset pacbio -o OUTDIR in1.fastq.gz [in2 ...]
  rebin_quality.py --thresholds 40:40,30:35,25:27,20:22,14:17,7:10,0:3 -o OUTDIR ...

  # nearest-value fallback (APPROXIMATE — only for schemes with no official ranges)
  rebin_quality.py --bins 2,12,23,37 -o OUTDIR /path/to/fastq_dir

Output: one `<stem>.binN.fastq.gz` per input (N = number of levels), gzip by the
output extension. Phred offset 33.

Schemes:
  --preset pacbio : PacBio's standard 7-bin QV map (range-based; the exact
                    production cascade [0,6]->3 [7,13]->10 [14,19]->17
                    [20,24]->22 [25,29]->27 [30,39]->35 [40,49]->40).
  Illumina binned deliveries (NovaSeq/NextSeq/i100) are also range-based; supply
  the exact cascade with --thresholds if you have it. --bins is nearest-value and
  will NOT reproduce a range-based scheme faithfully (e.g. it disagrees with the
  PacBio map at Q30/31/38/39), so prefer --preset/--thresholds.
"""

import argparse
import glob
import gzip
import os
import sys

PHRED_OFFSET = 33

# name -> descending (min_qscore, representative) cascade, matching an
# `if q >= min: return rep` chain. PacBio's production QV binning.
PRESETS = {
    "pacbio": [(40, 40), (30, 35), (25, 27), (20, 22), (14, 17), (7, 10), (0, 3)],
}


def lut_from_cascade(cascade):
    """cascade: descending [(min_q, rep), ...]; last entry must have min_q 0."""
    reps = sorted({r for _, r in cascade})
    def rep_for(q):
        for lo, r in cascade:  # descending: first threshold q clears wins
            if q >= lo:
                return r
        return cascade[-1][1]
    trans = {chr(q + PHRED_OFFSET): chr(rep_for(q) + PHRED_OFFSET) for q in range(94)}
    return str.maketrans(trans), reps


def lut_from_bins(bins):
    """Nearest-value fallback: map each Phred to the closest bin representative."""
    trans = {chr(q + PHRED_OFFSET): chr(min(bins, key=lambda b: (abs(b - q), b)) + PHRED_OFFSET)
             for q in range(94)}
    return str.maketrans(trans), sorted(bins)


def _open_r(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _open_w(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def rebin_file(src, dst, table):
    """Stream a FASTQ, translating only the quality line (every 4th)."""
    n = 0
    with _open_r(src) as fi, _open_w(dst) as fo:
        for i, line in enumerate(fi):
            if i % 4 == 3:  # quality line
                fo.write(line.rstrip("\n").translate(table) + "\n")
                n += 1
            else:
                fo.write(line)
    return n


def resolve_inputs(items):
    out = []
    for it in items:
        if os.path.isdir(it):
            for pat in ("*.fastq.gz", "*.fastq", "*.fq.gz", "*.fq"):
                out += sorted(glob.glob(os.path.join(it, pat)))
        else:
            out += sorted(glob.glob(it)) or [it]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", help="FASTQ file(s), glob(s), or a dir")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--preset", choices=sorted(PRESETS), help="named exact scheme")
    g.add_argument("--thresholds",
                   help="exact cascade 'min:rep,...' descending, e.g. 40:40,30:35,...,0:3")
    g.add_argument("--bins", help="APPROXIMATE nearest-value reps, e.g. 2,12,23,37")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--no-gzip", action="store_true", help="write plain .fastq")
    args = ap.parse_args()

    if args.preset:
        table, reps = lut_from_cascade(PRESETS[args.preset])
    elif args.thresholds:
        cascade = sorted(((int(a), int(b)) for a, b in
                          (p.split(":") for p in args.thresholds.split(","))),
                         reverse=True)
        if cascade[-1][0] != 0:
            sys.exit("--thresholds cascade must end at min 0 (e.g. ...,0:3)")
        table, reps = lut_from_cascade(cascade)
    else:
        bins = sorted(int(x) for x in args.bins.split(","))
        if len(bins) < 2:
            sys.exit("need at least 2 bins")
        table, reps = lut_from_bins(bins)
        print("NOTE: --bins is nearest-value (approximate); use --preset/--thresholds "
              "for a faithful range-based scheme.", file=sys.stderr)

    bins = reps
    os.makedirs(args.outdir, exist_ok=True)
    ext = "fastq" if args.no_gzip else "fastq.gz"

    files = resolve_inputs(args.inputs)
    if not files:
        sys.exit("no input FASTQ found")
    print(f"bins={bins}  ({len(bins)} levels)  ->  {args.outdir}")
    total = 0
    for src in files:
        stem = os.path.basename(src)
        for suf in (".gz", ".fastq", ".fq"):
            if stem.endswith(suf):
                stem = stem[: -len(suf)]
        dst = os.path.join(args.outdir, f"{stem}.bin{len(bins)}.{ext}")
        n = rebin_file(src, dst, table)
        total += n
        print(f"  {src} -> {dst}  ({n} reads)")
    print(f"done: {len(files)} file(s), {total} reads")


if __name__ == "__main__":
    main()
