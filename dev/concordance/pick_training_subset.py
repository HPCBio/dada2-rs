#!/usr/bin/env python3
"""pick_training_subset.py — choose a fixed set of samples to train error models on.

Once a run is large enough that `nbases` truncates, each tool draws its own
training reads and an arm comparison is measuring the draw as well as whatever
was under test (issue #205). Pinning the training set removes that: both tools
learn from exactly these samples, with `nbases` set high enough that neither
subsamples within them, and then denoise everything.

It is also closer to real use. Nobody trains on a whole run; the question is
whether a change survives a realistic training subset, not an exhaustive one.

Samples are taken in sorted order until the base budget is reached, so the
choice is deterministic and the manifest is reproducible from the same inputs.

Usage:
    pick_training_subset.py <filtered-dir> --bases 3e8 [--suffix F.fastq.gz] \
        [--out train_samples.txt]

Prints the chosen sample names (one per line) and a summary to stderr.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def bases_in(path: Path) -> tuple[int, int]:
    """(reads, sequence bases) in a FASTQ, counting sequence lines only."""
    reads = bases = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                reads += 1
                bases += len(line) - 1
    return reads, bases


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", type=Path, help="directory of filtered FASTQs")
    ap.add_argument("--bases", type=float, required=True,
                    help="target training bases, forward reads only (e.g. 3e8)")
    ap.add_argument("--suffix", default="F.fastq.gz",
                    help="forward-read suffix identifying a sample [default: F.fastq.gz]")
    ap.add_argument("--out", type=Path, help="write the manifest here (default: stdout)")
    args = ap.parse_args()

    fwds = sorted(p for p in args.dir.iterdir() if p.name.endswith(args.suffix))
    if not fwds:
        print(f"no *{args.suffix} in {args.dir}", file=sys.stderr)
        return 1

    chosen, total, reads_tot = [], 0, 0
    for p in fwds:
        if total >= args.bases:
            break
        reads, bases = bases_in(p)
        chosen.append(p.name[: -len(args.suffix)])
        total += bases
        reads_tot += reads

    text = "\n".join(chosen) + "\n"
    if args.out:
        args.out.write_text(text)
    else:
        sys.stdout.write(text)

    print(f"{len(chosen)} of {len(fwds)} samples, {reads_tot:,} reads, "
          f"{total:,} forward bases (target {args.bases:,.0f})", file=sys.stderr)
    if total < args.bases:
        print("warning: ran out of samples before reaching the target", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
