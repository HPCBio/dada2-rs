#!/usr/bin/env python3
"""check_read_lengths.py — do the reads reaching `dada` actually vary in length?

The pre-alignment screen runs on per-sample uniques, per orientation, BEFORE
merging. `kmer_dist8` normalises by `min(len) - k + 1`, so its length behaviour
only differs from a minimizer sketch's `min(total_a, total_b)` when the reads
being compared have DIFFERENT lengths.

That makes this a precondition check, not a curiosity. On a length-variable
amplicon such as ITS2 the interesting mechanism is live only if the reads fed to
`dada` are still variable-length. If filtering applied a fixed `truncLen`, every
read is the same length, the mechanism never engages, and a sweep on that data
cannot test it — the length variation would exist only after merging, by which
point the screen has already run.

Running this first avoids the failure this project keeps hitting: an instrument
that cannot exhibit the effect it is pointed at.

Usage:
  check_read_lengths.py <fastq.gz> [<fastq.gz> ...]
  check_read_lengths.py <dir>            # all *.fastq.gz within
"""

import glob
import gzip
import os
import statistics
import sys


def lengths(path, cap=200000):
    out = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(len(line.rstrip("\n")))
                if len(out) >= cap:
                    break
    return out


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    files = []
    for a in args:
        if os.path.isdir(a):
            files += sorted(glob.glob(os.path.join(a, "*.fastq.gz"))
                            + glob.glob(os.path.join(a, "*.fq.gz"))
                            + glob.glob(os.path.join(a, "*.fastq")))
        else:
            files.append(a)
    if not files:
        sys.exit("no FASTQ found")

    print(f"{'file':>44s} {'reads':>9s} {'min':>6s} {'median':>7s} {'max':>6s} {'distinct':>9s} {'%!=mode':>8s}")
    all_uniform = True
    for f in files:
        ls = lengths(f)
        if not ls:
            print(f"{os.path.basename(f)[-44:]:>44s} {'0':>9s}")
            continue
        mode = statistics.mode(ls)
        off = sum(1 for x in ls if x != mode)
        distinct = len(set(ls))
        if distinct > 1 and off / len(ls) > 0.01:
            all_uniform = False
        print(f"{os.path.basename(f)[-44:]:>44s} {len(ls):9,d} {min(ls):6d} "
              f"{int(statistics.median(ls)):7d} {max(ls):6d} {distinct:9d} {100*off/len(ls):7.2f}%")

    print()
    if all_uniform:
        print("VERDICT: reads are effectively FIXED-LENGTH.")
        print("  The k-mer screen's `min(len) - k + 1` normalisation never sees a")
        print("  length mismatch, so a minimizer sketch cannot differ from it for")
        print("  length reasons on this data. A sweep here measures the screens'")
        print("  general behaviour -- it does NOT test the length-variability")
        print("  hypothesis. If that is the question, the filtering step applied a")
        print("  fixed truncLen and needs to be re-run without one (the DADA2 ITS")
        print("  workflow omits truncLen for exactly this reason).")
    else:
        print("VERDICT: reads are VARIABLE-LENGTH.")
        print("  The length-normalisation difference between the two screens is")
        print("  live on this data, so a sweep here can test it. Watch net/gross in")
        print("  compare_counts_correlation.py: a length artefact should show up as")
        print("  SYSTEMATIC bias (net/gross well above 0.3), not as the cancelling")
        print("  reassignment seen on fixed-length 16S.")


if __name__ == "__main__":
    main()
