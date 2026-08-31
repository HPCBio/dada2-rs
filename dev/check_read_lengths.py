#!/usr/bin/env python3
"""check_read_lengths.py — how length-variable are the reads reaching `dada`?

Reports the read-length distribution per file, and specifically the fraction of
reads that differ from the modal length, which is the number that matters and the
one a min/median/max summary hides.

WHAT THIS IS *NOT* FOR. An earlier version of this script claimed that a length
difference makes the k-mer screen and a minimizer sketch behave differently,
because `kmer_dist8` normalises by `min(len) - k + 1` while a sketch normalises
by `min(total_a, total_b)`. **That distinction is false.** A length-L read has
exactly `L - k + 1` k-mers, so

    min(len_a, len_b) - k + 1  ==  min(#kmers_a, #kmers_b)

which is the *same* min-based normalisation the sketch uses. Both screens are
length-adaptive in the same way. The real difference between them is WHAT IS
COUNTED -- every k-mer versus a winnowed sample -- not how it is scaled.

Empirically the two agree with this: PacBio HiFi is 85% off-modal-length and
shows no systematic bias between the screens (net/gross 0.172, indistinguishable
from fixed-length Illumina's 0.150).

So read length is still worth knowing -- it drives alignment cost, band
sufficiency, and whether merging will overlap -- but do not expect it to predict
a screen difference.

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
    offs = []
    for f in files:
        ls = lengths(f)
        if not ls:
            print(f"{os.path.basename(f)[-44:]:>44s} {'0':>9s}")
            continue
        mode = statistics.mode(ls)
        off = sum(1 for x in ls if x != mode)
        distinct = len(set(ls))
        offs.append(100 * off / len(ls))
        print(f"{os.path.basename(f)[-44:]:>44s} {len(ls):9,d} {min(ls):6d} "
              f"{int(statistics.median(ls)):7d} {max(ls):6d} {distinct:9d} {100*off/len(ls):7.2f}%")

    print()
    if not offs:
        return
    med_off = statistics.median(offs)
    print(f"SUMMARY: median {med_off:.2f}% of reads differ from the modal length "
          f"(range {min(offs):.2f}-{max(offs):.2f}% across {len(offs)} files).")
    print()
    print("  For scale, on datasets already measured in this project:")
    print("    MiSeq SOP, truncLen applied :  0.00% off-mode  (fixed length)")
    print("    PacBio HiFi, quality-trimmed : 85.30% off-mode")
    print()
    if med_off < 0.5:
        print("  These reads are effectively FIXED-LENGTH. If a fixed truncLen was")
        print("  applied, note that it also caps any biological length variation:")
        print("  on a length-variable amplicon the variation then exists only")
        print("  post-merge, after the screen has already run.")
    else:
        print("  These reads are variable-length. That affects alignment cost and")
        print("  merge overlap. It does NOT by itself predict a difference between")
        print("  the k-mer and minimizer screens -- both normalise by the shorter")
        print("  sequence, and PacBio at 85% off-mode shows no bias between them.")


if __name__ == "__main__":
    main()
