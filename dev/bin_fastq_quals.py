#!/usr/bin/env python3
"""Re-bin FASTQ quality scores to the PacBio Revio 7-level scheme.

Emulates instrument-side quality binning so a full-resolution run (e.g. SequelIIe,
Q0-Q93) can be compared against a natively binned run (Revio) on an identical Q axis.

Apply this to RAW FASTQs, before derep. Dereplication averages quality across the
reads collapsed into a unique, so bin-then-derep and derep-then-bin do not agree;
natively binned data is bin-then-derep.

Note the Q50 branch is deliberately disabled (see --allow-q50): with it off, every
Q >= 40 collapses to Q40, reproducing the exact 7-level grid observed in Revio data
so the two error models share a Q axis and become cross-applicable. On SequelIIe
input that discards real resolution -- ~99.5% of its quality mass sits at Q >= 50
across ~54 distinct values -- which is intended, as a deliberately harsh stress test.

Usage:
    bin_fastq_quals.py in.fastq.gz out.fastq.gz
    bin_fastq_quals.py --allow-q50 in.fastq.gz out.fastq.gz
    ls *.fastq.gz | xargs -P8 -I{} bin_fastq_quals.py {} binned/{}
"""

import argparse
import gzip
import sys

PHRED_OFFSET = 33


def _open(path, mode):
    if path == "-":
        return sys.stdin.buffer if "r" in mode else sys.stdout.buffer
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def build_table(allow_q50: bool, max_q: int = 93) -> bytes:
    """Map every input Phred char to its binned Phred char, as a 256-byte table."""

    def bin_quality(q: int) -> int:
        if allow_q50 and q >= 50:
            return 50
        if q >= 40:
            return 40
        if q >= 30:
            return 35
        if q >= 25:
            return 27
        if q >= 20:
            return 22
        if q >= 14:
            return 17
        if q >= 7:
            return 10
        return 3

    table = bytearray(range(256))
    for q in range(max_q + 1):
        table[PHRED_OFFSET + q] = PHRED_OFFSET + bin_quality(q)
    return bytes(table)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", help="input FASTQ (.gz ok, '-' for stdin)")
    ap.add_argument("outfile", help="output FASTQ (.gz ok, '-' for stdout)")
    ap.add_argument("--allow-q50", action="store_true",
                    help="keep Q>=50 as a distinct Q50 level (8-level grid; "
                         "breaks Q-axis compatibility with Revio)")
    args = ap.parse_args()

    table = build_table(args.allow_q50)
    n_reads = 0
    changed_bases = 0
    total_bases = 0

    with _open(args.infile, "rb") as fin, _open(args.outfile, "wb") as fout:
        for i, line in enumerate(fin):
            if i % 4 == 3:
                qual = line.rstrip(b"\n")
                binned = qual.translate(table)
                changed_bases += sum(1 for a, b in zip(qual, binned) if a != b)
                total_bases += len(qual)
                n_reads += 1
                fout.write(binned + b"\n")
            else:
                fout.write(line)

    if total_bases:
        sys.stderr.write(
            "[bin_fastq_quals] %s: %d reads, %d bases, %.2f%% of bases re-binned\n"
            % (args.infile, n_reads, total_bases, 100.0 * changed_bases / total_bases)
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
