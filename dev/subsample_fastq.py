#!/usr/bin/env python3
"""Deterministic, pair-safe, NESTED subsampling of FASTQ files.

Built for depth-ladder experiments: rarefy a run to several fractions and ask
how an inference result changes as coverage per variant falls.

PREFER `seqkit` IF YOU HAVE IT. `seqkit sample -p F -s SEED` is faster and was
verified (seqkit 2.x) to be reproducible, pair-synced, and nested across
fractions when the seed is held fixed -- which is everything this script
provides. Use the same seed at every rung:

    seqkit sample -p 0.25 -s 11 in_1.fastq.gz -o out_1.fastq.gz
    seqkit sample -p 0.25 -s 11 in_2.fastq.gz -o out_2.fastq.gz

This script exists for two cases seqkit does not cover by construction:

  * seqkit's nesting is an emergent property of drawing one PRNG value per
    record in file order, not a documented guarantee -- re-verify it after any
    seqkit upgrade (see "verifying" below).
  * seqkit's pair sync is POSITIONAL. It assumes R1 and R2 are in identical
    order. If anything reorders or independently filters one mate, it
    desynchronises silently.

Three properties make this different from position-based samplers:

1. **Pair safety without seeds.** The keep/drop decision is a hash of the read
   NAME, not a random draw, so R1 and R2 make the same decision independently.
   Pairs cannot desynchronise even if the two files are processed separately,
   on different machines, or in different orders.

2. **The ladder is NESTED.** Every read kept at 12.5% is also kept at 25%, 50%
   and 100%. Rungs therefore differ only by additional removal, not by
   independent draws, which removes sampling noise from between-rung
   comparisons -- exactly what you want when the question is "how does the
   result move as depth falls" rather than "what is the variance at depth X".
   Independent draws at each rung would add noise that can easily exceed the
   effect being measured.

3. **One decompression pass for the whole ladder.** All fractions are written
   simultaneously.

Because decisions are name-derived, the same fraction is reproducible forever
without recording a seed, and two different datasets rarefied to the same
fraction are NOT correlated (names differ). Use --salt to draw a different
nested series from the same input.

Usage
-----
Whole run (batch mode):

    ./subsample_fastq.py --indir trimmed/ --outdir sub/ \\
        --fractions 0.5,0.25,0.125

One pair:

    ./subsample_fastq.py --r1 A_1.fastq.gz --r2 A_2.fastq.gz \\
        --outdir sub/ --fractions 0.5,0.25

Single-end: pass --r1 only, or --batch-single in batch mode.

Output layout mirrors the input basenames under one directory per fraction:

    sub/frac_0.500/A_1.fastq.gz
    sub/frac_0.500/A_2.fastq.gz
    sub/frac_0.250/A_1.fastq.gz
    ...

A TSV of realized counts is written to <outdir>/subsample_counts.tsv.

Note on rarefaction: subsampling reads reduces observed richness sublinearly,
so median coverage per surviving feature does NOT fall in proportion to depth.
Plan ladders against realized coverage, not against the nominal fraction.

Verifying a sampler (works for this script or for seqkit)
---------------------------------------------------------
    names() { gzcat "$1" | awk 'NR%4==1{print $1}' | sort; }
    # nesting: must be 0
    comm -23 <(names frac_0.125/S_1.fastq.gz) <(names frac_0.250/S_1.fastq.gz) | wc -l
    # the same check must be ABLE to fail, or it proves nothing:
    comm -23 <(names frac_0.250/S_1.fastq.gz) <(names frac_0.125/S_1.fastq.gz) | wc -l
    # pair sync: must be 0
    diff <(names frac_0.250/S_1.fastq.gz) <(names frac_0.250/S_2.fastq.gz) | wc -l
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import sys

# Reads whose hash falls below frac * 2**64 are kept.
_SCALE = float(1 << 64)


def _open_read(path: str) -> io.TextIOWrapper:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii")
    return open(path, "r", encoding="ascii")


def _open_write(path: str, level: int):
    if path.endswith(".gz"):
        return io.TextIOWrapper(
            gzip.open(path, "wb", compresslevel=level), encoding="ascii"
        )
    return open(path, "w", encoding="ascii")


def read_key(header: str, salt: bytes) -> int:
    """Hash a FASTQ header down to the mate-independent read identity.

    Takes the first whitespace-delimited token and strips a trailing /1 or /2.
    This covers both SRA-style ("@SRR123.1 1 length=250") and Illumina-native
    ("@inst:run:...:x:y 1:N:0:INDEX") naming, in both of which the first token
    is already identical between mates.
    """
    name = header.split(None, 1)[0]
    if name.startswith("@"):
        name = name[1:]
    if len(name) > 2 and name[-2] == "/" and name[-1] in "12":
        name = name[:-2]
    digest = hashlib.blake2b(salt + name.encode("ascii"), digest_size=8).digest()
    return int.from_bytes(digest, "big")


def iter_records(fh):
    """Yield 4-line FASTQ records as a single string, plus the header."""
    while True:
        header = fh.readline()
        if not header:
            return
        rest = fh.readline() + fh.readline() + fh.readline()
        yield header, header + rest


def subsample_pair(r1, r2, outputs, salt, thresholds):
    """Stream one (pair of) file(s), writing to every fraction that keeps it.

    `outputs` is a list of (r1_handle, r2_handle_or_None) parallel to
    `thresholds`. Returns (total, [kept_per_fraction]).
    """
    kept = [0] * len(thresholds)
    total = 0

    with _open_read(r1) as f1:
        f2 = _open_read(r2) if r2 else None
        try:
            g1 = iter_records(f1)
            g2 = iter_records(f2) if f2 else None

            for header, rec1 in g1:
                total += 1
                rec2 = None
                if g2 is not None:
                    try:
                        h2, rec2 = next(g2)
                    except StopIteration:
                        raise SystemExit(
                            f"error: {r2} ended early (R1 has more reads than R2) "
                            f"at record {total}"
                        )
                    # Pair-order check. Cheap, and catches the failure mode that
                    # silently corrupts every downstream result.
                    if read_key(h2, salt) != read_key(header, salt):
                        raise SystemExit(
                            f"error: read name mismatch at record {total}\n"
                            f"  R1: {header.strip()}\n  R2: {h2.strip()}\n"
                            "Files are not in the same order; re-sync before "
                            "subsampling."
                        )

                h = read_key(header, salt)
                for i, thr in enumerate(thresholds):
                    if h < thr:
                        kept[i] += 1
                        outputs[i][0].write(rec1)
                        if rec2 is not None:
                            outputs[i][1].write(rec2)

            if g2 is not None:
                try:
                    next(g2)
                except StopIteration:
                    pass
                else:
                    raise SystemExit(
                        f"error: {r2} has more reads than {r1}"
                    )
        finally:
            if f2 is not None:
                f2.close()

    return total, kept


def main():
    p = argparse.ArgumentParser(
        description="Deterministic nested subsampling of paired FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage\n-----\n", 1)[1],
    )
    p.add_argument("--indir", help="directory of paired FASTQ files (batch mode)")
    p.add_argument("--r1", help="single R1 file")
    p.add_argument("--r2", help="single R2 file (omit for single-end)")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument(
        "--fractions",
        required=True,
        help="comma-separated fractions in (0,1], e.g. 0.5,0.25,0.125",
    )
    p.add_argument("--r1-suffix", default="_1.fastq.gz")
    p.add_argument("--r2-suffix", default="_2.fastq.gz")
    p.add_argument(
        "--batch-single",
        action="store_true",
        help="batch mode: treat every --r1-suffix file as single-end",
    )
    p.add_argument(
        "--salt",
        default="",
        help="change this to draw a DIFFERENT nested series from the same input "
        "(default: empty, i.e. the canonical series)",
    )
    p.add_argument(
        "--compresslevel",
        type=int,
        default=1,
        help="gzip level for outputs (default 1: these are intermediates, and "
        "level 1 is several times faster for a few %% more size)",
    )
    args = p.parse_args()

    try:
        fractions = [float(x) for x in args.fractions.split(",") if x.strip()]
    except ValueError:
        p.error("--fractions must be comma-separated numbers")
    if not fractions:
        p.error("--fractions is empty")
    for f in fractions:
        if not (0.0 < f <= 1.0):
            p.error(f"fraction {f} outside (0, 1]")
    fractions.sort(reverse=True)

    salt = args.salt.encode("utf-8")
    thresholds = [int(f * _SCALE) for f in fractions]

    # Resolve the work list as (r1, r2_or_None) pairs.
    jobs = []
    if args.indir:
        if args.r1 or args.r2:
            p.error("--indir and --r1/--r2 are mutually exclusive")
        names = sorted(os.listdir(args.indir))
        for name in names:
            if not name.endswith(args.r1_suffix):
                continue
            r1 = os.path.join(args.indir, name)
            if args.batch_single:
                jobs.append((r1, None))
                continue
            mate = name[: -len(args.r1_suffix)] + args.r2_suffix
            r2 = os.path.join(args.indir, mate)
            if not os.path.exists(r2):
                raise SystemExit(f"error: no mate for {r1} (expected {r2})")
            jobs.append((r1, r2))
        if not jobs:
            raise SystemExit(
                f"error: no files ending in {args.r1_suffix!r} under {args.indir}"
            )
    elif args.r1:
        jobs.append((args.r1, args.r2))
    else:
        p.error("provide either --indir or --r1")

    frac_dirs = []
    for f in fractions:
        d = os.path.join(args.outdir, f"frac_{f:.3f}")
        os.makedirs(d, exist_ok=True)
        frac_dirs.append(d)

    rows = []
    for r1, r2 in jobs:
        outputs = []
        try:
            for d in frac_dirs:
                o1 = _open_write(
                    os.path.join(d, os.path.basename(r1)), args.compresslevel
                )
                o2 = (
                    _open_write(
                        os.path.join(d, os.path.basename(r2)), args.compresslevel
                    )
                    if r2
                    else None
                )
                outputs.append((o1, o2))

            total, kept = subsample_pair(r1, r2, outputs, salt, thresholds)
        finally:
            for o1, o2 in outputs:
                o1.close()
                if o2 is not None:
                    o2.close()

        sample = os.path.basename(r1)
        if sample.endswith(args.r1_suffix):
            sample = sample[: -len(args.r1_suffix)]
        for f, k in zip(fractions, kept):
            realized = (k / total) if total else 0.0
            rows.append((sample, f, total, k, realized))
        print(
            f"{sample}: {total} reads -> "
            + ", ".join(
                f"{f:.3f}:{k} ({100 * k / total:.2f}%)" if total else f"{f:.3f}:0"
                for f, k in zip(fractions, kept)
            ),
            file=sys.stderr,
        )

    tsv = os.path.join(args.outdir, "subsample_counts.tsv")
    with open(tsv, "w", encoding="ascii") as fh:
        fh.write("sample\tfraction\tinput_reads\tkept_reads\trealized_fraction\n")
        for sample, f, total, k, realized in rows:
            fh.write(f"{sample}\t{f:.3f}\t{total}\t{k}\t{realized:.6f}\n")

    grand = {}
    for _, f, total, k, _ in rows:
        t, kk = grand.get(f, (0, 0))
        grand[f] = (t + total, kk + k)
    print("", file=sys.stderr)
    for f in fractions:
        t, k = grand[f]
        print(
            f"frac {f:.3f}: {k}/{t} reads kept ({100 * k / t:.3f}%)"
            if t
            else f"frac {f:.3f}: no reads",
            file=sys.stderr,
        )
    print(f"\ncounts written to {tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
