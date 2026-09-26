#!/usr/bin/env python3
"""rtab_to_seqtab.py — convert an R DADA2 long-format table to a dada2-rs seqtab JSON.

Written to separate *input* differences from *implementation* differences at the
chimera step. Feeding R's pre-chimera table into `remove-bimera-denovo` (and our
table into R's `removeBimeraDenovo`) gives a 2x2 that a same-pipeline comparison
cannot: if both implementations agree on the same input, the gap is upstream; if
they disagree, it is the bimera search itself.

Input is write_reference.R's `write_long` schema -- sequence,sample,count, one
row per non-zero cell. Output is tagged `make-sequence-table` so
`remove-bimera-denovo` accepts it.

Usage:
    rtab_to_seqtab.py <r_long.csv> <out.json> [--hash md5|sha1]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_in")
    ap.add_argument("json_out")
    ap.add_argument("--hash", default="md5", choices=["md5", "sha1"])
    args = ap.parse_args()

    csv.field_size_limit(1 << 30)  # full-length 16S rows are long
    cells: dict[tuple[str, str], int] = {}
    samples: list[str] = []
    seqs: list[str] = []
    seen_s: set[str] = set()
    seen_q: set[str] = set()
    with open(args.csv_in, newline="") as fh:
        for row in csv.DictReader(fh):
            s, q, c = row["sample"], row["sequence"], int(row["count"])
            if s not in seen_s:
                seen_s.add(s)
                samples.append(s)
            if q not in seen_q:
                seen_q.add(q)
                seqs.append(q)
            cells[(s, q)] = cells.get((s, q), 0) + c

    # Column order follows R's makeSequenceTable: descending total abundance,
    # ties by sequence. Getting this wrong would not change the ASV set but
    # would make any row/column diff against R unreadable.
    tot = {q: 0 for q in seqs}
    for (_, q), c in cells.items():
        tot[q] += c
    seqs.sort(key=lambda q: (-tot[q], q))
    samples.sort()

    h = hashlib.md5 if args.hash == "md5" else hashlib.sha1
    out = {
        "dada2_rs_command": "make-sequence-table",
        "dada2_rs_version": "converted-from-R",
        "samples": samples,
        "sequences": seqs,
        "sequence_ids": [h(q.encode()).hexdigest() for q in seqs],
        "counts": [[cells.get((s, q), 0) for q in seqs] for s in samples],
    }
    with open(args.json_out, "w") as fh:
        json.dump(out, fh)
    print(
        f"wrote {args.json_out}: {len(seqs)} ASVs x {len(samples)} samples, "
        f"{sum(tot.values()):,} reads",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
