#!/usr/bin/env python3
"""noqual.py — quality-score-free error model, pure stdlib.

Mirrors R DADA2's `noqualErrfun(trans, pseudocount=1)`: aggregates each
transition row's counts across all quality columns, divides by the row sum
(plus pseudocount), and broadcasts the resulting per-transition rate across
every column. Diagonals are left to the dada2-rs init pass to enforce.

Wire format follows the dada2-rs --errfun external contract:
  argv[1] = path to input  trans TSV
  argv[2] = path to output err   TSV

Both files use R's read.table(..., row.names=1, header=TRUE,
check.names=FALSE) layout.

Usage:
  dada2-rs learn-errors *.fastq.gz \\
      --errfun external \\
      --errfun-cmd "python3 examples/external_errfun/noqual.py" \\
      -o err_noqual.json
"""

import sys

PSEUDOCOUNT = 1.0
ROW_LABELS = [
    "A2A", "A2C", "A2G", "A2T",
    "C2A", "C2C", "C2G", "C2T",
    "G2A", "G2C", "G2G", "G2T",
    "T2A", "T2C", "T2G", "T2T",
]


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        sys.stderr.write(f"usage: {argv[0]} <trans.tsv> <err.tsv>\n")
        return 2

    trans_path, err_path = argv[1], argv[2]

    with open(trans_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # header[0] is the row-name column header (often empty); rest are q-scores
        qual_cols = header[1:]
        nq = len(qual_cols)

        rows = {}
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            label, vals = fields[0], [int(x) for x in fields[1:]]
            if len(vals) != nq:
                sys.stderr.write(
                    f"row {label!r} has {len(vals)} cols, expected {nq}\n"
                )
                return 1
            rows[label] = vals

    missing = [lbl for lbl in ROW_LABELS if lbl not in rows]
    if missing:
        sys.stderr.write(f"missing rows in input: {missing}\n")
        return 1

    # Row sums per ref nt (sum of A2A + A2C + A2G + A2T over all q, etc.)
    err_rows = []
    for ref_nt in "ACGT":
        ref_row_labels = [f"{ref_nt}2{q}" for q in "ACGT"]
        ref_total = sum(sum(rows[lbl]) for lbl in ref_row_labels) + PSEUDOCOUNT
        for lbl in ref_row_labels:
            row_total = sum(rows[lbl])
            rate = row_total / ref_total
            err_rows.append((lbl, [rate] * nq))

    with open(err_path, "w") as out:
        out.write("\t" + "\t".join(qual_cols) + "\n")
        for label, vals in err_rows:
            out.write(label)
            for v in vals:
                out.write(f"\t{v:.10g}")
            out.write("\n")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
