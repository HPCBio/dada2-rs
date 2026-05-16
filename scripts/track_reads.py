#!/usr/bin/env python3
"""Summarize per-sample read counts across dada2-rs pipeline steps.

Reads a combination of JSON outputs and emits a TSV with one row per sample.

Inputs
------
--filter-and-trim   one or more filter-and-trim JSONs (one per sample)
--dada              one (single-end) or two (paired-end: R1, R2) make-sequence-table
                    JSONs built from the per-sample dada / dada-pooled outputs
--merge-pairs       a single merge-pairs JSON (paired-end only)
--seqtab            the final sequence table after remove-bimera-denovo

Output columns are emitted only when the corresponding input was supplied.
For paired-end runs (two --dada files) the denoised columns are split into
`denoisedF` and `denoisedR`; for single-end runs there is one `denoised`
column.  Samples present in some steps but not others get 0 in the missing
columns.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
from pathlib import Path
from typing import Any


def _open_json(path: Path) -> Any:
    """Read a JSON file, transparently handling .gz."""
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as raw:
            return json.load(io.TextIOWrapper(raw, encoding="utf-8"))
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _unwrap(doc: Any, expected_types: set[str]) -> dict:
    """Strip the {type, data} envelope our subcommands wrap outputs in."""
    if isinstance(doc, dict) and "type" in doc and "data" in doc:
        if doc["type"] not in expected_types:
            raise ValueError(
                f"unexpected JSON type '{doc['type']}'; expected one of {sorted(expected_types)}"
            )
        return doc["data"]
    return doc


def load_filter_trim(paths: list[Path], counts: dict[str, dict[str, int]]) -> None:
    for p in paths:
        data = _unwrap(_open_json(p), {"filter-and-trim"})
        sample = data["sample"]
        row = counts.setdefault(sample, {})
        row["input"] = int(data.get("reads_in", 0))
        row["filtered"] = int(data.get("reads_out", 0))


def _load_seqtab_rowsums(path: Path) -> dict[str, int]:
    """Read a make-sequence-table or remove-bimera-denovo JSON and return
    {sample: row-sum} per sample."""
    data = _unwrap(_open_json(path), {"make-sequence-table", "remove-bimera-denovo"})
    out: dict[str, int] = {}
    for sample, row_counts in zip(data.get("samples", []), data.get("counts", [])):
        out[sample] = int(sum(row_counts))
    return out


def load_dada_seqtabs(
    paths: list[Path], counts: dict[str, dict[str, int]]
) -> bool:
    """Populate denoised counts from one or two sequence-table JSONs.

    Returns True when paired-end (two paths supplied), False otherwise.
    """
    if not paths:
        return False
    if len(paths) > 2:
        raise SystemExit("--dada accepts at most 2 sequence-table JSONs (R1, R2)")

    paired = len(paths) == 2
    for idx, p in enumerate(paths):
        rowsums = _load_seqtab_rowsums(p)
        col = ("denoisedF", "denoisedR")[idx] if paired else "denoised"
        for sample, n in rowsums.items():
            counts.setdefault(sample, {})[col] = n
    return paired


def load_merge_pairs(path: Path | None, counts: dict[str, dict[str, int]]) -> None:
    if path is None:
        return
    data = _unwrap(_open_json(path), {"merge-pairs"})
    for entry in data.get("samples", []):
        sample = entry["sample"]
        row = counts.setdefault(sample, {})
        row["merged"] = int(entry.get("accepted_pairs", 0))


def load_seqtab(path: Path | None, counts: dict[str, dict[str, int]]) -> None:
    if path is None:
        return
    rowsums = _load_seqtab_rowsums(path)
    for sample, n in rowsums.items():
        counts.setdefault(sample, {})["nochim"] = n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize per-sample read counts across dada2-rs steps.",
    )
    ap.add_argument(
        "-f", "--filter-and-trim", nargs="+", type=Path, default=[], metavar="JSON",
        help="One or more filter-and-trim output JSONs (one per sample).",
    )
    ap.add_argument(
        "-d", "--dada", nargs="+", type=Path, default=[], metavar="JSON",
        help="One (single-end) or two (paired-end: R1 R2) make-sequence-table JSONs "
             "built from the dada / dada-pooled per-sample outputs.",
    )
    ap.add_argument(
        "-m", "--merge-pairs", type=Path, default=None, metavar="JSON",
        help="merge-pairs output JSON (paired-end only).  Omit for single-end runs.",
    )
    ap.add_argument(
        "-s", "--seqtab", type=Path, default=None, metavar="JSON",
        help="Final sequence-table JSON after remove-bimera-denovo.",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Write TSV to this file instead of stdout.",
    )
    args = ap.parse_args(argv)

    if not (args.filter_and_trim or args.dada or args.merge_pairs or args.seqtab):
        ap.error("supply at least one of --filter-and-trim / --dada / --merge-pairs / --seqtab")

    if args.merge_pairs is not None and len(args.dada) < 2:
        ap.error("--merge-pairs is only meaningful with paired-end --dada (two seq tables)")

    counts: dict[str, dict[str, int]] = {}
    load_filter_trim(args.filter_and_trim, counts)
    paired = load_dada_seqtabs(args.dada, counts)
    load_merge_pairs(args.merge_pairs, counts)
    load_seqtab(args.seqtab, counts)

    cols: list[str] = []
    if args.filter_and_trim:
        cols += ["input", "filtered"]
    if args.dada:
        cols += ["denoisedF", "denoisedR"] if paired else ["denoised"]
    if args.merge_pairs is not None:
        cols += ["merged"]
    if args.seqtab is not None:
        cols += ["nochim"]

    out: io.TextIOBase
    if args.output is None:
        out = sys.stdout
    else:
        out = args.output.open("w", encoding="utf-8")

    try:
        out.write("sample\t" + "\t".join(cols) + "\n")
        for sample in sorted(counts):
            row = counts[sample]
            vals = [str(row.get(c, 0)) for c in cols]
            out.write(sample + "\t" + "\t".join(vals) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
