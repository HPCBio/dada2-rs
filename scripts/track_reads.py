#!/usr/bin/env python3
"""Summarize per-sample read counts across dada2-rs pipeline steps.

Reads any combination of JSON outputs from filter-and-trim, dada / dada-pooled,
merge-pairs, and the final sequence table (make-sequence-table or
remove-bimera-denovo) and emits a TSV with one row per sample.

Inputs are gathered by --filter-and-trim / --dada / --merge-pairs / --seqtab.
Each flag accepts one or more JSON paths, or shell globs the user has
already expanded.  Steps the user did not run (e.g. merge-pairs for
single-end data, or a final seq table) can simply be omitted; samples that
appear in some steps but not others get 0 in the missing columns.

Output columns (only those for which inputs were supplied):
    sample, input, filtered, denoised, merged, nochim
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
        # filter-and-trim is single-sample: {sample, reads_in, reads_out}
        sample = data["sample"]
        row = counts.setdefault(sample, {})
        row["input"] = int(data.get("reads_in", 0))
        row["filtered"] = int(data.get("reads_out", 0))


def load_dada(paths: list[Path], counts: dict[str, dict[str, int]]) -> None:
    for p in paths:
        data = _unwrap(_open_json(p), {"dada"})
        # dada / dada-pooled emit one file per sample with {sample, total_reads, asvs:[...]}
        sample = data.get("sample") or p.stem
        total = data.get("total_reads")
        if total is None:
            total = sum(int(a.get("abundance", 0)) for a in data.get("asvs", []))
        row = counts.setdefault(sample, {})
        row["denoised"] = int(total)


def load_merge_pairs(paths: list[Path], counts: dict[str, dict[str, int]]) -> None:
    for p in paths:
        data = _unwrap(_open_json(p), {"merge-pairs"})
        # merge-pairs aggregates many samples in one file: {samples: [{sample, accepted_pairs, ...}]}
        for entry in data.get("samples", []):
            sample = entry["sample"]
            row = counts.setdefault(sample, {})
            row["merged"] = int(entry.get("accepted_pairs", 0))


def load_seqtab(paths: list[Path], counts: dict[str, dict[str, int]]) -> None:
    for p in paths:
        data = _unwrap(_open_json(p), {"make-sequence-table", "remove-bimera-denovo"})
        # SequenceTable: {samples:[...], counts:[[...] per sample]}
        for sample, row_counts in zip(data.get("samples", []), data.get("counts", [])):
            row = counts.setdefault(sample, {})
            row["nochim"] = int(sum(row_counts))


COLUMNS = [
    ("input", "filter_and_trim"),
    ("filtered", "filter_and_trim"),
    ("denoised", "dada"),
    ("merged", "merge_pairs"),
    ("nochim", "seqtab"),
]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Summarize per-sample read counts across dada2-rs steps.",
    )
    ap.add_argument(
        "--filter-and-trim", nargs="+", type=Path, default=[], metavar="JSON",
        help="One or more filter-and-trim output JSONs (one per sample).",
    )
    ap.add_argument(
        "--dada", nargs="+", type=Path, default=[], metavar="JSON",
        help="One or more dada or dada-pooled output JSONs (one per sample).",
    )
    ap.add_argument(
        "--merge-pairs", nargs="+", type=Path, default=[], metavar="JSON",
        help="One or more merge-pairs output JSONs.  Omit for single-end runs.",
    )
    ap.add_argument(
        "--seqtab", nargs="+", type=Path, default=[], metavar="JSON",
        help="Sequence-table JSON (make-sequence-table or remove-bimera-denovo).",
    )
    ap.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Write TSV to this file instead of stdout.",
    )
    args = ap.parse_args(argv)

    if not (args.filter_and_trim or args.dada or args.merge_pairs or args.seqtab):
        ap.error("supply at least one of --filter-and-trim / --dada / --merge-pairs / --seqtab")

    counts: dict[str, dict[str, int]] = {}
    load_filter_trim(args.filter_and_trim, counts)
    load_dada(args.dada, counts)
    load_merge_pairs(args.merge_pairs, counts)
    load_seqtab(args.seqtab, counts)

    # Emit only columns for which the user provided inputs.
    provided = {
        "filter_and_trim": bool(args.filter_and_trim),
        "dada": bool(args.dada),
        "merge_pairs": bool(args.merge_pairs),
        "seqtab": bool(args.seqtab),
    }
    cols = [c for c, src in COLUMNS if provided[src]]

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
