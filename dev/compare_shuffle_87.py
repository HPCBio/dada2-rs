#!/usr/bin/env python3
"""Collect the shuffle-phase #87/#139 diagnostics across runs into one table.

REVIVED for #139 (was retired in PR #125 when #87 closed). The projection line
it parses has been re-derived: #87 priced the work that carrying `compmax`
across buds would relocate at the *full-rescan* reconcile's scattered rate, and
#136 replaced that reconcile with an incremental update. The relocated cost now
splits into a sequential cluster-major walk and a much smaller scattered rescan,
so the `f` values here are NOT comparable with the 0.34-0.71 quoted in #87's
findings page -- those were full-rescan comps over build comps.

Pulls the `[dada] phase times`, `shuffle scan split`, `shuffle scan time` and
`#87 projection` lines out of each run's dada logs and lays them side by side,
so the cross-run comparison of the relocation fraction `f` is readable at a
glance. That cross-run comparison is the thing that decides #87: a single
run's WIN/LOSS verdict only covers that run's cluster count, whereas `f`
across platforms says whether a cluster-count-gated version could pay off.

The diagnostics come from PR #123 (branch perf/shuffle-scan-time-87); logs
from earlier builds are reported with blanks rather than skipped, so a stale
run is visible rather than silently absent.

Usage:
    dev/compare_shuffle_87.py tmp/full-pooling tmp/full-pooling-ITS ...
    dev/compare_shuffle_87.py --csv out.csv tmp/full-pooling*
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# The projection line is the payload; the rest is context for interpreting it.
RE_PHASE = re.compile(
    r"compare=([\d.]+)s .*?shuffle=([\d.]+)s\s+bud=([\d.]+)s\s+p_update=([\d.]+)s"
)
RE_SPLIT = re.compile(
    r"build=(\d+) \((\d+)% of scanned\) over (\d+) builds \(([\d.]+) comps/build\), "
    r"reconcile=(\d+)"
)
RE_TIME = re.compile(
    r"build=([\d.]+)s \(([\d.]+) ns/comp\)\s+reconcile=([\d.]+)s \(([\d.]+) ns/comp\)\s+"
    r"build=(\d+)% of shuffle time"
)
RE_PROJ = re.compile(
    r"post-bud reconcile=(\d+) over (\d+) buds, f=([\d.]+) vs break-even ([\d.]+) "
    r"→ (WIN|LOSS) \(projected ([+-][\d.]+)s on shuffle=([\d.]+)s\)"
)


def parse_log(path: Path) -> dict | None:
    """Extract one log's diagnostics. Returns None if it has no shuffle lines."""
    text = path.read_text(errors="replace")
    row: dict = {"log": str(path)}

    if m := RE_PHASE.search(text):
        row["compare_s"] = float(m.group(1))
        row["shuffle_s"] = float(m.group(2))
        row["p_update_s"] = float(m.group(4))
    if m := RE_SPLIT.search(text):
        row["comps_build"] = int(m.group(1))
        row["build_pct_comps"] = int(m.group(2))
        row["builds"] = int(m.group(3))
        row["comps_reconcile"] = int(m.group(5))
    if m := RE_TIME.search(text):
        row["build_s"] = float(m.group(1))
        row["build_ns"] = float(m.group(2))
        row["reconcile_s"] = float(m.group(3))
        row["reconcile_ns"] = float(m.group(4))
        row["build_pct_time"] = int(m.group(5))
    if m := RE_PROJ.search(text):
        row["comps_first_rec"] = int(m.group(1))
        row["f"] = float(m.group(3))
        row["breakeven"] = float(m.group(4))
        row["verdict"] = m.group(5)
        row["projected_s"] = float(m.group(6))

    # A log with none of the phase instrumentation is not a dada log at all.
    return row if len(row) > 1 else None


def collect(run_dirs: list[Path]) -> list[dict]:
    rows = []
    for d in run_dirs:
        logs = sorted(d.glob("dada/dada.R*.log")) or sorted(d.glob("**/dada.R*.log"))
        if not logs:
            print(f"warning: no dada logs under {d}", file=sys.stderr)
            continue
        for log in logs:
            if row := parse_log(log):
                row["run"] = d.name
                # R1/R2 from the filename; falls back to the stem.
                m = re.search(r"dada\.(R\d+)\.log$", log.name)
                row["read"] = m.group(1) if m else log.stem
                rows.append(row)
            else:
                print(f"warning: no shuffle diagnostics in {log}", file=sys.stderr)
    return rows


def fmt(row: dict, key: str, spec: str = "") -> str:
    v = row.get(key)
    if v is None:
        return "-"
    return format(v, spec) if spec else str(v)


def print_table(rows: list[dict]) -> None:
    hdr = [
        ("run", "run", 26),
        ("read", "read", 5),
        ("shuffle_s", "shuffle", 9),
        ("build_ns", "bld ns", 7),
        ("reconcile_ns", "rec ns", 7),
        ("ratio", "ratio", 6),
        ("build_pct_comps", "bld%cmp", 8),
        ("build_pct_time", "bld%time", 9),
        ("f", "f", 6),
        ("breakeven", "b/even", 7),
        ("verdict", "verdict", 8),
        ("projected_s", "proj s", 8),
    ]
    for r in rows:
        # ns/comp ratio: how much more a relocated comparison costs. This is
        # what makes the comp-count split overstate the build's share.
        if r.get("build_ns") and r.get("reconcile_ns"):
            r["ratio"] = r["reconcile_ns"] / r["build_ns"]

    print("  ".join(label.ljust(w) for _, label, w in hdr))
    print("  ".join("-" * w for _, _, w in hdr))
    specs = {
        "shuffle_s": ".1f",
        "build_ns": ".2f",
        "reconcile_ns": ".2f",
        "ratio": ".2f",
        "f": ".2f",
        "breakeven": ".2f",
        "projected_s": "+.1f",
    }
    for r in rows:
        print("  ".join(fmt(r, k, specs.get(k, "")).ljust(w) for k, _, w in hdr))

    verdicts = {r.get("verdict") for r in rows if r.get("verdict")}
    if not verdicts:
        print("\nNo #87 projection lines found — were these runs built from PR #123?")
        return

    fs = [r["f"] for r in rows if "f" in r]
    print(
        f"\nf range {min(fs):.2f}-{max(fs):.2f} across {len(fs)} logs; "
        f"verdicts: {', '.join(sorted(verdicts))}"
    )
    if verdicts == {"LOSS"}:
        # Worth stating plainly: a uniform LOSS closes the gating question too,
        # since no cluster count in the sweep gets f under break-even.
        print(
            "All LOSS: carrying compmax across buds relocates more work into the\n"
            "scattered reconcile than it removes from the sequential build."
        )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dirs", nargs="+", type=Path)
    ap.add_argument("--csv", type=Path, help="also write the table as CSV")
    args = ap.parse_args()

    rows = collect(args.run_dirs)
    if not rows:
        print("no dada logs found", file=sys.stderr)
        return 1

    print_table(rows)

    if args.csv:
        keys = sorted({k for r in rows for k in r})
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
