#!/usr/bin/env python3
"""bench_table.py — distill bench_pooled.py `summary.csv` files into Markdown.

The benchmark runs happen on the cluster (the datasets are large); this turns
their `summary.csv` output into Markdown tables you can paste straight into the
docs / README.

Two modes:

  Scorecard (default) — one row per run, end-to-end head-to-head:
    python3 bench_table.py \
        "MiSeq pooled=bench_true/summary.csv" \
        "MiSeq per-sample=bench_false/summary.csv" \
        "MiSeq pseudo=bench_pseudo/summary.csv"

  Per-step (one run) — the step-by-step breakdown:
    python3 bench_table.py --per-step bench_pseudo/summary.csv

`LABEL=path` lets you name each run; a bare `path` is labeled by its parent dir.
R columns compare against the R-single rows (fair end-to-end wall + overall RSS);
rust-only runs (no --run-r) just omit the R/speedup columns.
"""
import argparse
import csv
import sys
from pathlib import Path


def load(path):
    """Return {(stack, step): {wall_s, cpu_s, cores, maxrss_kb}} from a summary.csv."""
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows[(r["stack"], r["step"])] = r
    return rows


def fnum(rows, stack, step, col):
    r = rows.get((stack, step))
    if not r or r.get(col) in (None, ""):
        return None
    try:
        return float(r[col])
    except ValueError:
        return None


def mb(kb):
    return f"{kb / 1024:.0f}" if kb else "—"


def secs(s):
    return f"{s:.1f}" if s is not None else "—"


def ratio(num, den):
    return f"{num / den:.1f}×" if (num and den) else "—"


def label_for(spec):
    if "=" in spec:
        lab, path = spec.split("=", 1)
        return lab.strip(), Path(path.strip())
    p = Path(spec)
    return p.parent.name or p.stem, p


def scorecard(specs, rstack):
    # RSS ratio is dada2-rs ÷ R (fraction of the comparator's memory): <1× =
    # dada2-rs uses less, >1× = more. This reads as a reduction directly and
    # handles wins+increases in one column, so it runs OPPOSITE to Speedup
    # (R ÷ dada2-rs, where higher is better). See docs/results.md methodology.
    print("| Run | dada2-rs wall (s) | R wall (s) | Speedup | dada2-rs peak (MB) | R peak (MB) | Peak RSS (rs÷R) |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    for spec in specs:
        lab, path = label_for(spec)
        rows = load(path)
        rw = fnum(rows, "dada2-rs", "TOTAL", "wall_s")
        rr = fnum(rows, rstack, "TOTAL", "wall_s")
        rrss = fnum(rows, "dada2-rs", "TOTAL", "maxrss_kb")
        rrss_r = fnum(rows, rstack, "TOTAL", "maxrss_kb")
        print(f"| {lab} | {secs(rw)} | {secs(rr)} | {ratio(rr, rw)} | "
              f"{mb(rrss)} | {mb(rrss_r)} | {ratio(rrss, rrss_r)} |")


# Pipeline step order (illumina superset; pacbio steps are a subset).
STEP_ORDER = ["remove_primers", "filter", "learn_fwd", "learn_rev", "learn",
              "dada_fwd", "dada_rev", "dada", "merge", "make_table",
              "remove_bimera", "TOTAL"]


def per_step(spec, rstack):
    lab, path = label_for(spec)
    rows = load(path)
    steps = [s for (stk, s) in rows if stk == "dada2-rs"]
    steps = [s for s in STEP_ORDER if s in steps] + \
            [s for s in steps if s not in STEP_ORDER]
    print(f"**{lab}** — dada2-rs vs R (R-single end-to-end wall)\n")
    print("| Step | dada2-rs wall (s) | dada2-rs cores | dada2-rs peak (MB) | R wall (s) | Speedup |")
    print("|---|---:|---:|---:|---:|---:|")
    for s in steps:
        rw = fnum(rows, "dada2-rs", s, "wall_s")
        rc = fnum(rows, "dada2-rs", s, "cores")
        rrss = fnum(rows, "dada2-rs", s, "maxrss_kb")
        Rw = fnum(rows, rstack, s, "wall_s")
        cores = f"{rc:.1f}" if rc is not None else "—"
        print(f"| {s} | {secs(rw)} | {cores} | {mb(rrss)} | {secs(Rw)} | {ratio(Rw, rw)} |")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="+", metavar="[LABEL=]summary.csv")
    p.add_argument("--per-step", action="store_true",
                   help="emit a per-step table for a single run (instead of the scorecard)")
    p.add_argument("--r", choices=["single", "split"], default="single",
                   help="which R stack to compare against (default: single = fair end-to-end)")
    args = p.parse_args()
    rstack = f"R-{args.r}"

    if args.per_step:
        if len(args.runs) != 1:
            sys.exit("--per-step takes exactly one summary.csv")
        per_step(args.runs[0], rstack)
    else:
        scorecard(args.runs, rstack)


if __name__ == "__main__":
    main()
