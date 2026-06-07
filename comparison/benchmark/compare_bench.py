#!/usr/bin/env python3
"""Compare bench_pooled.py metrics across runs (peak RSS, wall, cores, ...).

bench_pooled.py writes a CSV for every run; this is the standalone diff tool for
those CSVs, the metrics counterpart to compare_asvs.py (which diffs ASV sets).
Run each configuration into its own --outdir, then compare the resulting CSVs
baseline-vs-many. bench_pooled.py needs no changes.

It reads the three CSV shapes bench_pooled.py emits, auto-detected by header
(override with --kind):
  * summary     : summary.csv      — per-step wall_s / cpu_s / cores / maxrss_kb,
                  one block per stack (dada2-rs, R-split, R-single). Keyed by
                  `step`; pick the stack with --stack (default dada2-rs).
  * sweep_jobs  : sweep_jobs.csv   — --sample-jobs-sweep. Keyed by sample_jobs.
  * sweep       : sweep.csv        — --thread-sweep. Keyed by threads.

For each metric it prints one table: rows = steps (or sweep points), columns =
runs, with a Δ% column per non-baseline run (vs the baseline). Lower-is-better
metrics (wall_s, maxrss_kb) and higher-is-better (cores, speedup, efficiency)
are both shown signed; read the sign with the metric in mind.

Usage:
  compare_bench.py --baseline LABEL=PATH --compare LABEL=PATH [--compare LABEL=PATH ...]

  PATH is the CSV itself OR a directory containing it (summary.csv by default,
  or the sweep CSV when --kind is given). LABEL is any short name (e.g. main,
  branch, k5, j4).

Examples:
  # peak RSS + wall, main build vs the deferred-division branch (issue #23)
  compare_bench.py --baseline main=bench_main --compare branch=bench_branch

  # just peak RSS, three runs
  compare_bench.py --metric maxrss_kb \\
      --baseline main=bench_main/summary.csv \\
      --compare b1=bench_b1/summary.csv --compare b2=bench_b2/summary.csv

  # compare a sample-jobs sweep between two builds, machine-readable out
  compare_bench.py --kind sweep_jobs --json diff.json \\
      --baseline main=bench_main --compare branch=bench_branch

Pure stdlib; no dependencies.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# Per CSV kind: the default filename, the key column (row identity), and the
# metric columns to diff. Headers are matched to auto-detect the kind.
KINDS = {
    "summary": {
        "file": "summary.csv",
        "key": "step",
        "metrics": ["wall_s", "cpu_s", "cores", "maxrss_kb"],
        "header": ["stack", "step", "wall_s", "cpu_s", "cores", "maxrss_kb"],
    },
    "sweep_jobs": {
        "file": "sweep_jobs.csv",
        "key": "sample_jobs",
        "metrics": ["wall_s", "cpu_s", "cores", "maxrss_kb", "speedup"],
        "header": ["sample_jobs", "threads_per_job", "wall_s", "cpu_s",
                   "cores", "maxrss_kb", "speedup"],
    },
    "sweep": {
        "file": "sweep.csv",
        "key": "threads",
        "metrics": ["wall_s", "cpu_s", "cores", "speedup", "efficiency"],
        "header": ["threads", "wall_s", "cpu_s", "cores", "speedup", "efficiency"],
    },
}

# Metrics where smaller is better — used only to annotate the Δ direction.
LOWER_IS_BETTER = {"wall_s", "cpu_s", "maxrss_kb"}


def parse_run(spec, default_file):
    """'label=path' -> (label, csv_path). path may be the CSV or its dir."""
    if "=" not in spec:
        sys.exit(f"--baseline/--compare expects LABEL=PATH, got: {spec!r}")
    label, raw = spec.split("=", 1)
    p = Path(raw)
    if p.is_dir():
        p = p / default_file
    if not p.is_file():
        sys.exit(f"{label}: not a file: {p}")
    return label, p


def detect_kind(path):
    with open(path, newline="") as f:
        header = next(csv.reader(f), [])
    for kind, spec in KINDS.items():
        if header == spec["header"]:
            return kind
    sys.exit(f"{path}: unrecognized CSV header {header}; pass --kind explicitly")


def load(path, kind, stack):
    """Return {key: {metric: float}} for the chosen kind (and stack, for summary).
    Keys are kept as strings; order of first appearance is preserved by the caller
    via the returned `order` list."""
    spec = KINDS[kind]
    rows, order = {}, []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if kind == "summary" and r.get("stack") != stack:
                continue
            key = r[spec["key"]]
            vals = {}
            for m in spec["metrics"]:
                v = r.get(m, "")
                vals[m] = float(v) if v not in ("", None) else None
            if key not in rows:
                order.append(key)
            rows[key] = vals
    return rows, order


def fmt(metric, v):
    if v is None:
        return "—"
    if metric == "maxrss_kb":
        mb = v / 1024
        return f"{mb/1024:.2f} GB" if mb >= 1024 else f"{mb:.0f} MB"
    if metric in ("wall_s", "cpu_s"):
        return f"{v:.1f}s"
    if metric in ("speedup",):
        return f"{v:.2f}×"
    if metric in ("efficiency",):
        return f"{v*100:.0f}%"
    return f"{v:.2f}"  # cores and anything else


def pct(base, v):
    if base in (None, 0) or v is None:
        return None
    return (v - base) / base * 100.0


def print_table(metric, keys, runs, base_label):
    """runs = list of (label, rows_dict); base_label is the baseline's label."""
    arrow = "↓ better" if metric in LOWER_IS_BETTER else "↑ better"
    print(f"\n### {metric}  ({arrow})")
    # Column layout: key | base value | (value, Δ%) per compare run.
    klen = max([len("step")] + [len(k) for k in keys])
    head = f"{'':<{klen}}  {base_label:>12}"
    for label, _ in runs[1:]:
        head += f"  {label:>12}{'Δ%':>8}"
    print(head)
    base_rows = runs[0][1]
    for k in keys:
        base_v = base_rows.get(k, {}).get(metric)
        line = f"{k:<{klen}}  {fmt(metric, base_v):>12}"
        for _, rows in runs[1:]:
            v = rows.get(k, {}).get(metric)
            d = pct(base_v, v)
            dstr = "—" if d is None else f"{d:+.1f}%"
            line += f"  {fmt(metric, v):>12}{dstr:>8}"
        print(line)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--baseline", required=True, metavar="LABEL=PATH",
                   help="reference run; deltas are computed against this one")
    p.add_argument("--compare", action="append", default=[], metavar="LABEL=PATH",
                   help="run to diff against the baseline (repeatable)")
    p.add_argument("--kind", choices=list(KINDS), default=None,
                   help="CSV shape (default: auto-detect from the baseline header)")
    p.add_argument("--stack", default="dada2-rs",
                   help="summary.csv only: which stack's rows to compare "
                        "(dada2-rs, R-split, R-single). Default dada2-rs")
    p.add_argument("--metric", action="append", default=None,
                   help="metric(s) to show (repeatable); default: all for the kind")
    p.add_argument("--json", metavar="OUT", help="also write the diff as JSON")
    args = p.parse_args()

    if not args.compare:
        sys.exit("need at least one --compare LABEL=PATH")

    # Resolve the baseline first so we can auto-detect the kind from its header.
    kind = args.kind
    base_label, base_path = parse_run(
        args.baseline, KINDS[kind]["file"] if kind else "summary.csv")
    if kind is None:
        kind = detect_kind(base_path)
        # Re-resolve in case the dir holds a non-default file for this kind.
        base_label, base_path = parse_run(args.baseline, KINDS[kind]["file"])

    spec = KINDS[kind]
    metrics = args.metric or spec["metrics"]
    for m in metrics:
        if m not in spec["metrics"]:
            sys.exit(f"metric {m!r} not in {kind} CSV (have: {', '.join(spec['metrics'])})")

    runs, key_order = [], []
    for label, path in [(base_label, base_path)] + [
            parse_run(c, spec["file"]) for c in args.compare]:
        rows, order = load(path, kind, args.stack)
        if not rows:
            where = f" (stack={args.stack})" if kind == "summary" else ""
            sys.exit(f"{label}: no rows in {path}{where}")
        runs.append((label, rows))
        for k in order:
            if k not in key_order:
                key_order.append(k)

    label_desc = ", ".join(l for l, _ in runs[1:])
    print("=" * 64)
    print(f"BENCH COMPARE — {kind}"
          + (f" (stack={args.stack})" if kind == "summary" else "")
          + f"\nbaseline: {base_label}   vs: {label_desc}")
    print("=" * 64)
    for m in metrics:
        print_table(m, key_order, runs, base_label)

    if args.json:
        out = {"kind": kind, "stack": args.stack if kind == "summary" else None,
               "baseline": base_label, "metrics": metrics,
               "runs": {label: rows for label, rows in runs}}
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f"\nWrote {args.json}")


if __name__ == "__main__":
    main()
