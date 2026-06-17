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
  branch, k5, j4). PATH may also be a comma-separated list of CSVs/dirs — repeated
  runs (reps) of the same config — which are reduced to the median per metric with
  the spread shown (`median ±N%`), so within-node noise is visible at a glance.

Examples:
  # peak RSS + wall, main build vs the deferred-division branch (issue #23)
  compare_bench.py --baseline main=bench_main --compare branch=bench_branch

  # reps: 3 runs per config, reported as median ±half-range
  compare_bench.py \\
      --baseline main=run1_main,run2_main,run3_main \\
      --compare branch=run1_br,run2_br,run3_br

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
import statistics
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
    """'label=path[,path...]' -> (label, [csv_path, ...]). Each path may be a CSV
    file or a directory containing the kind's CSV. Multiple comma-separated paths
    are repeated runs (reps) of the same configuration — reduced to the median per
    metric, with the spread reported so within-node noise is visible."""
    if "=" not in spec:
        sys.exit(f"--baseline/--compare expects LABEL=PATH[,PATH...], got: {spec!r}")
    label, raw = spec.split("=", 1)
    paths = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        p = Path(item)
        if p.is_dir():
            p = p / default_file
        if not p.is_file():
            sys.exit(f"{label}: not a file: {p}")
        paths.append(p)
    if not paths:
        sys.exit(f"{label}: no paths given")
    return label, paths


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


def load_reps(paths, kind, stack):
    """Load N rep CSVs for one label and reduce each metric to the MEDIAN across
    reps. Returns (rows, order, spread, nreps): rows[key][metric] = median, and
    spread[key][metric] = relative half-range ((max-min)/2/median) or None (single
    rep / zero median). Keys missing from a rep are simply absent from its median."""
    per_rep, order = [], []
    for path in paths:
        rows, od = load(path, kind, stack)
        per_rep.append(rows)
        for k in od:
            if k not in order:
                order.append(k)
    spec = KINDS[kind]
    rows, spread = {}, {}
    for k in order:
        rows[k], spread[k] = {}, {}
        for m in spec["metrics"]:
            vals = [r[k][m] for r in per_rep if k in r and r[k].get(m) is not None]
            if not vals:
                rows[k][m] = spread[k][m] = None
                continue
            med = statistics.median(vals)
            rows[k][m] = med
            spread[k][m] = ((max(vals) - min(vals)) / 2 / med
                            if len(vals) > 1 and med else None)
    return rows, order, spread, len(paths)


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
    """runs = list of (label, rows, spread, nreps); base_label is the baseline's.
    With reps, each cell shows `median ±spread%` and Δ% compares medians."""
    arrow = "↓ better" if metric in LOWER_IS_BETTER else "↑ better"
    print(f"\n### {metric}  ({arrow})")
    any_reps = any(nr > 1 for *_, nr in runs)
    vw = 18 if any_reps else 12   # wider value column when ±spread is shown
    klen = max([len("step")] + [len(k) for k in keys])

    def cell(rows, spread, k):
        v = rows.get(k, {}).get(metric)
        s = (spread.get(k) or {}).get(metric)
        txt = fmt(metric, v)
        if s is not None:
            txt += f" ±{s * 100:.0f}%"
        return txt, v

    head = f"{'':<{klen}}  {base_label:>{vw}}"
    for label, *_ in runs[1:]:
        head += f"  {label:>{vw}}{'Δ%':>8}"
    print(head)
    _, base_rows, base_spread, _ = runs[0]
    for k in keys:
        btxt, base_v = cell(base_rows, base_spread, k)
        line = f"{k:<{klen}}  {btxt:>{vw}}"
        for _, rows, spread, _ in runs[1:]:
            txt, v = cell(rows, spread, k)
            d = pct(base_v, v)
            dstr = "—" if d is None else f"{d:+.1f}%"
            line += f"  {txt:>{vw}}{dstr:>8}"
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
    base_label, base_paths = parse_run(
        args.baseline, KINDS[kind]["file"] if kind else "summary.csv")
    if kind is None:
        kind = detect_kind(base_paths[0])
        # Re-resolve in case the dir holds a non-default file for this kind.
        base_label, base_paths = parse_run(args.baseline, KINDS[kind]["file"])

    spec = KINDS[kind]
    metrics = args.metric or spec["metrics"]
    for m in metrics:
        if m not in spec["metrics"]:
            sys.exit(f"metric {m!r} not in {kind} CSV (have: {', '.join(spec['metrics'])})")

    runs, key_order = [], []
    for label, paths in [(base_label, base_paths)] + [
            parse_run(c, spec["file"]) for c in args.compare]:
        rows, order, spread, nreps = load_reps(paths, kind, args.stack)
        if not rows:
            where = f" (stack={args.stack})" if kind == "summary" else ""
            sys.exit(f"{label}: no rows in {', '.join(map(str, paths))}{where}")
        runs.append((label, rows, spread, nreps))
        for k in order:
            if k not in key_order:
                key_order.append(k)

    def label_n(label, nreps):
        return f"{label} (n={nreps})" if nreps > 1 else label

    label_desc = ", ".join(label_n(l, nr) for l, _, _, nr in runs[1:])
    print("=" * 64)
    print(f"BENCH COMPARE — {kind}"
          + (f" (stack={args.stack})" if kind == "summary" else "")
          + f"\nbaseline: {label_n(base_label, runs[0][3])}   vs: {label_desc}")
    if any(nr > 1 for *_, nr in runs):
        print("(values are medians across reps; ±N% = relative half-range)")
    print("=" * 64)
    for m in metrics:
        print_table(m, key_order, runs, base_label)

    if args.json:
        out = {"kind": kind, "stack": args.stack if kind == "summary" else None,
               "baseline": base_label, "metrics": metrics,
               "runs": {label: {"nreps": nr, "median": rows, "spread": spread}
                        for label, rows, spread, nr in runs}}
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f"\nWrote {args.json}")


if __name__ == "__main__":
    main()
