#!/usr/bin/env python3
"""Summarise a run_endmodel_sweep_289.sh output directory (#289, #51).

For each arm (<backend>.<end-model>), takes the median over replicates of the
step wall times, learn-errors rounds and seconds per round, and compares the
first replicate's outputs with the baseline arm (default nw.ends-free, R's
model): ASVs gained and lost, table cells that differ, reads moved, and the
change in error-model substitution counts.

Usage: compare_endmodel_arms.py <sweep-dir> [--baseline nw.ends-free] [--tsv out.tsv]
"""
import argparse
import csv
import glob
import json
import os
import statistics


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def step_times(rep_dir):
    path = os.path.join(rep_dir, "summary.csv")
    rows = [l for l in open(path) if not l.startswith("#")]
    out = {}
    for r in csv.DictReader(rows):
        if r["stack"] == "dada2-rs":
            out[r["step"]] = float(r["wall_s"])
    return out


def err_models(rep_dir):
    return sorted(glob.glob(os.path.join(rep_dir, "rust", "errors_*.json")))


def table(rep_dir):
    d = load_json(os.path.join(rep_dir, "rust", "seqtab_nochim.json"))
    cells = {}
    for i, s in enumerate(d["samples"]):
        for j, q in enumerate(d["sequences"]):
            c = d["counts"][i][j]
            if c:
                cells[(s, q)] = c
    return cells, set(d["sequences"])


def substitutions(err_path):
    t = load_json(err_path)["trans"]
    # trans rows are A2A, A2C, ... T2T; every fifth row (0, 5, 10, 15) is a match.
    return sum(sum(t[i]) for i in range(16) if i % 5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep")
    ap.add_argument("--baseline", default="nw.ends-free")
    ap.add_argument("--tsv")
    a = ap.parse_args()

    arms = sorted(d for d in os.listdir(a.sweep) if os.path.isdir(os.path.join(a.sweep, d)))
    if a.baseline not in arms:
        raise SystemExit(f"baseline arm {a.baseline} not found in {a.sweep}")
    base_dir = os.path.join(a.sweep, a.baseline, "rep1")
    base_cells, base_seqs = table(base_dir)
    base_subs = [substitutions(p) for p in err_models(base_dir)]

    rows = []
    for arm in arms:
        reps = sorted(glob.glob(os.path.join(a.sweep, arm, "rep*")))
        reps = [r for r in reps if os.path.exists(os.path.join(r, "summary.csv"))]
        if not reps:
            continue
        times = [step_times(r) for r in reps]
        steps = sorted({k for t in times for k in t})
        med = {k: statistics.median(t[k] for t in times if k in t) for k in steps}
        rounds = [load_json(p).get("iterations") for p in err_models(reps[0])]
        learn_steps = [k for k in steps if k.startswith("learn")]
        learn_s = sum(med[k] for k in learn_steps)
        per_round = learn_s / sum(rounds) if all(rounds) and rounds else float("nan")

        cells, seqs = table(reps[0])
        keys = set(base_cells) | set(cells)
        moved = sum(abs(base_cells.get(k, 0) - cells.get(k, 0)) for k in keys) / 2
        subs = [substitutions(p) for p in err_models(reps[0])]
        dsubs = ";".join(f"{(s - b) / b * 100:+.2f}%" for s, b in zip(subs, base_subs))
        rows.append({
            "arm": arm,
            "reps": len(reps),
            "total_s": round(med.get("TOTAL", float("nan")), 1),
            "learn_s": round(learn_s, 1),
            "rounds": "+".join(str(r) for r in rounds),
            "s_per_round": round(per_round, 2),
            "dada_s": round(sum(v for k, v in med.items() if k.startswith("dada")), 1),
            "asvs": len(seqs),
            "gained": len(seqs - base_seqs),
            "lost": len(base_seqs - seqs),
            "cells_differ": sum(base_cells.get(k, 0) != cells.get(k, 0) for k in keys),
            "reads_moved": int(moved),
            "reads": sum(cells.values()),
            "err_subs_vs_base": dsubs,
        })

    cols = list(rows[0])
    widths = {c: max(len(c), *(len(str(r[c])) for r in rows)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    for r in rows:
        print("  ".join(str(r[c]).ljust(widths[c]) for c in cols))
    print(f"\nbaseline: {a.baseline} rep1; times are medians over reps; "
          "table and error-model columns compare rep1.")
    if a.tsv:
        with open(a.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
