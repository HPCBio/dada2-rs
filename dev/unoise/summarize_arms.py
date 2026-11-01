#!/usr/bin/env python3
"""Summarize a run_arms.sh grid (#323): one row per arm, scored against a baseline.

For each finished arm under <out-dir> it reports, from the chimera-free table:

- ASVs before and after chimera removal, and reads kept;
- ASV set identity against the baseline arm: shared, only in the arm, only in
  the baseline, and Jaccard. Counts of ASVs alone are not compared, because two
  arms can call the same number of different sequences;
- Bray-Curtis dissimilarity of per-ASV read totals against the baseline, which
  catches arms that call the same set but move reads (--reassign can do only
  that);
- the summed dada phase timings from --metrics-json;
- with refeval.json present (REFERENCE set): TP, near, FP, references recovered,
  precision and recall.

--pair A:B adds a set-identity row for any two arms, e.g. sp_indel:unoise, the
comparison of the learned error model with the fixed rule when both count indels.

Usage: dev/unoise/summarize_arms.py <out-dir> [--baseline divisive]
           [--pair A:B ...] [--tsv out.tsv]
"""
import argparse
import csv
import json
import os
import sys

PHASES = ["compare", "shuffle", "bud", "p_update", "index_add"]


def load_table(path):
    """Per-sequence read totals (summed over samples) from a seqtab JSON."""
    with open(path) as fh:
        j = json.load(fh)
    seqs = [s.upper() for s in j["sequences"]]
    totals = [0] * len(seqs)
    for row in j["counts"]:
        for i, c in enumerate(row):
            totals[i] += c
    return dict(zip(seqs, totals))


def phase_seconds(path):
    """Sum of the dada phase timers over every run in a metrics JSON."""
    if not os.path.exists(path):
        return None

    def runs(o):
        if isinstance(o, dict):
            if isinstance(o.get("phases"), dict) and "shuffle" in o["phases"]:
                yield o
            for v in o.values():
                yield from runs(v)
        elif isinstance(o, list):
            for v in o:
                yield from runs(v)

    with open(path) as fh:
        rows = list(runs(json.load(fh)))
    return sum(r["phases"].get(k) or 0.0 for r in rows for k in PHASES)


def identity(a, b):
    sa, sb = set(a), set(b)
    shared = len(sa & sb)
    union = len(sa | sb)
    keys = sa | sb
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.values()) + sum(b.values())
    return {
        "shared": shared,
        "only_arm": len(sa - sb),
        "only_base": len(sb - sa),
        "jaccard": shared / union if union else 1.0,
        "bray_curtis": num / den if den else 0.0,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir")
    ap.add_argument("--baseline", default="divisive")
    ap.add_argument("--pair", action="append", default=[], metavar="A:B")
    ap.add_argument("--tsv", help="also write the arm table here")
    args = ap.parse_args()

    arms = sorted(
        d for d in os.listdir(args.out_dir)
        if os.path.exists(os.path.join(args.out_dir, d, ".done"))
    )
    if not arms:
        sys.exit(f"{args.out_dir}: no finished arms (no <arm>/.done)")
    tables = {a: load_table(os.path.join(args.out_dir, a, "nochim.json")) for a in arms}
    pre = {a: len(load_table(os.path.join(args.out_dir, a, "seqtab.json"))) for a in arms}
    base = tables.get(args.baseline)
    if base is None:
        print(f"note: baseline arm '{args.baseline}' not finished; no identity columns")

    rows = []
    for a in arms:
        o = os.path.join(args.out_dir, a)
        secs = [phase_seconds(os.path.join(o, f"{d}.metrics.json")) for d in ("fwd", "rev")]
        secs = [s for s in secs if s is not None]
        r = {
            "arm": a,
            "asvs_prechim": pre[a],
            "asvs": len(tables[a]),
            "reads": sum(tables[a].values()),
            "dada_s": f"{sum(secs):.2f}" if secs else "",
        }
        if base is not None:
            r.update(identity(tables[a], base))
        ref = os.path.join(o, "refeval.json")
        if os.path.exists(ref):
            with open(ref) as fh:
                e = json.load(fh)
            for k in ("tp", "tp_exact", "near", "fp", "recovered_refs", "n_ref",
                      "precision", "recall"):
                r[k] = e.get(k, "")
        rows.append(r)

    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]

    def fmt(v):
        return f"{v:.4f}" if isinstance(v, float) else str(v)

    widths = {c: max(len(c), *(len(fmt(r.get(c, ""))) for r in rows)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    for r in rows:
        print("  ".join(fmt(r.get(c, "")).ljust(widths[c]) for c in cols))
    if args.tsv:
        with open(args.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    for p in args.pair:
        a, _, b = p.partition(":")
        if a not in tables or b not in tables:
            print(f"pair {p}: arm not finished")
            continue
        i = identity(tables[a], tables[b])
        print(
            f"pair {a} vs {b}: shared {i['shared']}, only {a} {i['only_arm']}, "
            f"only {b} {i['only_base']}, Jaccard {i['jaccard']:.4f}, "
            f"Bray-Curtis {i['bray_curtis']:.4f}"
        )


if __name__ == "__main__":
    main()
