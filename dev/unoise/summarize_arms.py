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
- the summed dada phase timings, and the dada step's wall time and peak memory,
  from --metrics-json. The vsearch arm has no phases; its wall time and peak
  memory are those of --cluster_unoise, from time_fwd.txt (bash's `time`, used
  where GNU time is missing, records no memory);
- with refeval.json present (REFERENCE set): TP, near, FP, references recovered,
  precision and recall.

--sample NAME scores one sample instead: the refeval columns come from
<arm>/sample_NAME/, written by dev/unoise/eval_sample.sh (e.g. MiSeq SOP's
Mock), with the sample's ASVs and reads, the reads in FP ASVs (in a pooled run
mostly 1-2 read ASVs bled in from other samples, so fp alone overstates the
error) and tp_chim_lost, the TP ASVs that chimera removal dropped.

--pair A:B adds a set-identity row for any two arms, e.g. sp_indel:unoise, the
comparison of the learned error model with the fixed rule when both count indels.

Usage: dev/unoise/summarize_arms.py <out-dir> [--baseline divisive]
           [--sample NAME] [--pair A:B ...] [--tsv out.tsv]
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


def run_cost(path):
    """Wall seconds and peak RSS (MB) of one dada run, from its metrics JSON."""
    if not os.path.exists(path):
        return None, None
    with open(path) as fh:
        j = json.load(fh)
    peaks = (j.get("pipeline") or {}).get("peak_rss_mb") or {}
    return j.get("wall_seconds"), max(peaks.values(), default=None)


def time_file(path):
    """Wall seconds and peak RSS (MB) from a time_<dir>.txt, whichever timer
    wrote it: GNU time -v, BSD time -l or bash's `time` (no memory)."""
    if not os.path.exists(path):
        return None, None
    wall = peak = None
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if t.startswith("Elapsed (wall clock)"):
                secs = 0.0
                for part in t.rsplit(" ", 1)[1].split(":"):
                    secs = secs * 60 + float(part)
                wall = secs
            elif t.startswith("Maximum resident set size (kbytes)"):
                peak = int(t.rsplit(" ", 1)[1]) / 1024
            elif t.endswith("maximum resident set size"):
                peak = int(t.split()[0]) / 2**20
            elif t.split()[1:2] == ["real"]:
                wall = float(t.split()[0])
            elif t.startswith("bash time:"):
                f = t.split()
                wall = float(f[f.index("wall_seconds") + 1])
    return wall, peak


def fp_reads(path):
    """Reads in FP ASVs, from the ;size= in a per-ASV TSV's query_id."""
    n = 0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["class"] == "FP":
                n += next((int(f[5:]) for f in r["query_id"].split(";")
                           if f.startswith("size=")), 0)
    return n


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
    ap.add_argument("--sample", help="score this sample (eval_sample.sh output)")
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
        # The dada step's wall time and peak memory (the vsearch arm: its
        # --cluster_unoise step, from the time file).
        costs = [run_cost(os.path.join(o, f"{d}.metrics.json")) for d in ("fwd", "rev")]
        costs = [c for c in costs if c[0] is not None]
        if not costs:
            costs = [c for c in (time_file(os.path.join(o, "time_fwd.txt")),) if c[0] is not None]
        walls = [w for w, _ in costs]
        peaks = [p for _, p in costs if p is not None]
        r = {
            "arm": a,
            "asvs_prechim": pre[a],
            "asvs": len(tables[a]),
            "reads": sum(tables[a].values()),
            "dada_s": f"{sum(secs):.2f}" if secs else "",
            "wall_s": f"{sum(walls):.2f}" if walls else "",
            "peak_mb": f"{max(peaks):.0f}" if peaks else "",
        }
        if base is not None:
            r.update(identity(tables[a], base))
        rd = os.path.join(o, f"sample_{args.sample}") if args.sample else o
        ref = os.path.join(rd, "refeval.json")
        if args.sample and not os.path.exists(ref):
            sys.exit(f"{ref} missing: run dev/unoise/eval_sample.sh first")
        if os.path.exists(ref):
            with open(ref) as fh:
                e = json.load(fh)
            if args.sample:
                with open(os.path.join(rd, "nochim.fa")) as fh:
                    sizes = [int(l.split("size=")[1]) for l in fh if l.startswith(">")]
                r["s_asvs"], r["s_reads"] = len(sizes), sum(sizes)
            for k in ("tp", "tp_exact", "near", "fp", "recovered_refs", "n_ref",
                      "precision", "recall"):
                r[k] = e.get(k, "")
            if args.sample:
                r["fp_reads"] = fp_reads(os.path.join(rd, "refeval_asv.tsv"))
                pj = os.path.join(rd, "refeval_prechim.json")
                if os.path.exists(pj):
                    with open(pj) as fh:
                        cb = json.load(fh).get("chimera_breakdown") or {}
                    r["tp_chim_lost"] = cb.get("tp_chimeric", "")
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
