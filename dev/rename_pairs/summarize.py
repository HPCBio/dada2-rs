#!/usr/bin/env python3
"""How far each rename pair's losing sibling is from budding (issue #246).

Reads find_pairs.py's table and trace_pairs.sh's per-sample runs. For each
pair, in the baseline arm and in the arm it flipped in:

- reproduced: the trace run reports the same name as the floor run did
  (baseline: base_seq is an ASV and arm_seq is not; arm: the reverse). A pair
  that is not reproduced was traced under different inputs or flags, and its
  numbers say nothing about the floor run.
- the sibling: whichever of the pair is not its cluster's centre. Its reads,
  e_reads, and the budding p-value pA, and how many orders of magnitude
  pA x nraw sits above OMEGA_A ("short"; <= 0 means it would bud).

pA is recomputed as the budding test computes it: `calc_pA` conditioned on
the sequence being present, P(X >= reads | E) / (1 - e^-E). It is NOT the
trace's `pval`, which is the final OMEGA_C pass's unconditioned value
(prior = true, as in R). For a doubleton that is smaller by about a factor of
E, and an earlier trace read it as the budding p. The trace's `pval` is kept
as a check: recomputed unconditioned from reads and e_reads, it must agree,
and disagreements are counted.

Assumes no --prior-file and no singleton detection: a raw with a prior buds on
the unconditioned value.

Usage:
  summarize.py <pairs.tsv> <trace-dir> [--base insertion] [-o summary.tsv]
"""
import argparse
import csv
import json
import math
import os
import sys

COLUMNS = ["sample", "base_id", "arm_id", "arm", "edit", "base_reads", "arm_reads",
           "reproduced", "same_cluster", "cluster_reads",
           "sib_reads_base", "e_reads_base", "pA_nraw_base", "short_base",
           "sib_reads_arm", "e_reads_arm", "pA_nraw_arm", "short_arm", "short_min"]


def calc_pa(n, e, prior=False):
    """calc_pA: P(X >= n | E), conditioned on X >= 1 unless prior."""
    if e <= 0:
        return 0.0
    if e < n:
        tail = sum(math.exp(-e + k * math.log(e) - math.lgamma(k + 1)) for k in range(n, n + 200))
    else:
        tail = 1 - sum(math.exp(-e + k * math.log(e) - math.lgamma(k + 1)) for k in range(n))
    if prior:
        return tail
    norm = 1 - math.exp(-e)
    if norm < 1e-7:  # calc_pA's TAIL_APPROX_CUTOFF
        norm = e - 0.5 * e * e
    return tail / norm


class Run:
    """One (sample, arm) trace run: its ASVs and its cluster trace."""

    def __init__(self, d, tag):
        self.asvs = {a["sequence"].upper() for a in json.load(open(f"{d}/{tag}.json"))["asvs"]}
        t = json.load(open(f"{d}/{tag}.trace.json"))
        if t.get("trace_no_members") or t.get("trace_min_abund", 1) > 1:
            sys.exit(f"ERROR: {d}/{tag}.trace.json has incomplete members; re-run without trace filters")
        self.nraw, self.omega_a = t["nraw"], t["omega_a"]
        seqs = [s.upper() for s in t["sequences"]]
        self.where = {}  # seq -> (cluster, member record)
        for c in t["clusters"]:
            for m in c.get("members") or []:
                self.where[seqs[m["raw_seq_id"]]] = (c, m)
        self.centre = {c["id"]: seqs[c["center_seq_id"]] for c in t["clusters"]}


def sibling(run, a, b):
    """The non-centre member of {a, b} in run, with its numbers; None if both or neither is a centre."""
    ca, ma = run.where.get(a, (None, None))
    cb, mb = run.where.get(b, (None, None))
    if ca is None or cb is None:
        return None
    same = ca["id"] == cb["id"]
    for seq, c, m in ((a, ca, ma), (b, cb, mb)):
        if run.centre[c["id"]] != seq:
            p = calc_pa(m["abundance"], m["e_reads"])
            short = math.inf if p >= 1 else (-math.inf if p <= 0 else
                                               math.log10(p * run.nraw / run.omega_a))
            return {"same": same, "cluster_reads": c["abundance"], "reads": m["abundance"],
                    "e_reads": m["e_reads"], "pval": m["pval"], "pA_nraw": p * run.nraw, "short": short}
    return None


def fmt(x, spec):
    return "" if x is None else format(x, spec)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pairs")
    ap.add_argument("trace_dir")
    ap.add_argument("--base", default="insertion", help="baseline arm directory name")
    ap.add_argument("-o", "--output")
    a = ap.parse_args()

    runs, rows, p_mismatch = {}, [], 0

    def run(smp, tag):
        if (smp, tag) not in runs:
            runs[smp, tag] = Run(os.path.join(a.trace_dir, smp), tag)
        return runs[smp, tag]

    for p in csv.DictReader(open(a.pairs), delimiter="\t"):
        s, t, arm = p["base_seq"], p["arm_seq"], p["arms"].split(",")[0]
        rb, ra = run(p["sample"], a.base), run(p["sample"], arm)
        reproduced = s in rb.asvs and t not in rb.asvs and t in ra.asvs and s not in ra.asvs
        sb, sa = sibling(rb, s, t), sibling(ra, s, t)
        for sib in (sb, sa):
            if sib and sib["e_reads"] > 0:
                ref = calc_pa(sib["reads"], sib["e_reads"], prior=True)
                if not math.isclose(sib["pval"], ref, rel_tol=1e-6, abs_tol=1e-300):
                    p_mismatch += 1
        shorts = [x["short"] for x in (sb, sa) if x]
        rows.append({
            "sample": p["sample"], "base_id": p["base_id"], "arm_id": p["arm_id"], "arm": arm,
            "edit": p["edit"], "base_reads": p["base_reads"], "arm_reads": p["arm_reads"],
            "reproduced": "yes" if reproduced else "NO",
            "same_cluster": "yes" if sb and sa and sb["same"] and sa["same"] else "NO",
            "cluster_reads": fmt(sb and sb["cluster_reads"], "d"),
            "sib_reads_base": fmt(sb and sb["reads"], "d"),
            "e_reads_base": fmt(sb and sb["e_reads"], ".2e"),
            "pA_nraw_base": fmt(sb and sb["pA_nraw"], ".1e"),
            "short_base": fmt(sb and sb["short"], ".1f"),
            "sib_reads_arm": fmt(sa and sa["reads"], "d"),
            "e_reads_arm": fmt(sa and sa["e_reads"], ".2e"),
            "pA_nraw_arm": fmt(sa and sa["pA_nraw"], ".1e"),
            "short_arm": fmt(sa and sa["short"], ".1f"),
            "short_min": fmt(min(shorts) if shorts else None, ".1f"),
        })

    out = open(a.output, "w") if a.output else sys.stdout
    print("\t".join(COLUMNS), file=out)
    for r in rows:
        print("\t".join(r[c] for c in COLUMNS), file=out)

    # The summary goes to stderr so a -o-less run still yields a clean table.
    err = sys.stderr
    ok = [r for r in rows if r["reproduced"] == "yes" and r["same_cluster"] == "yes"]
    omega = next(iter(runs.values())).omega_a if runs else float("nan")
    print(f"\n{len(rows)} pair(s); reproduced in one cluster in both arms: {len(ok)}", file=err)
    for label, bad in (("not reproduced", lambda r: r["reproduced"] != "yes"),
                       ("not in one cluster", lambda r: r["reproduced"] == "yes" and r["same_cluster"] != "yes")):
        n = sum(map(bad, rows))
        if n:
            print(f"  {label}: {n} (excluded below)", file=err)
    if p_mismatch:
        print(f"  WARNING: {p_mismatch} trace pval(s) disagree with the unconditioned calc_pA "
              "recomputed from reads and e_reads; the trace or this script is wrong", file=err)
    bins = [(-math.inf, 0, "would bud (<= 0)"), (0, 1, "0-1"), (1, 3, "1-3"), (3, 6, "3-6"),
            (6, 10, "6-10"), (10, 20, "10-20"), (20, math.inf, "> 20")]
    print(f"\norders of magnitude short of OMEGA_A = {omega:.0e}, closer sibling of each pair:", file=err)
    for lo, hi, label in bins:
        sel = [r for r in ok if lo < float(r["short_min"]) <= hi or (hi == 0 and float(r["short_min"]) <= 0)]
        if sel:
            reads = sorted(max(int(r["base_reads"]), int(r["arm_reads"])) for r in sel)
            print(f"  {label:>16}: {len(sel):4} pair(s), ASV reads {reads[0]}-{reads[-1]}", file=err)


if __name__ == "__main__":
    main()
