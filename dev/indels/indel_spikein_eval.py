#!/usr/bin/env python3
"""indel_spikein_eval.py — score dada output on indel_spikein.py samples (issue #291).

For each synthetic sample, the truth is a base truth set (e.g. a mock's
reference alleles, or its denoised ASVs) plus that sample's variants.
`dada2-rs reference-eval` classifies every ASV against it; this script joins
the per-reference recovery to `truth.tsv` and reports:

  - plain variants recovered, by size, position, type and reads moved;
  - homopolymer controls recovered (the target is 0: they must not split);
  - false positives (ASVs beyond --near-diffs of every truth sequence).

A variant counts as recovered only when an ASV matches it exactly.

Usage:
  indel_spikein_eval.py --spikes spikes/bei --dada-dir spikes/bei/dada_default \\
      --base-truth truth/bei/bei_union_27F-1492R.fasta --out spikes/bei/eval_default
"""

import argparse
import collections
import csv
import os
import subprocess
import sys


def read_fasta(path):
    out, name = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            out.append([name, ""])
        elif name:
            out[-1][1] += line
    return out


def bins(n):
    for hi, label in ((1, "1"), (4, "2-4"), (19, "5-19"), (99, "20-99"), (499, "100-499")):
        if n <= hi:
            return label
    return ">=500"


def table(rows, key, order=None):
    g = collections.defaultdict(list)
    for r in rows:
        g[key(r)].append(r["recovered"])
    keys = order if order else sorted(g)
    return "  ".join(f"{k}:{sum(g[k])}/{len(g[k])}" for k in keys if k in g)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--spikes", required=True, help="indel_spikein.py --out-dir")
    ap.add_argument("--dada-dir", required=True, help="dada --output-dir over the spike FASTQs")
    ap.add_argument("--base-truth", required=True, help="FASTA of true sequences without spikes")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bin", default="dada2-rs")
    ap.add_argument("--near-diffs", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    truth = list(csv.DictReader(open(os.path.join(a.spikes, "truth.tsv")), delimiter="\t"))
    base = read_fasta(a.base_truth)
    os.makedirs(a.out, exist_ok=True)
    recovered, fp = {}, collections.Counter()
    for sample in sorted({t["sample"] for t in truth}):
        ref = os.path.join(a.out, f"{sample}.truth.fasta")
        with open(ref, "w") as fh:
            for name, seq in base:
                fh.write(f">{name}\n{seq}\n")
            for t in truth:
                if t["sample"] == sample:
                    fh.write(f">{t['variant']}\n{t['sequence']}\n")
        per_ref = os.path.join(a.out, f"{sample}.per_ref.tsv")
        per_asv = os.path.join(a.out, f"{sample}.per_asv.tsv")
        subprocess.run(
            [a.bin, "reference-eval", "--asvs", os.path.join(a.dada_dir, f"{sample}.json"),
             "--reference", ref, "--band", "-1", "--threads", str(a.threads),
             "--near-diffs", str(a.near_diffs), "--per-ref", per_ref, "--per-asv", per_asv,
             "-o", os.path.join(a.out, f"{sample}.json")],
            check=True, stderr=subprocess.DEVNULL,
        )
        for r in csv.DictReader(open(per_ref), delimiter="\t"):
            recovered[(sample, r["ref_id"])] = r["recovered"] == "true"
        for r in csv.DictReader(open(per_asv), delimiter="\t"):
            if r["class"] == "FP":
                fp[sample] += 1

    for t in truth:
        t["recovered"] = recovered[(t["sample"], t["variant"])]
    plain = [t for t in truth if t["context"] == "plain"]
    hp = [t for t in truth if t["context"] == "hp"]
    moved = ["1", "2-4", "5-19", "20-99", "100-499", ">=500"]
    print(f"plain variants recovered: {sum(t['recovered'] for t in plain)}/{len(plain)}")
    print("  by size:    ", table(plain, lambda t: int(t["size"])))
    print("  by position:", table(plain, lambda t: t["position"], ["5p", "mid", "3p"]))
    print("  by type:    ", table(plain, lambda t: t["type"]))
    print("  by reads moved:", table(plain, lambda t: bins(int(t["reads_moved"])), moved))
    small = [t for t in plain if int(t["size"]) <= 10]
    print("  size<=10 by reads moved:", table(small, lambda t: bins(int(t["reads_moved"])), moved))
    print(f"homopolymer controls recovered (target 0): {sum(t['recovered'] for t in hp)}/{len(hp)}")
    print(f"false positives (> {a.near_diffs} edits from all truth): {sum(fp.values())}")
    with open(os.path.join(a.out, "variants.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(truth)


if __name__ == "__main__":
    sys.exit(main())
