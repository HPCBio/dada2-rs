#!/usr/bin/env python3
"""Summarise a DADA2RS_MEMBER_ORDER experiment (issue #157).

Each arm directory holds one run's per-sample dada JSON (dada, dada-pseudo or
dada-pooled output). Per arm, ASVs are summed over samples, then compared with
the baseline arm and pairwise among the non-baseline arms.

The shuffle arms form the null: their pairwise churn is the spread of the
member-order floor. The sorted arm is placed against it.

If an arm directory has the run's --verbose log as log.txt, birth order is
compared too: the sequence of raws that bud. It is more sensitive than the
final table, since two tied births can swap and still produce the same ASVs.

Usage:
  compare_member_order.py <baseline-dir> <arm-dir> [<arm-dir> ...]
Arm names are the directory basenames.
"""
import itertools
import re
import json
import os
import statistics
import sys


def load(arm_dir):
    """seq -> total abundance over all per-sample JSONs in arm_dir."""
    tot = {}
    for f in sorted(os.listdir(arm_dir)):
        if not f.endswith(".json") or f.startswith("_"):
            continue
        d = json.load(open(os.path.join(arm_dir, f)))
        if "asvs" not in d:
            continue
        for a in d["asvs"]:
            s = a["sequence"].upper()
            tot[s] = tot.get(s, 0) + a["abundance"]
    return tot


def births(arm_dir):
    """Raw indices in the order they budded, from the run's --verbose log."""
    p = os.path.join(arm_dir, "log.txt")
    if not os.path.exists(p):
        return None
    return re.findall(r"Division \((?:naive|prior)\): Raw (\d+)", open(p).read())


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else None


def nearest(seq, pool):
    """Smallest equal-length Hamming distance from seq to any sequence in pool."""
    ds = [h for p in pool if (h := hamming(seq, p)) is not None]
    return min(ds) if ds else None


def compare(a, b):
    """Churn between two arms' pooled tables."""
    only_a, only_b = set(a) - set(b), set(b) - set(a)
    moved = sum(abs(a.get(s, 0) - b.get(s, 0)) for s in set(a) | set(b)) // 2
    near = [nearest(s, set(a)) for s in only_b]
    return {
        "n_a": len(a), "n_b": len(b),
        "only_a": len(only_a), "only_b": len(only_b),
        "churn": len(only_a) + len(only_b),
        "reads_moved": moved,
        "reads_total": sum(a.values()),
        "only_b_max_abund": max((b[s] for s in only_b), default=0),
        "only_b_h1": sum(1 for h in near if h == 1),
        "only_b_h_le2": sum(1 for h in near if h is not None and h <= 2),
    }


def main():
    base_dir, arm_dirs = sys.argv[1], sys.argv[2:]
    base = load(base_dir)
    arms = {os.path.basename(d.rstrip("/")): load(d) for d in arm_dirs}
    name = os.path.basename(base_dir.rstrip("/"))

    print(f"baseline {name}: {len(base)} ASVs, {sum(base.values())} reads\n")
    print(f"{'arm':12} {'ASVs':>5} {'churn':>6} {'-base':>6} {'+arm':>5} "
          f"{'reads moved':>12} {'max new abund':>14} {'new H1':>7} {'new H<=2':>8}")
    for arm, t in arms.items():
        c = compare(base, t)
        print(f"{arm:12} {c['n_b']:5} {c['churn']:6} {c['only_a']:6} {c['only_b']:5} "
              f"{c['reads_moved']:7} ({100 * c['reads_moved'] / c['reads_total']:.3f}%) "
              f"{c['only_b_max_abund']:14} {c['only_b_h1']:7} {c['only_b_h_le2']:8}")

    b0 = births(base_dir)
    if b0 is not None:
        print(f"\nbirth order ({len(b0)} births in baseline):")
        for arm, d in zip(arms, arm_dirs):
            b = births(d)
            if b is None:
                continue
            same_set = sorted(b) == sorted(b0)
            moved = sum(x != y for x, y in zip(b0, b)) + abs(len(b) - len(b0))
            print(f"  {arm:12} same set {same_set}, positions differing {moved}")

    shuffles = [a for a in arms if a.startswith("shuffle")]
    if len(shuffles) >= 2:
        pair = [compare(arms[x], arms[y])["churn"] for x, y in itertools.combinations(shuffles, 2)]
        vs_base = [compare(base, arms[s])["churn"] for s in shuffles]
        print(f"\nnull (shuffle arms): churn vs baseline {vs_base}, "
              f"pairwise median {statistics.median(pair)} range {min(pair)}-{max(pair)}")
        if "sorted" in arms:
            vs_sorted = [compare(arms["sorted"], arms[s])["churn"] for s in shuffles]
            print(f"sorted vs baseline {compare(base, arms['sorted'])['churn']}; "
                  f"sorted vs each shuffle {vs_sorted}")


if __name__ == "__main__":
    main()
