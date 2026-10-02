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
Births are compared per partition: the whole log for dada-pooled, each sample's
run for dada (split on run_member_order.sh's "### sample" markers, or on each
run's leading [derep] line). dada-pseudo is skipped: its round-1 progress lines
carry no sample tag, so with concurrent samples they cannot be attributed.

Usage:
  compare_member_order.py <baseline-dir> <arm-dir> [<arm-dir> ...] [--id HASH ...]
Arm names are the directory basenames. --id (repeatable) takes a short id as
printed in the churned-ASV list, or any prefix of it, and prints that ASV's
sequence and its per-sample abundance in every arm.
"""
import hashlib
import itertools
import re
import json
import os
import statistics
import sys


def commands(arm_dir):
    """The set of dada2_rs_command tags among arm_dir's per-sample JSONs."""
    tags = set()
    for f in os.listdir(arm_dir):
        if f.endswith(".json") and not f.startswith("_"):
            tag = json.load(open(os.path.join(arm_dir, f))).get("dada2_rs_command")
            if tag:
                tags.add(tag)
    return tags


def load_per_sample(arm_dir):
    """sample file -> {seq: abundance} for each per-sample JSON in arm_dir."""
    out = {}
    for f in sorted(os.listdir(arm_dir)):
        if not f.endswith(".json") or f.startswith("_"):
            continue
        d = json.load(open(os.path.join(arm_dir, f)))
        if "asvs" in d:
            out[f[: -len(".json")]] = {a["sequence"].upper(): a["abundance"] for a in d["asvs"]}
    return out


def show_ids(ids, dirs, names):
    """Per-sample abundance of the ASVs whose short id starts with each of ids."""
    per = {n: load_per_sample(d) for n, d in zip(names, dirs)}
    seqs = {s for arm in per.values() for smp in arm.values() for s in smp}
    for want in ids:
        hits = sorted(s for s in seqs if short_id(s).startswith(want))
        if not hits:
            print(f"\n--id {want}: no ASV with that id in any arm")
            continue
        for seq in hits:
            print(f"\n--id {want} = {short_id(seq)}  length {len(seq)}\n{seq}")
            samples = sorted({smp for arm in per.values() for smp, t in arm.items() if seq in t})
            print(f"  {'sample':32} " + " ".join(f"{n:>10}" for n in names))
            for smp in samples:
                row = " ".join(f"{per[n].get(smp, {}).get(seq, 0):>10}" for n in names)
                print(f"  {smp:32} {row}")


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


BIRTH = re.compile(r"Division \((?:naive|prior)\): Raw (\d+)")


def births(arm_dir, command):
    """Per-partition lists of raw indices in budding order, or None.

    dada-pooled: one partition, the whole log. dada: one per sample run.
    dada-pseudo: None (not attributable; see the module docstring).
    """
    p = os.path.join(arm_dir, "log.txt")
    if not os.path.exists(p) or command == "dada-pseudo":
        return None
    text = open(p).read()
    if command == "dada-pooled":
        return [BIRTH.findall(text)]
    marker = "### sample " if "### sample " in text else "[derep] "
    parts = text.split(marker)[1:]
    return [BIRTH.findall(part) for part in parts]


def short_id(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:8]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else None


K = 10


def kmers(seq):
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def edit_end_free(a, b, band=12):
    """Banded unit-cost edit distance, free trailing overhang on either side.

    The free end matters for fixed-length reads: after an indel, one read runs
    on past the other's end. Hamming turns a single indel into dozens of
    mismatches (an AAC repeat unit scored 38); this scores it 3.
    """
    n, m = len(a), len(b)
    big = n + m
    prev = [j if j <= band else big for j in range(m + 1)]
    best = prev[m]
    for i in range(1, n + 1):
        cur = [big] * (m + 1)
        if i <= band:
            cur[0] = i
        for j in range(max(1, i - band), min(m, i + band) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        best = min(best, cur[m])
        prev = cur
    d = min(best, min(prev))
    # Outside the band the DP is only an upper bound; report it as "> band".
    return d if d <= band else band + 1


class EditIndex:
    """Nearest edit distance to a pool, via a k-mer prefilter."""

    def __init__(self, pool, top=5):
        self.top = top
        self.index = {}
        for s in pool:
            for k in kmers(s):
                self.index.setdefault(k, []).append(s)

    def nearest(self, seq):
        shared = {}
        for k in kmers(seq):
            for s in self.index.get(k, ()):
                shared[s] = shared.get(s, 0) + 1
        cands = sorted(shared, key=shared.get, reverse=True)[: self.top]
        return min((edit_end_free(seq, c) for c in cands if c != seq), default=None)


def nearest(seq, pool):
    """Smallest equal-length Hamming distance from seq to any sequence in pool."""
    ds = [h for p in pool if (h := hamming(seq, p)) is not None]
    return min(ds) if ds else None


def compare(a, b):
    """Churn between two arms' pooled tables."""
    only_a, only_b = set(a) - set(b), set(b) - set(a)
    moved = sum(abs(a.get(s, 0) - b.get(s, 0)) for s in set(a) | set(b)) // 2
    near = [nearest(s, set(a)) for s in only_b]
    idx = EditIndex(set(a)) if only_b else None
    near_edit = [idx.nearest(s) for s in only_b]
    return {
        "n_a": len(a), "n_b": len(b),
        "only_a": len(only_a), "only_b": len(only_b),
        "churn": len(only_a) + len(only_b),
        "reads_moved": moved,
        "reads_total": sum(a.values()),
        "only_a_max_abund": max((a[s] for s in only_a), default=0),
        "only_b_max_abund": max((b[s] for s in only_b), default=0),
        "only_b_h1": sum(1 for h in near if h == 1),
        "only_b_h_le2": sum(1 for h in near if h is not None and h <= 2),
        "only_b_e_le3": sum(1 for e in near_edit if e is not None and e <= 3),
    }


def main():
    argv, ids = sys.argv[1:], []
    while "--id" in argv:
        i = argv.index("--id")
        if i + 1 >= len(argv):
            sys.exit("--id needs a value")
        ids.append(argv[i + 1].lower())
        del argv[i:i + 2]
    if len(argv) < 2:
        sys.exit(__doc__)
    base_dir, arm_dirs = argv[0], argv[1:]
    # Every arm must come from the same command, or the comparison is between
    # modes, not member orders (three "modes" once all ran pooled, unnoticed).
    tags = {d: commands(d) for d in [base_dir, *arm_dirs]}
    if any(len(t) != 1 for t in tags.values()) or len(set().union(*tags.values())) != 1:
        for d, t in tags.items():
            print(f"  {d}: {sorted(t)}", file=sys.stderr)
        sys.exit("ERROR: arms were not all produced by one command; refusing to compare")
    print(f"command: {next(iter(tags[base_dir]))}")
    base = load(base_dir)
    arms = {os.path.basename(d.rstrip("/")): load(d) for d in arm_dirs}
    name = os.path.basename(base_dir.rstrip("/"))

    print(f"baseline {name}: {len(base)} ASVs, {sum(base.values())} reads\n")
    print(f"{'arm':12} {'ASVs':>5} {'churn':>6} {'-base':>6} {'+arm':>5} "
          f"{'reads moved':>18} {'max lost':>9} {'max new':>8} {'new H1':>7} {'new H<=2':>8} {'new E<=3':>8}")
    for arm, t in arms.items():
        c = compare(base, t)
        print(f"{arm:12} {c['n_b']:5} {c['churn']:6} {c['only_a']:6} {c['only_b']:5} "
              f"{c['reads_moved']:9} ({100 * c['reads_moved'] / c['reads_total']:.3f}%) "
              f"{c['only_a_max_abund']:9} {c['only_b_max_abund']:8} "
              f"{c['only_b_h1']:7} {c['only_b_h_le2']:8} {c['only_b_e_le3']:8}")

    # Which ASVs churn, by id, so the same flip across arms is visible.
    print("\nchurned ASVs (- lost from baseline, + new in arm; abundance; nearest to the other side")
    print("  as H = Hamming, equal length only, and E = edit distance with a free end):")
    base_idx = EditIndex(set(base))
    for arm, t in arms.items():
        lost, new = set(base) - set(t), set(t) - set(base)
        if not (lost or new):
            continue
        arm_idx = EditIndex(set(t))
        rows = [("-", s, base[s], nearest(s, set(t)), arm_idx.nearest(s)) for s in lost]
        rows += [("+", s, t[s], nearest(s, set(base)), base_idx.nearest(s)) for s in new]
        fmt_e = lambda e: "E>12" if e is not None and e > 12 else f"E{e}"
        desc = ", ".join(f"{sign}{short_id(s)}:{n} (H{h} {fmt_e(e)})" for sign, s, n, h, e in
                         sorted(rows, key=lambda r: -r[2]))
        print(f"  {arm:12} {desc}")

    command = next(iter(tags[base_dir]))
    b0 = births(base_dir, command)
    if command == "dada-pseudo":
        print("\nbirth order: skipped for dada-pseudo (round-1 lines are not attributable to samples)")
    elif b0 is not None:
        print(f"\nbirth order ({sum(map(len, b0))} births in {len(b0)} partition(s) in baseline):")
        for arm, d in zip(arms, arm_dirs):
            b = births(d, command)
            if b is None:
                continue
            if len(b) != len(b0):
                print(f"  {arm:12} partition count differs ({len(b)} vs {len(b0)}); not comparable")
                continue
            set_diff = sum(sorted(x) != sorted(y) for x, y in zip(b0, b))
            order_diff = sum(x != y for x, y in zip(b0, b))
            print(f"  {arm:12} partitions with a different birth set {set_diff}, "
                  f"with the same set in a different order {order_diff - set_diff}")

    if ids:
        show_ids(ids, [base_dir, *arm_dirs], [name, *arms])

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
