#!/usr/bin/env python3
"""List the rename pairs between member-order arms, per sample (issue #246).

A rename pair is one organism reported under two names: in one sample, an ASV
that only the baseline arm reports and an ASV that only the other arm reports,
within --max-edit edits of each other. Arms come from
dev/run_member_order.sh with MODE=per-sample (`dada`), so each pair can be
re-run on its own sample by trace_pairs.sh. Pooled and pseudo arms are refused:
their clusters span samples, and a one-sample trace would not reproduce them.

Pairing is one-to-one within a sample, nearest first. An ASV that churns with
no partner within --max-edit (absorbed into an ASV both arms keep, or a split)
is counted in the summary on stderr but not listed.

A pair that flips in several arms is listed once, with every such arm; the
first is the one trace_pairs.sh runs.

Usage:
  find_pairs.py <baseline-dir> <arm-dir> [<arm-dir> ...] [--max-edit N] [-o pairs.tsv]
Arm names are directory basenames as run_member_order.sh writes them
(shuffle_2 for shuffle:2). Without -o the table goes to stdout.
"""
import argparse
import os
import sys

# Reuse the floor tool's loaders and distances, so a pair here is a pair there.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from compare_member_order import commands, edit_end_free, hamming, load_per_sample, short_id  # noqa: E402

COLUMNS = ["sample", "base_id", "arm_id", "edit", "hamming", "base_reads", "arm_reads",
           "arms", "base_seq", "arm_seq"]


def pairs_in_sample(base, arm, max_edit):
    """One-to-one (lost, new, edit) pairs within max_edit, plus unpaired counts."""
    lost = [s for s in base if s not in arm]
    new = [s for s in arm if s not in base]
    cands = sorted(
        ((e, -min(base[a], arm[b]), a, b) for a in lost for b in new
         if (e := edit_end_free(a, b)) <= max_edit),
        key=lambda c: c[:2])
    used_a, used_b, out = set(), set(), []
    for e, _, a, b in cands:
        if a in used_a or b in used_b:
            continue
        used_a.add(a)
        used_b.add(b)
        out.append((a, b, e))
    return out, len(lost) - len(used_a), len(new) - len(used_b)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("baseline")
    ap.add_argument("arms", nargs="+")
    ap.add_argument("--max-edit", type=int, default=3)
    ap.add_argument("-o", "--output")
    a = ap.parse_args()

    dirs = [a.baseline, *a.arms]
    for d in dirs:
        tags = commands(d)
        if tags != {"dada"}:
            sys.exit(f"ERROR: {d} holds {sorted(tags) or 'no dada JSON'}; "
                     "need per-sample `dada` output (run_member_order.sh MODE=per-sample)")
    base = load_per_sample(a.baseline)
    pairs = {}  # (sample, base_seq, arm_seq) -> row
    print(f"{'arm':12} {'pairs':>6} {'unpaired -base':>15} {'unpaired +arm':>14}", file=sys.stderr)
    for d in a.arms:
        name = os.path.basename(d.rstrip("/"))
        arm = load_per_sample(d)
        if set(arm) != set(base):
            sys.exit(f"ERROR: {d} and {a.baseline} hold different samples")
        n = un_a = un_b = 0
        for smp in sorted(base):
            found, ua, ub = pairs_in_sample(base[smp], arm[smp], a.max_edit)
            un_a += ua
            un_b += ub
            for s, t, e in found:
                n += 1
                key = (smp, s, t)
                if key in pairs:
                    pairs[key]["arms"].append(name)
                    continue
                h = hamming(s, t)
                pairs[key] = {
                    "sample": smp, "base_id": short_id(s), "arm_id": short_id(t),
                    "edit": e, "hamming": "" if h is None else h,
                    "base_reads": base[smp][s], "arm_reads": arm[smp][t],
                    "arms": [name], "base_seq": s, "arm_seq": t,
                }
        print(f"{name:12} {n:6} {un_a:15} {un_b:14}", file=sys.stderr)

    rows = sorted(pairs.values(), key=lambda r: (-max(r["base_reads"], r["arm_reads"]), r["sample"]))
    print(f"{len(rows)} distinct pair(s) within edit {a.max_edit}", file=sys.stderr)
    out = open(a.output, "w") if a.output else sys.stdout
    print("\t".join(COLUMNS), file=out)
    for r in rows:
        r = dict(r, arms=",".join(r["arms"]))
        print("\t".join(str(r[c]) for c in COLUMNS), file=out)


if __name__ == "__main__":
    main()
