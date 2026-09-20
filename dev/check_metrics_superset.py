#!/usr/bin/env python3
"""check_metrics_superset.py -- does --metrics-json carry every number --verbose prints?

Issue #162 migrates the machine-readable metrics out of `--verbose` prose and
into `--metrics-json`. Before any prose line is removed, the JSON has to be
shown to be a SUPERSET of it -- otherwise a quantity that
`docs/findings/data/*.txt` archives as primary data silently stops being
collected, and nobody notices until an analysis needs it.

This checks that claim mechanically. It is deliberately a claim about
CONTENT, not formatting: for each prose topic it names the JSON path that must
carry it, and reports the ones that nothing carries.

Usage:
    # one run
    dada2-rs dada ... --verbose --metrics-json m.json --metrics-attribution 2> v.txt
    dev/check_metrics_superset.py v.txt m.json

    # a whole sweep: phase_split.txt holds EVERY arm, so each is checked
    # against its own metrics/<arm>.json
    dev/check_metrics_superset.py <sweep-out-dir>

    # also compare the NUMBERS, not just which fields exist
    dev/check_metrics_superset.py --numeric <sweep-out-dir>

Checking a multi-arm phase_split.txt against a single arm's JSON reports false
gaps -- the minimizer arms' prose has no home in a k-mer arm's document, and
rightly so. Pass the sweep directory and it splits by `===== <arm>` instead.

Exit status is 1 when a topic the prose reports has no JSON home, so this can
gate the follow-up PR that strips the prose.
"""

import json
import os
import re
import sys

# Prose topic -> the JSON path that must carry it (checked against runs[0]).
# `None` means "deliberately not migrated": it stays in --verbose as run shape
# or a warning, per the tiering decided on #162.
TOPICS = {
    # --- migrated: these must exist in JSON before the prose can go ---
    "compare attribution": "compare.attribution",
    "compare split": "compare.split",
    "map parallel efficiency": "compare.map_parallel_efficiency",
    "phase times": "phases",
    "shuffle phases": "shuffle",
    "shuffle scan split": "shuffle.comps_build",
    "shuffle scan time": "shuffle.build",
    "shuffle redundancy": "shuffle.moves",
    "bud redundancy": "bud",
    "p-update churn": "p_update",
    "reconcile incremental (#136)": "reconcile",
    "reconcile rescan-necessity (#136)": "reconcile",
    "move pruning (#132)": "move_pruning",
    "#87 carry (#139)": "carry_87",
    "#87 projection (#139)": "carry_87",
    # --- deliberately staying in --verbose (run shape / warnings) ---
    "cpu allocation": None,
    "alignment backend": None,
    "tuning gates": None,
    "resident Raw footprint": "footprint",
    "minimizer index": "index",
    "minimizer sketch": None,
    "minimizer pooled diversity": None,
    "kmer8 fill": None,
    "kmer8 pooled diversity": None,
    "warning": None,
    "wrote": None,
}


def dig(obj, path):
    """Follow a dotted path; return (found, value)."""
    cur = obj
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return False, None
        cur = cur[part]
    return True, cur


def prose_topics(lines):
    """Topic headers the run actually printed, in order of first appearance."""
    seen = []
    for line in lines:
        m = re.match(r"^\[dada[^\]]*\]\s{1,3}(\S.*)$", line)
        if not m:
            continue
        body = m.group(1)
        # Indented continuation rows belong to the header above them.
        if line.startswith("[dada]   ") or line.startswith("[dada]     "):
            continue
        for topic in TOPICS:
            if body.startswith(topic) and topic not in seen:
                seen.append(topic)
                break
    return seen


# Prose row -> JSON field, inside the `compare attribution` block. The prose
# prints 2 decimals, so a summed comparison carries +/-0.005s of rounding per
# block; the tolerance below is derived from that, not guessed.
ATTRIB_ROWS = {
    "map": "map",
    "reduction": "reduction",
    "store": "store",
    "free": "free",
    "setup": "setup",
    "unattributed": "unattributed",
}


def prose_attrib_sums(lines):
    """Sum each compare-attribution row across every block in this arm.

    Order-independent on purpose: a concurrent per-sample run emits prose in
    completion order while the JSON is sorted by sample, so pairing block N
    with runs[N] would compare different samples and call the mismatch a bug.
    """
    sums, blocks = {k: 0.0 for k in ATTRIB_ROWS}, 0
    inside = False
    for line in lines:
        if "compare attribution (of" in line:
            inside, blocks = True, blocks + 1
            continue
        if inside:
            m = re.match(r"^\[dada\]\s+(\w+)\s+([0-9.]+)s", line)
            if m and m.group(1) in ATTRIB_ROWS:
                sums[m.group(1)] += float(m.group(2))
                continue
            # The block ends at the first line that is not one of its rows.
            if not re.match(r"^\[dada\]\s{3,}", line):
                inside = False
    return sums, blocks


def numeric_check(lines, doc, label):
    """Do the prose and the JSON agree on the numbers, not just the fields?"""
    prose, blocks = prose_attrib_sums(lines)
    if blocks == 0:
        print("    (no compare-attribution blocks in this arm's prose)")
        return 0
    runs = doc["runs"]
    jsn = {k: 0.0 for k in ATTRIB_ROWS}
    n = 0
    for r in runs:
        a = r.get("compare", {}).get("attribution")
        if not a:
            continue
        n += 1
        for prose_key, json_key in ATTRIB_ROWS.items():
            jsn[prose_key] += float(a.get(json_key) or 0.0)

    print(f"    prose blocks {blocks}  vs  json runs {n}")

    # Printed for EVERY arm, not only mismatching ones: an earlier version
    # only showed these on a mismatch, which made a shortfall look specific to
    # the one arm that happened to trip the block check.
    #
    # NOTE these do not all have the same expected count. `compare attribution`,
    # `phase times` and `compare split` are emitted once per run_dada under
    # --verbose. `resident Raw footprint` sits in the FRESH-BUILD branch of
    # dada_uniques_cached, and the ASV summary comes from denoise_and_serialize,
    # which not every path uses. Read a shortfall in those two as a question,
    # not as proof of loss.
    if True:
        counts = {
            "compare attribution": sum(1 for l in lines if "compare attribution (of" in l),
            "phase times": sum(1 for l in lines if "] phase times" in l),
            "compare split": sum(1 for l in lines if "compare split (of" in l),
            "resident Raw footprint": sum(
                1 for l in lines if "] resident Raw footprint" in l
            ),
            "ASV summary": sum(1 for l in lines if " ASV(s) from " in l),
        }
        per_run = ("compare attribution", "phase times", "compare split")
        print(f"    per-topic header counts (json runs = {n}):")
        for k, v in counts.items():
            if v == n:
                flag = ""
            elif k in per_run:
                flag = "   <-- SHORT (one per run_dada, so this is a real gap)"
            else:
                flag = "   <-- fewer (not one-per-run; see note in source)"
            print(f"        {k:24s} {v:4d}{flag}")

    if blocks != n:
        print("    *** block count differs -- sums are not comparable")
        return 1

    # 2-decimal prose rounding, worst case, summed over the blocks.
    tol = 0.005 * blocks
    bad = 0
    print(f"    {'row':14s} {'prose':>12s} {'json':>12s} {'diff':>10s}   (tol {tol:.3f}s)")
    for k in ATTRIB_ROWS:
        d = abs(prose[k] - jsn[k])
        flag = "" if d <= tol else "  <-- MISMATCH"
        if flag:
            bad += 1
        print(f"    {k:14s} {prose[k]:12.3f} {jsn[k]:12.3f} {d:10.3f}{flag}")
    if bad:
        print("    *** the prose and the JSON disagree beyond prose rounding.")
        print("        They read the same accumulators, so a gap here is a real bug.")
    return bad


def split_arms(lines):
    """Split a sweep phase_split.txt into {arm: lines}. One section if unmarked."""
    arms, cur, name = {}, [], None
    for line in lines:
        m = re.match(r"^=====\s+(\S+)", line)
        if m:
            if name is not None:
                arms[name] = cur
            name, cur = m.group(1), []
            continue
        cur.append(line)
    if name is None:
        return None
    arms[name] = cur
    return arms


def check_one(lines, doc, label):
    """Returns the number of topics with no JSON home."""
    if not doc.get("runs"):
        print(f"{label}: metrics JSON has no runs[]")
        return 1
    runs = doc["runs"]

    print(f"schema_version : {doc.get('schema_version')}")
    print(f"measure_level  : {doc.get('measure_level')}")
    print(f"runs           : {len(doc['runs'])}")
    print(f"prose lines    : {sum(1 for l in lines if l.startswith('[dada'))}")
    print()

    topics = prose_topics(lines)
    if not topics:
        print(f"{label}: no [dada] topic headers found -- was --verbose passed?")
        return 1

    # Checked across EVERY run, not just the first. A per-sample sweep emits one
    # run per sample, and a field present in run 0 but absent in run 17 is a
    # partial migration -- which reads as success if you only look at the head.
    carried, staying, missing, partial = [], [], [], []
    for topic in topics:
        path = TOPICS[topic]
        if path is None:
            staying.append(topic)
            continue
        present = sum(1 for r in runs if dig(r, path)[1] is not None)
        if present == len(runs):
            carried.append((topic, path))
        elif present == 0:
            missing.append((topic, path))
        else:
            partial.append((topic, path, present, len(runs)))

    print(f"CARRIED BY JSON ({len(carried)}) -- prose may be removed")
    for topic, path in carried:
        print(f"    {topic:36s} -> runs[].{path}")
    print()
    print(f"STAYING IN --verbose ({len(staying)}) -- run shape and warnings")
    for topic in staying:
        print(f"    {topic}")
    print()

    if partial:
        print(f"*** PARTIAL ({len(partial)}) -- present in SOME runs only")
        for topic, path, n, tot in partial:
            print(f"    {topic:36s} -> runs[].{path}  ({n}/{tot} runs)")
        print()

    if missing:
        print(f"*** NOT CARRIED ({len(missing)}) -- removing these lines WOULD LOSE DATA")
        for topic, path in missing:
            print(f"    {topic:36s} -- expected runs[].{path}")
        print()
        print("Add these to the schema before stripping the prose, or record an")
        print("explicit decision that the quantity is not worth keeping.")
        return len(missing) + len(partial)

    if partial:
        return len(partial)

    print(f"Every migrated topic has a JSON home in all {len(runs)} run(s).")
    print("Safe to strip the prose.")
    return 0


def main():
    args = sys.argv[1:]
    numeric = "--numeric" in args
    args = [a for a in args if a != "--numeric"]

    # Sweep-directory mode: check each arm against its own document.
    if len(args) == 1:
        root = args[0]
        split_path = os.path.join(root, "phase_split.txt")
        if not os.path.isfile(split_path):
            sys.exit(f"{split_path}: not found (expected a sweep output directory)")
        arms = split_arms(open(split_path, errors="replace").read().splitlines())
        if arms is None:
            sys.exit(f"{split_path} has no `===== <arm>` markers; pass <verbose> <json>")
        bad, checked = 0, 0
        for arm, lines in arms.items():
            mpath = os.path.join(root, "metrics", f"{arm}.json")
            print("=" * 68)
            print(f"ARM {arm}")
            print("=" * 68)
            if not os.path.isfile(mpath):
                print(f"    no metrics/{arm}.json -- arm skipped or run before #162\n")
                continue
            checked += 1
            doc = json.load(open(mpath))
            bad += check_one(lines, doc, arm)
            if numeric:
                print("NUMERIC CROSS-CHECK (compare attribution, summed over blocks)")
                bad += numeric_check(lines, doc, arm)
            print()
        if checked == 0:
            sys.exit("no arm had a metrics JSON; was the sweep run on this branch?")
        print(f"{checked} arm(s) checked, {bad} unmapped topic(s)")
        return 1 if bad else 0

    if len(args) != 2:
        sys.exit(__doc__)

    lines = open(args[0], errors="replace").read().splitlines()
    arms = split_arms(lines)
    if arms is not None:
        print(
            f"NOTE: {args[0]} holds {len(arms)} arms "
            f"({', '.join(list(arms)[:4])}{' ...' if len(arms) > 4 else ''}).\n"
            "      Checking ALL of their prose against ONE arm's JSON reports false\n"
            "      gaps -- a k-mer arm has no minimizer index, and should not.\n"
            "      Pass the sweep directory instead to check each arm on its own.\n"
        )
    doc = json.load(open(args[1]))
    bad = check_one(lines, doc, args[1])
    if numeric:
        print("NUMERIC CROSS-CHECK (compare attribution, summed over blocks)")
        bad += numeric_check(lines, doc, args[1])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
