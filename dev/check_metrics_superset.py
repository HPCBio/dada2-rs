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
    dada2-rs dada ... --verbose --metrics-json m.json --metrics-attribution 2> v.txt
    dev/check_metrics_superset.py v.txt m.json

Exit status is 1 when a topic the prose reports has no JSON home, so this can
gate the follow-up PR that strips the prose.
"""

import json
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


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    verbose_path, json_path = sys.argv[1], sys.argv[2]

    lines = open(verbose_path, errors="replace").read().splitlines()
    doc = json.load(open(json_path))

    if not doc.get("runs"):
        sys.exit("metrics JSON has no runs[]")
    run = doc["runs"][0]

    print(f"schema_version : {doc.get('schema_version')}")
    print(f"measure_level  : {doc.get('measure_level')}")
    print(f"runs           : {len(doc['runs'])}")
    print(f"prose lines    : {sum(1 for l in lines if l.startswith('[dada'))}")
    print()

    topics = prose_topics(lines)
    if not topics:
        sys.exit("no [dada] topic headers found -- was --verbose passed?")

    carried, staying, missing = [], [], []
    for topic in topics:
        path = TOPICS[topic]
        if path is None:
            staying.append(topic)
            continue
        ok, val = dig(run, path)
        if ok and val is not None:
            carried.append((topic, path))
        else:
            missing.append((topic, path))

    print(f"CARRIED BY JSON ({len(carried)}) -- prose may be removed")
    for topic, path in carried:
        print(f"    {topic:36s} -> runs[].{path}")
    print()
    print(f"STAYING IN --verbose ({len(staying)}) -- run shape and warnings")
    for topic in staying:
        print(f"    {topic}")
    print()

    if missing:
        print(f"*** NOT CARRIED ({len(missing)}) -- removing these lines WOULD LOSE DATA")
        for topic, path in missing:
            print(f"    {topic:36s} -- expected runs[].{path}")
        print()
        print("Add these to the schema before stripping the prose, or record an")
        print("explicit decision that the quantity is not worth keeping.")
        return 1

    print("Every migrated topic has a JSON home. Safe to strip the prose.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
