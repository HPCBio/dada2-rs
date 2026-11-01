#!/usr/bin/env python3
"""Share of denoising time a single-pass partition could remove (#323).

Reads --metrics-json files from `dada`, `dada-pooled` or `dada-pseudo` and sums
the per-run phase timers. Single pass keeps every comparison and drops the
shuffle, the bud search and the p-value refresh, so the sum of those three is
the most it can save. Shares are of the summed phases: `loop_total` excludes
the initial compare, so it is not used as the denominator.

A share is only meaningful with the thread count it was measured at: the
removable phases are serial, so their share grows with threads.

Usage: dev/phase_share.py <metrics.json> [...]
"""
import json
import sys

PHASES = ["compare", "shuffle", "bud", "p_update", "index_add"]
REMOVABLE = ["shuffle", "bud", "p_update"]


def runs(o):
    """Yield every object carrying a phases block (one per run_dada call)."""
    if isinstance(o, dict):
        if isinstance(o.get("phases"), dict) and "shuffle" in o["phases"]:
            yield o
        for v in o.values():
            yield from runs(v)
    elif isinstance(o, list):
        for v in o:
            yield from runs(v)


def main(paths):
    if not paths:
        sys.exit(__doc__)
    for path in paths:
        with open(path) as fh:
            rows = list(runs(json.load(fh)))
        if not rows:
            print(f"{path}: no phases block found")
            continue
        t = {k: sum(r["phases"].get(k) or 0.0 for r in rows) for k in PHASES}
        total = sum(t.values())
        rem = sum(t[k] for k in REMOVABLE)
        threads = rows[0].get("run", {}).get("threads", "?")

        def pct(x):
            return f"{x / total:.1%}" if total else "n/a"

        print(
            f"{path}: runs={len(rows)} threads={threads} total={total:.2f}s "
            f"compare={t['compare']:.2f}s ({pct(t['compare'])}) "
            f"shuffle={t['shuffle']:.2f}s bud={t['bud']:.2f}s "
            f"p_update={t['p_update']:.2f}s -> removable {rem:.2f}s ({pct(rem)})"
        )


if __name__ == "__main__":
    main(sys.argv[1:])
