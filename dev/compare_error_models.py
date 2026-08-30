#!/usr/bin/env python3
"""compare_error_models.py — is a learned error model still moving?

Given learn-errors JSONs in ascending order of whatever was varied (usually
`--nbases`), report how far each model sits from the next. The question a
convergence ladder answers is **"has err_out stopped moving?"**, not "did the ASV
count change" — an ASV count can sit still while the model underneath it drifts,
which is the trap this repo has hit repeatedly with `n_asv`.

Distances are computed on log10(rate), because error rates span many orders of
magnitude and a relative change at 1e-7 matters as much as one at 1e-2. Only the
off-diagonal (substitution) entries are compared; the diagonal is forced to a
constant by the fitter and carries no information.

Usage:
  compare_error_models.py LABEL=model.json [LABEL=model.json ...]
"""

import json
import math
import sys

# Rows of the 16-row matrix that are self-transitions (A2A, C2C, G2G, T2T).
DIAGONAL = (0, 5, 10, 15)


def load(path):
    j = json.load(open(path))
    rows = j["err_out"]
    off = [r for i, r in enumerate(rows) if i not in DIAGONAL]
    trans = sum(sum(r) for r in j["trans"]) if "trans" in j else None
    return {
        "rows": off,
        "nq": j.get("nq", len(rows[0])),
        "iters": j.get("iterations"),
        "converged": j.get("converged"),
        "trans": trans,
    }


def dist(a, b):
    """Max and RMS difference in log10(rate), aligned by QUALITY COLUMN.

    `nq` is `max_q + 1` and therefore varies with the data: a model fitted on
    fewer samples can top out at Q39 (nq=40) where another reaches Q40 (nq=41).
    Comparing must align column j of one matrix to column j of the other and
    truncate to the shorter.

    Flattening each matrix and zipping the flat lists -- the obvious thing, and
    what this script did first -- shifts every row by the nq difference, and the
    shift accumulates down the matrix. It manufactured apparent 100x model
    changes that fell *exactly* on nq boundaries and vanished once aligned. A
    comparison whose largest differences all sit on a metadata boundary is
    measuring the metadata.
    """
    floor = 1e-12
    nq = min(len(a["rows"][0]), len(b["rows"][0]))
    diffs = []
    for ra, rb in zip(a["rows"], b["rows"]):
        for x, y in zip(ra[:nq], rb[:nq]):
            diffs.append(
                abs(math.log10(max(x, floor)) - math.log10(max(y, floor)))
            )
    rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    return max(diffs), rms


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    models = []
    for a in args:
        label, _, path = a.partition("=")
        models.append((label, load(path)))

    print(f"{'model':>14s} {'nq':>4s} {'iters':>5s} {'conv':>5s} {'transitions':>14s}")
    for label, m in models:
        t = f"{m['trans']:,}" if m["trans"] is not None else "-"
        print(f"{label:>14s} {m['nq']:4d} {m['iters']!s:>5s} {str(m['converged']):>5s} {t:>14s}")

    print()
    print("Step-to-step movement in log10(rate), off-diagonal entries only.")
    print("A ladder has converged when these approach zero; if the last step is")
    print("as large as the first, the model is still moving and the top of the")
    print("ladder is not a converged model.")
    print()
    print(f"{'step':>28s} {'max Δlog10':>12s} {'rms Δlog10':>12s} {'max fold':>10s}")
    for (la, ma), (lb, mb) in zip(models, models[1:]):
        mx, rms = dist(ma, mb)
        print(f"{la + ' -> ' + lb:>28s} {mx:12.4f} {rms:12.4f} {10 ** mx:9.2f}x")

    # Distance of every rung from the top rung, which is the best available
    # estimate of "converged" on this dataset.
    top_label, top = models[-1]
    print()
    print(f"{'distance from ' + top_label:>28s} {'max Δlog10':>12s} {'rms Δlog10':>12s} {'max fold':>10s}")
    for label, m in models[:-1]:
        mx, rms = dist(m, top)
        print(f"{label:>28s} {mx:12.4f} {rms:12.4f} {10 ** mx:9.2f}x")


if __name__ == "__main__":
    main()
