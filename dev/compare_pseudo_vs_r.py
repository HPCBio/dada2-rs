#!/usr/bin/env python3
"""Which of our pseudo-pooling behaviours does R's `pool="pseudo"` resemble?

Issue #100. We believe R re-fits the error model between its two pseudo rounds
(dada.R:371-378) while we do not. That belief is a code reading. This is an
empirical test of it that does NOT require patching R:

  * if R re-estimates, its ASV table should sit CLOSER to our
    `--reestimate-err-between-rounds` arm than to our default;
  * if it does not, R should sit closer to our default.

The test is deliberately RELATIVE. Our output and R's differ for many unrelated
reasons (alignment kernel, tie-breaks, floating point), so absolute agreement is
not the signal and equality is not expected. What discriminates is the *direction*:
which arm is nearer, and by how much relative to the gap between our own two arms.

    d(R, default)  vs  d(R, reest)      compared against    d(default, reest)

If our two arms differ from each other by far less than either differs from R,
this test has no power on that dataset and says so rather than reporting a winner.

CONFOUND YOU MUST CONTROL: our default loess surface is `direct`; R DADA2's
`loessErrfun` uses R's default `interpolate`. That difference alone moves ~1 read
per sample on a 362-sample benchmark -- the same order as the effect being
measured here. Run our arms with `--loess-preset r-dada2` (or, for bit-parity,
`--errfun external --errfun-cmd "Rscript examples/external_errfun/loess_reference.R"`)
so the error-model fit is not a second difference. A sequence table records no
error-model provenance, so this script cannot verify it -- it prints the reminder
unconditionally and leaves the check to you.

Inputs are sequence tables, not per-sample dada output, because R's per-sample
round-2 output is never exposed -- one `dada()` call produces the whole pseudo run.
That also means the ASV-provenance measure (see
docs/findings/pseudo-pooling-priors-vs-error-model.md) cannot be computed for R;
table-level agreement is what is available.

Usage:
  compare_pseudo_vs_r.py --r-table seqtab_R.csv \
      --default seqtab.default.json --reest seqtab.reest.json [--min-abundance N]

  --r-table     R long CSV (`sequence,sample,count`), as written by
                dev/benchmark/run_dada2_pooled.R or dev/concordance/write_reference.R.
                Prefer the PRE-chimera table (`seqtab_R.csv`): chimera removal can
                mask or redistribute a denoising difference.
  --default     our seqtab JSON from a normal dada-pseudo run
  --reest       our seqtab JSON from --reestimate-err-between-rounds

Pure stdlib.
"""

import argparse
import csv
import gzip
import json
import sys


def load_r_csv(path):
    """{sequence: total count} from a long CSV."""
    opener = gzip.open if path.endswith(".gz") else open
    out = {}
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh):
            seq = row["sequence"].strip().upper()
            out[seq] = out.get(seq, 0) + int(row["count"])
    return out


def load_rs_json(path):
    """{sequence: total count} from a make-sequence-table / remove-bimera JSON."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        d = json.load(fh)
    if "sequences" not in d:
        for v in d.values():
            if isinstance(v, dict) and "sequences" in v:
                d = v
                break
    seqs = [s.upper() for s in d["sequences"]]
    counts = d.get("counts") or []
    totals = [0] * len(seqs)
    for row in counts:
        for j, c in enumerate(row):
            totals[j] += c
    return dict(zip(seqs, totals)), d


def agreement(a, b):
    """Set + count agreement between two {seq: total} maps."""
    sa, sb = set(a), set(b)
    shared = sa & sb
    union = sa | sb
    jaccard = len(shared) / len(union) if union else 1.0
    within2 = sum(
        1 for s in shared
        if a[s] and b[s] and 0.5 <= a[s] / b[s] <= 2.0
    )
    return {
        "jaccard": jaccard,
        "shared": len(shared),
        "only_a": len(sa - sb),
        "only_b": len(sb - sa),
        "counts_within_2x": within2 / len(shared) if shared else 1.0,
        # distance: 1 - jaccard, the quantity compared between arms
        "distance": 1.0 - jaccard,
    }


def filter_abund(m, n):
    return {s: c for s, c in m.items() if c >= n} if n else m


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--r-table", required=True)
    ap.add_argument("--default", required=True)
    ap.add_argument("--reest", required=True)
    ap.add_argument("--min-abundance", type=int, default=0,
                    help="Ignore ASVs with total count below this (advisory second view)")
    args = ap.parse_args()

    r = load_r_csv(args.r_table)
    dft, _ = load_rs_json(args.default)
    rst, _ = load_rs_json(args.reest)

    print("=" * 72)
    print("Does R's pool=\"pseudo\" resemble our default or our re-estimated arm?")
    print("=" * 72)
    print(f"\nASVs: R {len(r)}   default {len(dft)}   reest {len(rst)}")

    views = [("all ASVs", 0)]
    if args.min_abundance > 1:
        views.append((f"total count >= {args.min_abundance}", args.min_abundance))
    for label, n in views:
        R, D, S = (filter_abund(m, n) for m in (r, dft, rst))
        a_rd, a_rs, a_ds = agreement(R, D), agreement(R, S), agreement(D, S)
        print(f"\n--- {label} ---")
        for nm, a in (("R vs default", a_rd), ("R vs reest", a_rs),
                      ("default vs reest (our own gap)", a_ds)):
            print(f"  {nm:<32} jaccard {a['jaccard']:.4f}  "
                  f"shared {a['shared']:<6} +{a['only_b']:<5} -{a['only_a']:<5} "
                  f"counts within 2x {a['counts_within_2x']:.3f}")

        gap = abs(a_rd["distance"] - a_rs["distance"])
        own = a_ds["distance"]
        print()
        if own == 0:
            print("  NO POWER: our two arms produced identical tables here; this")
            print("  dataset cannot discriminate. Use data where the arms differ.")
        elif gap < 0.1 * own:
            print(f"  NO POWER: |d(R,default) - d(R,reest)| = {gap:.5f} is under 10% of")
            print(f"  our own arm gap ({own:.5f}). R is not meaningfully nearer either.")
        else:
            nearer = "default" if a_rd["distance"] < a_rs["distance"] else "reest"
            print(f"  R is closer to our {nearer.upper()} arm "
                  f"(difference {gap:.5f} vs own gap {own:.5f}).")
            if nearer == "reest":
                print("  -> CONSISTENT with R re-fitting the error model between rounds.")
            else:
                print("  -> AGAINST the re-fit hypothesis; R looks priors-only here.")
            print("  Relative test only: unrelated differences (alignment kernel,")
            print("  tie-breaks) displace both arms from R and are not controlled by it.")

    # A sequence table carries no error-model provenance (samples / sequences /
    # counts only), so this cannot be checked from the inputs -- state it every
    # time instead of pretending to verify it.
    print("\n  CHECK BEFORE TRUSTING THE VERDICT: were both of our arms run with an")
    print("  R-matched error model (--loess-preset r-dada2, or --errfun external with")
    print("  loess_reference.R)? Our default loess surface is 'direct' while R's")
    print("  loessErrfun uses 'interpolate'; that alone moves ~1 read/sample on a")
    print("  362-sample run -- the same order as the effect under test. If the arms")
    print("  used our default surface, both are displaced from R for an unrelated")
    print("  reason and the direction above is not trustworthy.")
    print("\n" + "=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
