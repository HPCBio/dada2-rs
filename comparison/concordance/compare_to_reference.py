#!/usr/bin/env python3
"""compare_to_reference.py — concordance guardrail: dada2-rs vs an R DADA2 reference.

CI sanity check (NOT an exact-match test). Compares a dada2-rs post-chimera
sequence table against a *static* R DADA2 reference CSV (checked into the repo,
produced once by `write_reference.R`; CI never runs R) on two axes:

  1. ASV-set agreement — recall (fraction of R's ASVs that dada2-rs also found),
     precision (fraction of dada2-rs ASVs that R also found), and Jaccard.
  2. Per-ASV count agreement — Pearson correlation of log10 total counts over the
     shared ASVs, plus the fraction of shared ASVs whose counts agree within 2×.

These will never be exactly equal (different alignment kernels, tie-breaks,
floating point), so the gate is threshold-based and meant to catch *regressions /
divergence*, not enforce identity. Defaults are warn-only (always exit 0); pass
`--gate` to fail when a metric falls below its threshold.

Inputs:
  --rs PATH         dada2-rs seqtab JSON (make-sequence-table or, preferred, the
                    chimera-filtered remove-bimera-denovo output). `.json[.gz]`.
  --reference PATH  R reference CSV, long format: columns `sequence,sample,count`
                    (header required; zero-count rows may be omitted). Sequences
                    are upper-cased and compared exactly.

Comparison is on TOTAL count per ASV (summed across samples) by default, which is
robust to per-sample assignment noise; `--per-sample` additionally reports the
per-(ASV,sample) count correlation as advisory.

Usage:
  compare_to_reference.py --rs seqtab.nochim.json --reference reference/illumina.csv
  compare_to_reference.py --rs out.json --reference ref.csv --gate \
      --min-recall 0.85 --min-precision 0.80 --min-count-corr 0.95 \
      --min-abundance 2 --summary "$GITHUB_STEP_SUMMARY" --json report.json

Pure stdlib; no dependencies.
"""

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, newline="")


def load_rs_seqtab(path):
    """dada2-rs seqtab JSON -> {SEQUENCE: total_count}. Accepts make-sequence-table
    or remove-bimera-denovo (same schema: sequences[] + counts[sample][seq])."""
    with _open(path) as fh:
        d = json.load(fh)
    tag = d.get("dada2_rs_command")
    if tag not in ("make-sequence-table", "remove-bimera-denovo"):
        sys.exit(
            f"{path}: dada2_rs_command={tag!r}; expected a seqtab "
            "(make-sequence-table or remove-bimera-denovo)"
        )
    seqs = d.get("sequences", [])
    totals = [0] * len(seqs)
    for row in d.get("counts", []):  # counts[sample][seq]
        for j, c in enumerate(row):
            if j < len(totals):
                totals[j] += c
    out = {}
    for s, t in zip(seqs, totals):
        if s:
            out[s.upper()] = out.get(s.upper(), 0) + t
    return out


def load_reference_csv(path):
    """R reference CSV (long: sequence,sample,count) -> {SEQUENCE: total_count}."""
    out = {}
    with _open(path) as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        for need in ("sequence", "count"):
            if need not in cols:
                sys.exit(
                    f"{path}: reference CSV must have a '{need}' column "
                    f"(found: {reader.fieldnames})"
                )
        scol, ccol = cols["sequence"], cols["count"]
        for r in reader:
            seq = (r.get(scol) or "").strip().upper()
            if not seq:
                continue
            try:
                c = float(r.get(ccol) or 0)
            except ValueError:
                continue
            if c > 0:
                out[seq] = out.get(seq, 0) + c
    if not out:
        sys.exit(f"{path}: reference CSV produced no ASVs")
    return out


def filter_min_abund(m, min_abundance):
    if min_abundance <= 1:
        return m
    return {s: c for s, c in m.items() if c >= min_abundance}


def pearson(xs, ys):
    """Pearson correlation; None if undefined (n<2 or zero variance)."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy)


def compare(rs, ref):
    rs_set, ref_set = set(rs), set(ref)
    shared = rs_set & ref_set
    only_rs = rs_set - ref_set
    only_ref = ref_set - rs_set
    union = rs_set | ref_set

    recall = len(shared) / len(ref_set) if ref_set else 0.0
    precision = len(shared) / len(rs_set) if rs_set else 0.0
    jaccard = len(shared) / len(union) if union else 0.0

    # Count agreement over shared ASVs, on log10 of totals.
    lx = [math.log10(ref[s]) for s in shared]
    ly = [math.log10(rs[s]) for s in shared]
    count_corr = pearson(lx, ly)
    within_2x = (
        sum(1 for s in shared if 0.5 <= (rs[s] / ref[s]) <= 2.0) / len(shared)
        if shared
        else 0.0
    )

    rs_reads, ref_reads = sum(rs.values()), sum(ref.values())
    return {
        "rs_asvs": len(rs_set),
        "ref_asvs": len(ref_set),
        "shared": len(shared),
        "only_rs": len(only_rs),
        "only_ref": len(only_ref),
        "recall": recall,
        "precision": precision,
        "jaccard": jaccard,
        "count_corr_log10": count_corr,
        "within_2x_frac": within_2x,
        "rs_total_reads": rs_reads,
        "ref_total_reads": ref_reads,
        "reads_ratio_rs_over_ref": (rs_reads / ref_reads) if ref_reads else None,
    }


def fmt(v):
    return "n/a" if v is None else f"{v:.3f}"


def render(c, thresholds, gate):
    lines = []
    lines.append(f"ASVs: dada2-rs {c['rs_asvs']}  vs  R {c['ref_asvs']}  "
                 f"(shared {c['shared']}, only-rs {c['only_rs']}, only-R {c['only_ref']})")
    lines.append(f"reads retained (rs/R): {fmt(c['reads_ratio_rs_over_ref'])}")
    checks = [
        ("recall (R ASVs found)", "recall", thresholds.get("recall")),
        ("precision (rs ASVs in R)", "precision", thresholds.get("precision")),
        ("count corr (log10)", "count_corr_log10", thresholds.get("count_corr")),
    ]
    ok = True
    for label, key, thr in checks:
        val = c[key]
        if thr is None:
            lines.append(f"  {label}: {fmt(val)}")
            continue
        passed = val is not None and val >= thr
        ok = ok and passed
        mark = "PASS" if passed else "FAIL"
        lines.append(f"  {label}: {fmt(val)}  (>= {thr:.2f})  [{mark}]")
    lines.append(f"  within-2x count agreement: {fmt(c['within_2x_frac'])}  (advisory)")
    lines.append(f"  jaccard: {fmt(c['jaccard'])}  (advisory)")
    status = "PASS" if ok else ("FAIL" if gate else "FAIL (advisory — not gating)")
    lines.append(f"==> {status}")
    return "\n".join(lines), ok


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rs", required=True, help="dada2-rs seqtab JSON (.json[.gz])")
    p.add_argument("--reference", required=True, help="R reference CSV (sequence,sample,count)")
    p.add_argument("--label", default="run", help="label for the report header")
    p.add_argument("--min-abundance", type=int, default=1,
                   help="drop ASVs with total count below this on BOTH sides before "
                        "comparing (filters singleton noise). Default 1 (no filter)")
    p.add_argument("--min-recall", type=float, default=None,
                   help="threshold: fraction of R's ASVs that dada2-rs must also find")
    p.add_argument("--min-precision", type=float, default=None,
                   help="threshold: fraction of dada2-rs ASVs that R also found")
    p.add_argument("--min-count-corr", type=float, default=None,
                   help="threshold: Pearson r of log10 total counts over shared ASVs")
    p.add_argument("--gate", action="store_true",
                   help="exit non-zero when a thresholded metric fails (default: "
                        "warn-only, always exit 0)")
    p.add_argument("--summary", metavar="FILE",
                   help="append a markdown summary to FILE (e.g. $GITHUB_STEP_SUMMARY)")
    p.add_argument("--json", metavar="OUT", help="write the full report as JSON")
    args = p.parse_args(argv)

    rs = filter_min_abund(load_rs_seqtab(args.rs), args.min_abundance)
    ref = filter_min_abund(load_reference_csv(args.reference), args.min_abundance)
    c = compare(rs, ref)

    thresholds = {
        "recall": args.min_recall,
        "precision": args.min_precision,
        "count_corr": args.min_count_corr,
    }
    body, ok = render(c, thresholds, args.gate)
    header = (f"Concordance vs R DADA2 — {args.label}"
              + (f"  [min_abundance={args.min_abundance}]" if args.min_abundance > 1 else ""))
    print(header)
    print(body)

    if args.summary:
        with open(args.summary, "a") as fh:
            fh.write(f"### {header}\n\n```\n{body}\n```\n\n")

    if args.json:
        Path(args.json).write_text(json.dumps(
            {"label": args.label, "min_abundance": args.min_abundance,
             "thresholds": thresholds, "gate": args.gate, "metrics": c}, indent=2))

    if args.gate and not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
