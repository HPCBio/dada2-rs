#!/usr/bin/env python3
"""Compare ASV sets between dada2-rs runs (A/B or baseline-vs-many).

Generalizes the ASV-tracing comparison used in the issue #15 k-mer-size
investigation: load each run's ASVs into a {sequence: abundance} map, diff each
run against a baseline, stratify the churned sets by abundance, and (for
equal-length ASVs) report the nearest-neighbor Hamming distance of each churned
ASV to the baseline set — the signal that distinguishes benign cluster
*fragmentation* (Hamming-1/2 from an existing ASV) from genuinely novel sequences.

It is schema-tolerant across the artifacts dada2-rs emits, auto-detected per input:
  * dada / dada-pooled : a DIRECTORY of per-sample `*.json` (asvs[] pooled), or a
                         single such file. Carries birth_type.
  * make-sequence-table / remove-bimera-denovo: a seqtab `.json`
                         (sequences[] + per-sample counts[][]; the latter is the
                         chimera-filtered final table).
  * merge-pairs        : a merged `.json` (samples[].merged[], accept-filtered).

Usage:
  compare_asvs.py --baseline LABEL=PATH --compare LABEL=PATH [--compare LABEL=PATH ...]

  PATH is a file (seqtab/merged JSON) or a directory (per-sample dada JSONs).
  LABEL is any short name (e.g. k5, k7, runA).

Examples:
  # dada-level, baseline k5 vs k6/k7/k8 (directories of per-sample dada JSONs)
  compare_asvs.py --baseline k5=dada_k5 --compare k6=dada_k6 --compare k7=dada_k7

  # final chimera-filtered tables (seqtab JSONs)
  compare_asvs.py --baseline k5=seqtab.nonchim.k5.json --compare k8=seqtab.nonchim.k8.json

  # merged amplicons, drop low-abundance noise, machine-readable out
  compare_asvs.py --baseline a=merged_a.json --compare b=merged_b.json \
      --min-abundance 10 --json report.json

Pure stdlib; no dependencies.
"""

import argparse
import json
import os
import statistics
import sys
from glob import glob


# ---------------------------------------------------------------------------
# Loading: each loader returns {sequence: {"abundance": int, "birth": str|None}}
# pooled across samples (abundances summed; first-seen birth_type kept).
# ---------------------------------------------------------------------------

def _accumulate(acc, seq, abundance, birth=None):
    e = acc.get(seq)
    if e is None:
        acc[seq] = {"abundance": int(abundance), "birth": birth}
    else:
        e["abundance"] += int(abundance)
        if e["birth"] is None and birth is not None:
            e["birth"] = birth


def _load_one_json(path, acc):
    """Add one JSON file's ASVs to acc. Returns the detected kind."""
    with open(path) as fh:
        d = json.load(fh)
    tag = d.get("dada2_rs_command") if isinstance(d, dict) else None

    if tag in ("dada", "dada-pooled"):
        for a in d.get("asvs", []):
            s = a.get("sequence")
            if s:
                _accumulate(acc, s, a.get("abundance", 0), a.get("birth_type"))
        return "dada"

    if tag in ("make-sequence-table", "remove-bimera-denovo"):
        # Both use the seqtab schema (sequences[] + samples[] x counts[][]);
        # remove-bimera-denovo is the chimera-filtered final table.
        seqs = d.get("sequences", [])
        counts = d.get("counts", [])  # samples x sequences
        totals = [0] * len(seqs)
        for row in counts:
            for j, c in enumerate(row):
                if j < len(totals):
                    totals[j] += c
        for s, t in zip(seqs, totals):
            if s:
                _accumulate(acc, s, t)
        return "seqtab"

    if tag == "merge-pairs":
        for samp in d.get("samples", []):
            for m in samp.get("merged", []):
                s = m.get("sequence")
                if s and m.get("accept", False):
                    _accumulate(acc, s, m.get("abundance", 0))
        return "merged"

    raise ValueError(
        f"{path}: unrecognized dada2_rs_command tag {tag!r} "
        f"(expected dada, dada-pooled, make-sequence-table, remove-bimera-denovo, or merge-pairs)"
    )


def load_run(path):
    """Load a run from a file or a directory of per-sample dada JSONs.
    Returns (asv_map, kind)."""
    acc = {}
    if os.path.isdir(path):
        files = sorted(glob(os.path.join(path, "*.json"))
                       + glob(os.path.join(path, "*.json.gz")))
        if not files:
            raise ValueError(f"{path}: directory contains no *.json files")
        kinds = set()
        for f in files:
            if f.endswith(".gz"):
                raise ValueError(f"{f}: gzipped input not supported; gunzip first")
            kinds.add(_load_one_json(f, acc))
        kind = kinds.pop() if len(kinds) == 1 else "mixed"
    else:
        kind = _load_one_json(path, acc)
    return acc, kind


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def filter_min_abund(asv_map, min_abundance):
    if min_abundance <= 1:
        return asv_map
    return {s: e for s, e in asv_map.items() if e["abundance"] >= min_abundance}


def abundance_summary(asv_map, seqs):
    """Stratify a subset of sequences by abundance."""
    ab = sorted((asv_map[s]["abundance"] for s in seqs), reverse=True)
    if not ab:
        return {"n": 0, "lt10": 0, "ge10": 0, "max": 0, "median": 0}
    return {
        "n": len(ab),
        "singletons": sum(1 for a in ab if a == 1),
        "lt10": sum(1 for a in ab if a < 10),
        "ge10": sum(1 for a in ab if a >= 10),
        "max": max(ab),
        "median": int(statistics.median(ab)),
    }


def hamming(a, b):
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))


def nearest_in(seq, pool_list):
    """Min Hamming distance from seq to any equal-length sequence in pool_list,
    plus that neighbor. Returns (dist, neighbor) or (None, None) if no equal-length."""
    best_d, best = None, None
    for t in pool_list:
        h = hamming(seq, t)
        if h is not None and (best_d is None or h < best_d):
            best_d, best = h, t
            if h == 0:
                break
    return best_d, best


def compare(baseline, base_map, other_label, other_map, do_hamming, n_examples):
    bset, oset = set(base_map), set(other_map)
    only_base = bset - oset
    only_other = oset - bset
    shared = bset & oset
    result = {
        "label": other_label,
        "baseline_n": len(bset),
        "other_n": len(oset),
        "shared": len(shared),
        f"only_{baseline}": len(only_base),
        f"only_{other_label}": len(only_other),
        "churn": len(only_base) + len(only_other),
        "abundance": {
            f"only_{baseline}": abundance_summary(base_map, only_base),
            f"only_{other_label}": abundance_summary(other_map, only_other),
        },
    }
    if do_hamming:
        # For the churned-in-other set (the "new at this run" ASVs), find nearest
        # baseline neighbor — the fragmentation signal.
        base_list = list(bset)
        examples = []
        for s in sorted(only_other, key=lambda s: -other_map[s]["abundance"])[:n_examples]:
            d, nb = nearest_in(s, base_list)
            examples.append({
                "abundance": other_map[s]["abundance"],
                "birth": other_map[s]["birth"],
                "min_hamming_to_baseline": d,
                "neighbor_abundance": base_map[nb]["abundance"] if nb else None,
            })
        # Distribution of nearest-neighbor distances over ALL churned-in-other.
        dists = []
        for s in only_other:
            d, _ = nearest_in(s, base_list)
            if d is not None:
                dists.append(d)
        result["hamming"] = {
            "n_with_equal_length_neighbor": len(dists),
            "hamming1": sum(1 for d in dists if d == 1),
            "le2": sum(1 for d in dists if d <= 2),
            "examples": examples,
        }
    return result


# ---------------------------------------------------------------------------
# CLI + reporting
# ---------------------------------------------------------------------------

def parse_label_path(s):
    if "=" not in s:
        raise argparse.ArgumentTypeError(f"expected LABEL=PATH, got {s!r}")
    label, path = s.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError(f"expected LABEL=PATH, got {s!r}")
    return label, path


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Compare ASV sets between dada2-rs runs (baseline vs one or more).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--baseline", required=True, type=parse_label_path,
                   metavar="LABEL=PATH", help="Baseline run (file or dada dir).")
    p.add_argument("--compare", required=True, action="append", type=parse_label_path,
                   metavar="LABEL=PATH", help="Run(s) to compare vs baseline. Repeatable.")
    p.add_argument("--min-abundance", type=int, default=1,
                   help="Drop ASVs below this pooled abundance before diffing (default 1).")
    p.add_argument("--no-hamming", action="store_true",
                   help="Skip nearest-neighbor Hamming analysis (faster; needed if lengths vary widely).")
    p.add_argument("--examples", type=int, default=8,
                   help="Number of top-abundance churned ASVs to detail per comparison (default 8).")
    p.add_argument("--json", metavar="OUT", help="Also write the full report as JSON.")
    args = p.parse_args(argv)

    blabel, bpath = args.baseline
    base_map, bkind = load_run(bpath)
    base_map = filter_min_abund(base_map, args.min_abundance)

    report = {
        "baseline": {"label": blabel, "path": bpath, "kind": bkind, "n_asv": len(base_map)},
        "min_abundance": args.min_abundance,
        "comparisons": [],
    }

    for olabel, opath in args.compare:
        omap, okind = load_run(opath)
        omap = filter_min_abund(omap, args.min_abundance)
        if okind != bkind:
            print(f"WARNING: {olabel} kind={okind} differs from baseline kind={bkind}; "
                  f"comparing anyway by sequence.", file=sys.stderr)
        report["comparisons"].append(
            compare(blabel, base_map, olabel, omap, not args.no_hamming, args.examples)
        )

    _print_report(report, blabel)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"\nWrote {args.json}")
    return 0


def _print_report(report, blabel):
    b = report["baseline"]
    print(f"\nBaseline {b['label']} ({b['kind']}): {b['n_asv']} ASVs"
          + (f"  [min_abundance={report['min_abundance']}]" if report["min_abundance"] > 1 else ""))
    print("=" * 64)
    for c in report["comparisons"]:
        ol = c["label"]
        print(f"\n{ol} ({c['other_n']} ASVs) vs {blabel} ({c['baseline_n']}):")
        print(f"  shared={c['shared']}  only_{blabel}={c[f'only_{blabel}']}  "
              f"only_{ol}={c[f'only_{ol}']}  churn={c['churn']}")
        ab_b = c["abundance"][f"only_{blabel}"]
        ab_o = c["abundance"][f"only_{ol}"]
        print(f"  only_{blabel} abundance: n={ab_b['n']} <10={ab_b['lt10']} "
              f">=10={ab_b['ge10']} max={ab_b['max']} median={ab_b['median']}")
        print(f"  only_{ol} abundance: n={ab_o['n']} <10={ab_o['lt10']} "
              f">=10={ab_o['ge10']} max={ab_o['max']} median={ab_o['median']}")
        if "hamming" in c:
            h = c["hamming"]
            print(f"  nearest-baseline Hamming (only_{ol}): "
                  f"Hamming-1={h['hamming1']} <=2={h['le2']} "
                  f"of {h['n_with_equal_length_neighbor']} equal-length")
            if h["examples"]:
                print(f"  top only_{ol} ASVs (fragmentation check):")
                for e in h["examples"]:
                    print(f"    abund={e['abundance']:>6}  birth={str(e['birth']):<10}  "
                          f"minH={e['min_hamming_to_baseline']}  "
                          f"neighbor_abund={e['neighbor_abundance']}")


if __name__ == "__main__":
    sys.exit(main())
