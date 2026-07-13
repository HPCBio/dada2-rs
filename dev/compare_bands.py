#!/usr/bin/env python3
"""Compare dada2-rs pooled runs across --band values (the BAND_SIZE A/B).

Companion to run_band_sweep.sh. Given the per-band pooled records
(`pooled_band<B>.json[.gz]`, written by `dada-pooled --pooled-record`) and the
per-band error models (`errors_band<B>.json`), it evaluates the two pass
criteria for "a tighter band is safe":

  1. ASV-set identity  — the set of ASV sequences is unchanged vs the baseline
     band (always required). This is the concordance-guardrail gate (project #35).
  2. Abundance impact — the fraction of reads that change assignment stays under
     --max-reads-moved-pct (default 0.01%). Partition churn is reported per
     unique with a MULTI-READ breakdown and center-aware labels (failed /
     absorbed / center), so a lost ASV ("center->failed") or an offset-variant
     reshuffle is called out explicitly, but a benign Hamming-1 member
     reassignment that leaves the ASV set intact does not fail the run.

The learned error matrix (`err_out`) max-abs delta vs baseline is reported too,
but as INFORMATIONAL context, not a gate.

Everything is keyed by SEQUENCE, not array index, so it is robust to any
reordering of the merged unique table across runs. Pure stdlib.

Usage:
  compare_bands.py --baseline 16=OUTDIR/pooled_band16.json.gz \
      --compare 8=OUTDIR/pooled_band8.json.gz \
      --compare 4=OUTDIR/pooled_band4.json.gz \
      [--errors 16=OUTDIR/errors_band16.json --errors 8=...] \
      [--json report.json]
"""

import argparse
import gzip
import json
import sys


def _load(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        d = json.load(fh)
    # Unwrap a {"tag":..., "data":{...}} envelope if present.
    if isinstance(d, dict) and "data" in d and isinstance(d["data"], dict):
        d = d["data"]
    return d


def load_pooled(path):
    """Return (asv_set, assign) for one pooled record.

    asv_set : frozenset of ASV sequences (the inferred centers).
    assign  : {unique_sequence -> center_sequence or None}, the partition.
    counts  : {unique_sequence -> pooled count}.
    """
    d = _load(path)
    asvs = d.get("asvs")
    uniques = d.get("uniques")
    amap = d.get("map")
    if asvs is None or uniques is None or amap is None:
        sys.exit(f"{path}: not a dada-pooled record (need asvs/uniques/map)")
    if len(uniques) != len(amap):
        sys.exit(f"{path}: uniques ({len(uniques)}) != map ({len(amap)})")
    asv_seq = [a["sequence"].upper() for a in asvs]
    asv_set = frozenset(asv_seq)
    assign, counts = {}, {}
    for u, m in zip(uniques, amap):
        s = u["sequence"].upper()
        counts[s] = u.get("count", 1)
        assign[s] = asv_seq[m] if (m is not None and 0 <= m < len(asv_seq)) else None
    return asv_set, assign, counts


def err_maxdelta(path_a, path_b):
    """Max abs difference between two runs' learned err_out matrices."""
    a, b = _load(path_a), _load(path_b)
    ea, eb = a.get("err_out"), b.get("err_out")
    if ea is None or eb is None:
        return None
    if len(ea) != len(eb) or any(len(x) != len(y) for x, y in zip(ea, eb)):
        return float("inf")  # shape mismatch => not comparable
    m = 0.0
    for ra, rb in zip(ea, eb):
        for x, y in zip(ra, rb):
            m = max(m, abs(x - y))
    return m


def diff(base, other):
    """Diff two (asv_set, assign, counts) triples. Returns a result dict."""
    b_asv, b_assign, b_counts = base
    o_asv, o_assign, o_counts = other
    added = sorted(o_asv - b_asv)
    removed = sorted(b_asv - o_asv)

    # Partition churn: uniques present in both whose center assignment changed.
    shared = b_assign.keys() & o_assign.keys()
    churn = [s for s in shared if b_assign[s] != o_assign[s]]
    multi = [s for s in churn if b_counts.get(s, 1) > 1]
    singl = [s for s in churn if b_counts.get(s, 1) == 1]

    def state(seq, assign):
        """Human label for a unique's map state. A null map is unassigned/failed;
        a unique mapped to ITSELF is a center (an ASV)."""
        c = assign[seq]
        if c is None:
            return "failed"
        if c == seq:
            return "center"
        return "absorbed"

    def classify(s):
        return f"{state(s, b_assign)}->{state(s, o_assign)}"

    multi_detail = sorted(
        ({"seq": s[:24] + "...", "count": b_counts.get(s, 1), "change": classify(s)}
         for s in multi),
        key=lambda r: -r["count"],
    )

    # Abundance impact: total reads carried by uniques whose assignment changed.
    total_reads = sum(b_counts.values())
    reads_moved = sum(b_counts.get(s, 1) for s in churn)
    reads_moved_pct = (100.0 * reads_moved / total_reads) if total_reads else 0.0

    return {
        "asv_added": added,
        "asv_removed": removed,
        "asv_identical": not added and not removed,
        "churn_total": len(churn),
        "churn_multiread": len(multi),
        "churn_singleton": len(singl),
        "reads_moved": reads_moved,
        "reads_moved_pct": reads_moved_pct,
        "multiread_detail": multi_detail[:20],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, metavar="LABEL=POOLED",
                    help="baseline band, e.g. 16=out/pooled_band16.json.gz")
    ap.add_argument("--compare", action="append", default=[], metavar="LABEL=POOLED",
                    help="a tightened band's pooled record; repeatable")
    ap.add_argument("--errors", action="append", default=[], metavar="LABEL=ERRJSON",
                    help="optional error model per label for the err_out delta")
    ap.add_argument("--max-reads-moved-pct", type=float, default=0.01, metavar="PCT",
                    help="pass ceiling on the %% of reads that change assignment "
                         "(default 0.01); ASV-set identity is always required")
    ap.add_argument("--json", metavar="PATH", help="write machine-readable report")
    args = ap.parse_args()

    def split(spec):
        if "=" not in spec:
            sys.exit(f"expected LABEL=PATH, got {spec!r}")
        return spec.split("=", 1)

    b_label, b_path = split(args.baseline)
    base = load_pooled(b_path)
    err_paths = dict(split(s) for s in args.errors)

    print(f"baseline {b_label}: {len(base[0])} ASVs, {len(base[1])} uniques\n")
    report = {"baseline": b_label, "n_asv": len(base[0]), "comparisons": []}
    all_pass = True

    for spec in args.compare:
        label, path = split(spec)
        other = load_pooled(path)
        r = diff(base, other)
        r["label"] = label
        r["n_asv"] = len(other[0])
        if b_label in err_paths and label in err_paths:
            r["err_out_maxdelta"] = err_maxdelta(err_paths[b_label], err_paths[label])

        # Pass = ASV catalogue unchanged AND abundance impact under the ceiling.
        # (Raw multi-read churn is reported but not gated: a Hamming-1 member
        # reassignment that leaves the ASV set intact is benign.)
        ok = r["asv_identical"] and r["reads_moved_pct"] <= args.max_reads_moved_pct
        r["pass"] = ok
        all_pass &= ok
        status = "PASS" if ok else "FAIL"
        print(f"{label} vs {b_label}: {status}")
        print(f"  ASVs: {r['n_asv']}  added={len(r['asv_added'])} "
              f"removed={len(r['asv_removed'])}  identical={r['asv_identical']}")
        print(f"  reads moved: {r['reads_moved']} ({r['reads_moved_pct']:.4f}%)  "
              f"[ceiling {args.max_reads_moved_pct}%]")
        print(f"  partition churn: total={r['churn_total']} "
              f"multi-read={r['churn_multiread']} singleton={r['churn_singleton']}")
        if "err_out_maxdelta" in r:
            d = r["err_out_maxdelta"]
            print(f"  err_out max|delta| = {d:.3e}  (informational)")
        for m in r["multiread_detail"]:
            print(f"    MULTI-READ CHURN count={m['count']:>5} {m['change']:>16}  {m['seq']}")
        print()
        report["comparisons"].append(r)

    report["all_pass"] = all_pass
    print("=" * 48)
    print("OVERALL: " + ("PASS — tighter band(s) preserve the ASV set within "
                         f"{args.max_reads_moved_pct}% read drift"
                          if all_pass else
                          "FAIL — an ASV was gained/lost or read drift exceeded "
                         f"{args.max_reads_moved_pct}%"))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"\nwrote {args.json}")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
