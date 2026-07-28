#!/usr/bin/env python3
"""Strict comparison of `dada-pseudo` against a manual two-round run (issue #100).

Driven by run_pseudo_equivalence.sh. In dada2-rs the two paths should be
IDENTICAL -- dada-pseudo's round 2 is `mark_priors` + `dada_uniques`, exactly
what `dada --prior` does, with a single error model resolved once and shared by
both rounds. So this deliberately compares for equality rather than for
"similar enough": any difference is a bug on our side, and reporting it loudly
is the point.

Three things are checked, in the order that makes a failure diagnosable:

  1. prior sets      -- if these differ, nothing downstream is comparable, and
                        the fault is in prior SELECTION, not denoising
  2. per-sample ASVs -- sequence -> abundance, exactly, per sample
  3. per-sample maps -- the unique -> cluster assignment, exactly

Round 1 is loaded only to report how much the priors actually changed; a run
where round 2 == round 1 proves nothing (it means the priors were inert), so
that case is called out rather than passed silently.

Usage:
  compare_pseudo_equivalence.py --priors-a A.fa --priors-b B.fa \
      --pseudo DIR --manual DIR [--round1 DIR] [--arm-c DIR]

  compare_pseudo_equivalence.py --priors-only A.fa B.fa
      Quiet mode: exit 0 if the two prior sets are identical, 1 otherwise.
      Used by the driver to decide whether arm C is worth running.

Exit status: 0 if every comparison is identical, 1 otherwise.
"""

import argparse
import gzip
import json
import sys
from pathlib import Path


def read_fasta_seqs(path):
    """Sequence set from a FASTA. Ids are ignored: the two arms name priors
    differently (`prior1..` vs a SHA1 id) but the sequences are what matter."""
    seqs = set()
    cur = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.add("".join(cur).upper())
                    cur = []
            else:
                cur.append(line)
    if cur:
        seqs.add("".join(cur).upper())
    return seqs


def load_sample_dir(path):
    """{sample: {"asvs": {seq: abundance}, "map": [...], "total_reads": n}}.

    Keyed by the `sample` field inside each JSON rather than by filename, so a
    naming difference between arms cannot masquerade as a content difference.
    """
    out = {}
    for f in sorted(Path(path).glob("*.json*")):
        opener = gzip.open if f.name.endswith(".gz") else open
        with opener(f, "rt") as fh:
            d = json.load(fh)
        # Files are tagged ({"dada": {...}}); unwrap a single-key envelope.
        if len(d) == 1 and isinstance(next(iter(d.values())), dict):
            inner = next(iter(d.values()))
            if "asvs" in inner:
                d = inner
        name = d.get("sample", f.stem)
        out[name] = {
            "asvs": {a["sequence"].upper(): a["abundance"] for a in d.get("asvs", [])},
            "map": d.get("map", []),
            "total_reads": d.get("total_reads"),
        }
    return out


def diff_asvs(a, b, label_a, label_b, limit=5):
    """Human-readable diff of two {seq: abundance} maps."""
    lines = []
    only_a = set(a) - set(b)
    only_b = set(b) - set(a)
    shared = set(a) & set(b)
    changed = {s: (a[s], b[s]) for s in shared if a[s] != b[s]}

    if only_a:
        lines.append(f"      {len(only_a)} ASV(s) only in {label_a}:")
        for s in sorted(only_a, key=lambda s: -a[s])[:limit]:
            lines.append(f"        abund={a[s]:<8} {s[:60]}...")
    if only_b:
        lines.append(f"      {len(only_b)} ASV(s) only in {label_b}:")
        for s in sorted(only_b, key=lambda s: -b[s])[:limit]:
            lines.append(f"        abund={b[s]:<8} {s[:60]}...")
    if changed:
        lines.append(f"      {len(changed)} shared ASV(s) with different abundance:")
        for s, (x, y) in sorted(changed.items(), key=lambda kv: -abs(kv[1][0] - kv[1][1]))[:limit]:
            lines.append(f"        {label_a}={x:<8} {label_b}={y:<8} {s[:50]}...")
    return lines


def summarize(a, b, label_a, label_b):
    """Aggregate scoreboard across samples.

    Reports reads recovered as well as ASV counts: a change in partitioning can
    move reads between clusters, or add clusters, without recovering any reads
    that were previously dropped. Those are different claims and the ASV count
    alone cannot distinguish them.
    """
    lines = []
    names = sorted(set(a) & set(b))
    asv_a = sum(len(a[n]["asvs"]) for n in names)
    asv_b = sum(len(b[n]["asvs"]) for n in names)
    reads_a = sum(a[n]["total_reads"] or 0 for n in names)
    reads_b = sum(b[n]["total_reads"] or 0 for n in names)
    gained = sum(len(set(b[n]["asvs"]) - set(a[n]["asvs"])) for n in names)
    lost = sum(len(set(a[n]["asvs"]) - set(b[n]["asvs"])) for n in names)
    moved = sum(
        abs(a[n]["asvs"][s] - b[n]["asvs"][s])
        for n in names
        for s in set(a[n]["asvs"]) & set(b[n]["asvs"])
        if a[n]["asvs"][s] != b[n]["asvs"][s]
    )
    ndiff = sum(1 for n in names if a[n]["asvs"] != b[n]["asvs"] or a[n]["map"] != b[n]["map"])

    w = max(len(label_a), len(label_b), 6)
    lines.append(f"    samples differing        : {ndiff}/{len(names)}")
    lines.append(f"    total ASVs {label_a:>{w}}      : {asv_a}")
    lines.append(f"    total ASVs {label_b:>{w}}      : {asv_b}  ({asv_b - asv_a:+d})")
    lines.append(f"    reads recovered {label_a:>{w}} : {reads_a}")
    lines.append(f"    reads recovered {label_b:>{w}} : {reads_b}  ({reads_b - reads_a:+d})")
    lines.append(f"    ASVs only in {label_b:<{w}}    : {gained}")
    lines.append(f"    ASVs only in {label_a:<{w}}    : {lost}")
    lines.append(f"    reads reassigned among shared ASVs: {moved}")
    if reads_a == reads_b and (gained or lost or moved):
        lines.append(
            "    NOTE: identical reads recovered -- the difference is partitioning"
        )
        lines.append(
            "          (reads moved between clusters), NOT extra reads rescued."
        )
    return lines


def compare_dirs(pseudo, manual, label_a, label_b):
    """Returns (ok, report_lines)."""
    lines = []
    ok = True

    names = sorted(set(pseudo) | set(manual))
    missing = [n for n in names if n not in pseudo or n not in manual]
    if missing:
        ok = False
        lines.append(f"  SAMPLE SET DIFFERS: {missing}")

    for name in names:
        if name not in pseudo or name not in manual:
            continue
        pa, ma = pseudo[name], manual[name]
        same_asvs = pa["asvs"] == ma["asvs"]
        same_map = pa["map"] == ma["map"]
        if same_asvs and same_map:
            lines.append(f"  {name}: identical ({len(pa['asvs'])} ASVs, {pa['total_reads']} reads)")
            continue

        ok = False
        lines.append(f"  {name}: DIFFERS")
        lines.append(
            f"      ASVs {len(pa['asvs'])} ({label_a}) vs {len(ma['asvs'])} ({label_b}); "
            f"reads {pa['total_reads']} vs {ma['total_reads']}"
        )
        if not same_asvs:
            lines.extend(diff_asvs(pa["asvs"], ma["asvs"], label_a, label_b))
        if not same_map:
            n = sum(1 for x, y in zip(pa["map"], ma["map"]) if x != y)
            extra = abs(len(pa["map"]) - len(ma["map"]))
            lines.append(f"      map differs at {n} position(s)" + (f" (+{extra} length)" if extra else ""))
    return ok, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--priors-only", nargs=2, metavar=("A", "B"))
    ap.add_argument("--priors-a")
    ap.add_argument("--priors-b")
    ap.add_argument("--pseudo")
    ap.add_argument("--manual")
    ap.add_argument("--round1")
    ap.add_argument("--arm-c")
    ap.add_argument(
        "--expect-differences",
        action="store_true",
        help="Comparing two deliberately different configurations (e.g. "
        "--reestimate-err-between-rounds on vs off) rather than testing "
        "equivalence. Differences are then the finding, not a bug, and the "
        "exit status reports 0 when they are present.",
    )
    args = ap.parse_args()

    if args.priors_only:
        a, b = (read_fasta_seqs(p) for p in args.priors_only)
        return 0 if a == b else 1

    for req in ("priors_a", "priors_b", "pseudo", "manual"):
        if not getattr(args, req):
            ap.error(f"--{req.replace('_', '-')} is required")

    overall = True
    print("=" * 72)
    print("dada-pseudo vs manual two-round equivalence (issue #100)")
    print("=" * 72)

    # ---- 1. prior sets -------------------------------------------------
    pa = read_fasta_seqs(args.priors_a)
    pb = read_fasta_seqs(args.priors_b)
    print(f"\n[1] prior sets: {len(pa)} (dada-pseudo) vs {len(pb)} (manual)")
    if pa == pb:
        print("    identical")
    else:
        overall = False
        print(f"    DIFFER: {len(pa - pb)} only in pseudo, {len(pb - pa)} only in manual")
        print("    -> prior SELECTION differs; see arm C to isolate denoising")
        for s in sorted(pa - pb)[:3]:
            print(f"      pseudo-only: {s[:60]}...")
        for s in sorted(pb - pa)[:3]:
            print(f"      manual-only: {s[:60]}...")

    # ---- 2/3. per-sample ASVs and maps ---------------------------------
    pseudo = load_sample_dir(args.pseudo)
    manual = load_sample_dir(args.manual)
    print(f"\n[2] round-2 output: dada-pseudo vs manual ({len(pseudo)} sample(s))")
    ok, lines = compare_dirs(pseudo, manual, "pseudo", "manual")
    overall &= ok
    print("\n".join(lines))
    print("\n    -- aggregate --")
    print("\n".join(summarize(pseudo, manual, "pseudo", "manual")))

    if args.arm_c:
        armc = load_sample_dir(args.arm_c)
        print("\n[2b] arm C: manual round 2 re-run with dada-pseudo's OWN priors")
        okc, linesc = compare_dirs(pseudo, armc, "pseudo", "armC")
        print("\n".join(linesc))
        if okc:
            print("    -> denoising agrees given identical priors; the fault is prior selection")
        else:
            print("    -> differs even with identical priors; the fault is in round 2 itself")

    # ---- did the priors do anything at all? ----------------------------
    if args.round1:
        r1 = load_sample_dir(args.round1)
        inert = all(
            r1.get(n, {}).get("asvs") == pseudo[n]["asvs"] for n in pseudo if n in r1
        )
        print("\n[3] sanity: did priors change anything vs round 1?")
        if inert:
            print("    NO -- round 2 == round 1 in every sample. The priors were inert,")
            print("    so this run does NOT exercise the prior path. Use more samples,")
            print("    or a dataset with shared low-abundance variants, before")
            print("    concluding anything from a PASS above.")
        else:
            deltas = [
                (n, len(pseudo[n]["asvs"]) - len(r1[n]["asvs"])) for n in sorted(pseudo) if n in r1
            ]
            print("    yes -- ASV count change per sample (round 2 - round 1):")
            for n, d in deltas:
                print(f"      {n}: {d:+d}")

    print("\n" + "=" * 72)
    if args.expect_differences:
        print("RESULT:", "no differences found -- the configurations agree" if overall
              else "DIFFERENCES FOUND -- see the aggregate block above")
        print("=" * 72)
        return 0
    print("RESULT:", "IDENTICAL -- the two paths agree" if overall
          else "DIFFERENCES FOUND -- see above (this is a bug on our side)")
    print("=" * 72)
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
