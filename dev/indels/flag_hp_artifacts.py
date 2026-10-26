#!/usr/bin/env python3
"""flag_hp_artifacts.py — flag homopolymer-length artifacts in a sequence table.

A post-denoising step in the spirit of chimera removal (issue #291): an ASV
that differs from a more abundant ASV only in homopolymer run lengths may be a
PCR or sequencing artifact (slippage) or a real variant (e.g. an rRNA copy
that differs in one run). The denoiser judges one unique against one centre;
the evidence that separates the two lies across the run:

  tracking     an artifact is drawn from its parent's molecules, so it is a
               near-constant fraction of the parent in every sample (Pearson
               dispersion ~1, no samples with the variant but not the parent)
  recurrence   an edit in the same direction at the same site (base, run
               length and context) recurs in other, unrelated parents at a
               similar or higher fraction; edit size is ignored, so a double
               slip (+2) is explained by single slips (+1) at that site
  context rate the fraction is within what that context shows elsewhere in the
               run (--rates: per-run slippage table from hp_variant_counts.py,
               which also counts variants the denoiser absorbed; looked up
               by site, like recurrence)
  fold         the parent is much more abundant than the variant

Candidates are found without alignment: ASVs that differ only in run lengths
have the same sequence once every run is collapsed to one base. Within such a
group each ASV's parent is the most abundant ASV within --max-edits run-length
changes of it: the dominant sequence an artifact would come from. A base
changed to match its neighbour (TGG -> TTG) also leaves the collapsed sequence
unchanged, as +1 and -1 on adjacent runs; such pairs are substitutions, and
so is a +1/-1 pair within --pair-window bases (a base moved between two runs
of it); those ASVs are not candidates, nor are ASVs with other substitutions.

Limitation: an artifact family founded by a non-homopolymer indel (e.g. an
inserted A between Gs) has a different collapsed sequence from its true
parent, so its members are judged against each other; they come out
unresolved or real-candidate (kept), never artifact.

Labels (provisional; every piece of evidence is a column):
  artifact        parent >= --min-fold more abundant, every edit's context
                  explained (recurs, or within --max-ratio of the context
                  rate), and tracking does not argue against it
  real-candidate  varies independently of its parent (orphans, or
                  dispersion > 10), or far above its context's rate (or
                  above --min-unexplained with no context evidence)
  unresolved      otherwise; single-sample runs lack the tracking evidence

Inputs:  a dada2-rs sequence table (make-sequence-table or
         remove-bimera-denovo JSON, optionally .gz).
Outputs: --out TSV, one row per candidate. With --out-seqtab, a table with
         artifacts removed (--action remove) or their reads added to their
         parent (--action merge); the default only flags.

Usage:
  flag_hp_artifacts.py nochim.json --out hp_flags.tsv [--rates hp_runs.tsv]
      [--context 6] [--min-fold 4] [--max-ratio 4] [--max-edits 2]
      [--action flag|remove|merge --out-seqtab filtered.json]
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import statistics
import sys
from pathlib import Path


def read_json(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as raw:
            return json.load(io.TextIOWrapper(raw, encoding="utf-8"))
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def runs(seq: str) -> list[tuple[str, int, int]]:
    """(base, start, length) for every run, length 1 included."""
    out, i = [], 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        out.append((seq[i], i, j - i + 1))
        i = j + 1
    return out


def tracking(v: list[int], p: list[int]) -> dict:
    pairs = [(x, y) for x, y in zip(v, p) if x + y > 0]
    nv, tot = sum(x for x, _ in pairs), sum(x + y for x, y in pairs)
    frac = nv / tot if tot else 0.0
    chi2 = sum((x - (x + y) * frac) ** 2 / ((x + y) * frac * (1 - frac))
               for x, y in pairs if 0 < frac < 1)
    k = len(pairs)
    return dict(frac=frac, samples=k, dispersion=chi2 / (k - 1) if k > 1 else float("nan"),
                orphans=sum(1 for x, y in pairs if x >= 2 and y == 0))


def site(key):
    """An edit's site: direction, base, run length and context, not size."""
    return (key[0][0],) + key[1:]


def load_rates(path: Path, k: int, min_exact: int):
    """Context key -> [(parent md5, frac)] from an hp_variant_counts.py table."""
    rates = collections.defaultdict(list)
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["frac"] in ("", "nan") or int(r["parent_exact"]) < min_exact:
                continue
            if r.get("to_end") and int(r["to_end"]) < 10:
                continue
            key = (r["edit"], r["base"], int(r["run_len"]), r["left"][-k:], r["right"][:k])
            rates[site(key)].append((r.get("parent_md5", ""), float(r["frac"])))
    return rates


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("seqtab", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--rates", type=Path, help="hp_variant_counts.py table for this run")
    ap.add_argument("--context", type=int, default=6, help="bases either side defining a context")
    ap.add_argument("--min-fold", type=float, default=4.0, help="parent / variant abundance for an artifact")
    ap.add_argument("--max-ratio", type=float, default=4.0, help="fraction / context rate still explained")
    ap.add_argument("--max-edits", type=int, default=2, help="run-length changes from the parent")
    ap.add_argument("--pair-window", type=int, default=10,
                    help="a +1/-1 pair this close is substitutions, not slippage")
    ap.add_argument("--min-exact", type=int, default=300, help="--rates rows: minimum error-free parent reads")
    ap.add_argument("--min-unexplained", type=float, default=0.1,
                    help="fraction above which an edit without context evidence counts as real")
    ap.add_argument("--action", choices=["flag", "remove", "merge"], default="flag")
    ap.add_argument("--out-seqtab", type=Path)
    a = ap.parse_args()
    if a.action != "flag" and not a.out_seqtab:
        ap.error("--action remove/merge needs --out-seqtab")

    t = read_json(a.seqtab)
    seqs = t["sequences"]
    counts = {s: [row[i] for row in t["counts"]] for i, s in enumerate(seqs)}
    total = {s: sum(c) for s, c in counts.items()}
    multi = len(t["samples"]) > 1
    rates = load_rates(a.rates, a.context, a.min_exact) if a.rates else None

    groups = collections.defaultdict(list)
    rl = {}
    for s in seqs:
        r = runs(s)
        rl[s] = r
        groups["".join(b for b, _, _ in r)].append(s)

    cands = []
    for members in groups.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda s: -total[s])
        for i, s in enumerate(members[1:], 1):
            for p in members[:i]:  # most abundant first
                if total[p] <= total[s]:
                    continue
                idx = [j for j, (pr, sr) in enumerate(zip(rl[p], rl[s])) if pr[2] != sr[2]]
                delta = {j: rl[s][j][2] - rl[p][j][2] for j in idx}
                if any(delta[i] * delta[j] < 0 and abs(rl[p][j][1] - rl[p][i][1]) <= a.pair_window
                       for i in idx for j in idx if i < j):
                    continue  # a +1/-1 pair close together: substitutions
                nd = sum(abs(d) for d in delta.values())
                if nd <= a.max_edits:
                    cands.append((s, p, nd, [(rl[p][j], rl[s][j]) for j in idx]))
                    break

    # Context of each edit, keyed on the parent sequence.
    k = a.context
    def edits_of(p, diffs):
        out = []
        for (b, start, n), (_, _, m) in diffs:
            key = (f"{m - n:+d}", b, n, p[max(0, start - k):start], p[start + n:start + n + k])
            out.append(key)
        return out
    within = collections.defaultdict(list)  # site -> [(parent, frac)] among candidates
    info = {}
    for s, p, nd, diffs in cands:
        tr = tracking(counts[s], counts[p])
        keys = edits_of(p, diffs)
        info[s] = (p, nd, keys, tr)
        for key in keys:
            within[site(key)].append((p, tr["frac"]))

    rows = []
    for s, (p, nd, keys, tr) in info.items():
        pmd5 = hashlib.md5(p.encode()).hexdigest()
        ev = []
        explained_all, far_above = True, False
        for key in keys:
            others = [f for q, f in within[site(key)] if q != p]
            ctx = [f for (q, f) in rates.get(site(key), []) if q != pmd5] if rates is not None else []
            ctx_rate = statistics.median(ctx) if ctx else None
            explained = (len(others) >= 1 and max(others) * a.max_ratio >= tr["frac"]) or \
                        (ctx_rate is not None and len(ctx) >= 2 and tr["frac"] <= a.max_ratio * max(ctx_rate, 1e-6))
            if ctx_rate is not None and len(ctx) >= 2 and tr["frac"] > a.max_ratio * max(ctx_rate, 1e-6) and not others:
                far_above = True
            if not explained and not others and not ctx and tr["frac"] >= a.min_unexplained:
                far_above = True
            explained_all &= explained
            ev.append(f"{key[0]}{key[1]}{key[2]} {key[3]}|{key[4]} recur={len(others)} "
                      f"ctx={'-' if ctx_rate is None else f'{ctx_rate:.4g}'}(n={len(ctx)})")
        fold = total[p] / total[s]
        independent = multi and (tr["orphans"] > 0 or (tr["dispersion"] == tr["dispersion"] and tr["dispersion"] > 10))
        if independent or far_above:
            label = "real-candidate"
        elif fold >= a.min_fold and explained_all:
            label = "artifact"
        else:
            label = "unresolved"
        rows.append(dict(label=label, reads=total[s], parent_reads=total[p], fold=f"{fold:.1f}",
                         frac=f"{tr['frac']:.4g}", samples=tr["samples"],
                         dispersion=f"{tr['dispersion']:.3g}" if multi else "NA",
                         orphans=tr["orphans"] if multi else "NA", edits=nd, evidence="; ".join(ev),
                         parent=p, sequence=s))
    rows.sort(key=lambda r: -r["reads"])
    with a.out.open("w", newline="") as fh:
        fields = ["label", "reads", "parent_reads", "fold", "frac", "samples", "dispersion", "orphans",
                  "edits", "evidence", "parent", "sequence"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    by = collections.Counter(r["label"] for r in rows)
    print(f"{len(seqs)} ASVs, {len(rows)} homopolymer-length candidates: " +
          ", ".join(f"{k} {n} ({sum(r['reads'] for r in rows if r['label'] == k)} reads)" for k, n in by.most_common()) +
          ("" if multi else "; one sample, so no tracking evidence"), file=sys.stderr)

    if a.action != "flag":
        drop = {r["sequence"]: r["parent"] for r in rows if r["label"] == "artifact"}
        keep = [i for i, s in enumerate(seqs) if s not in drop]
        idx = {s: i for i, s in enumerate(seqs)}
        new_counts = []
        for row in t["counts"]:
            row = list(row)
            if a.action == "merge":
                for s, p in drop.items():
                    row[idx[p]] += row[idx[s]]
            new_counts.append([row[i] for i in keep])
        out = dict(t)
        out["sequences"] = [seqs[i] for i in keep]
        if "sequence_ids" in t:
            out["sequence_ids"] = [t["sequence_ids"][i] for i in keep]
        out["counts"] = new_counts
        a.out_seqtab.write_text(json.dumps(out))
        print(f"{a.action}: {len(drop)} artifact ASVs -> {a.out_seqtab}", file=sys.stderr)


if __name__ == "__main__":
    main()
