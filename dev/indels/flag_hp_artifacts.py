#!/usr/bin/env python3
"""flag_hp_artifacts.py — flag homopolymer-length and damage artifacts in a sequence table.

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

Candidates are found without alignment, on each ASV's run-collapsed sequence
(every run reduced to one base):
  run lengths  ASVs that differ only in run lengths share a collapsed sequence.
  one new run  ASVs whose collapsed sequences are one insertion apart differ by
               one indel outside a run (an A inserted between Gs: GGGG -> GGAGG)
               plus any run-length changes. Found by hashing every one-base
               deletion of each collapsed sequence; a deletion that joins two
               equal neighbours (GAG -> G) is collapsed again.
Each ASV's parent is the most abundant ASV, from either search, within
--max-edits edited bases of it: the dominant sequence an artifact would come
from. Edits of opposite sign within --pair-window bases are a substitution
(TGG -> TTG is +1/-1 on adjacent runs; a base moved between two runs of it
likewise), so such pairs are not candidates, nor are ASVs with other
substitutions. A new run's site is keyed like a run's, with run length 0.

Damage. Cytosine deamination in the template reads as C>T (G>A on the other
strand). A damaged molecule is copied faithfully, so every read of it reports
the change at high quality and the denoiser keeps it at a fraction of a percent
of its parent; on PacBio mocks these are 96-100% of the single-substitution
variants below 1% of a parent, at 0.03-0.2%, spread along the gene. A
candidate differs from a more abundant ASV by one C>T or G>A substitution, at
or below --damage-max-frac of it. Found without alignment: two same-length
sequences one substitution apart share their first or their second half. Damage
hits many sites of the same template, so the other such candidates of the same
parent ("siblings") are the evidence, as recurrence is for slippage.

Labels (provisional; every piece of evidence is a column):
  artifact        parent >= --min-fold more abundant, every edit's context
                  explained (recurs, or within --max-ratio of the context
                  rate), and tracking does not argue against it
  real-candidate  varies independently of its parent (orphans, or
                  dispersion > 10), or far above its context's rate (or
                  above --min-unexplained with no context evidence)
  damage          a damage candidate with >= --damage-min-siblings such
                  candidates of its parent (itself included), and tracking
                  does not argue against it
  unresolved      otherwise; single-sample runs lack the tracking evidence

Inputs:  a dada2-rs sequence table (make-sequence-table or
         remove-bimera-denovo JSON, optionally .gz).
Outputs: --out TSV, one row per candidate. With --out-seqtab, a table with
         artifact and damage ASVs removed (--action remove) or their reads
         added to their parent (--action merge); the default only flags.

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


def collapsed(r: list[tuple[str, int, int]]) -> str:
    return "".join(b for b, _, _ in r)


_MOD = (1 << 61) - 1
_BASE = 1_000_003


def prefix_hashes(c: str) -> list[int]:
    h = [0] * (len(c) + 1)
    for i, ch in enumerate(c):
        h[i + 1] = (h[i] * _BASE + ord(ch)) % _MOD
    return h


def one_insertion_pairs(collapsed_of: dict[str, str]) -> list[tuple[str, str, int, int]]:
    """(longer, shorter, i, j) for ASV pairs whose collapsed sequences differ
    by one inserted character: removing longer[i:j] (j = i + 2 when the removed
    character's neighbours are equal and join) leaves shorter."""
    by_hash = collections.defaultdict(list)
    for s, c in collapsed_of.items():
        by_hash[prefix_hashes(c)[-1], len(c)].append(s)
    maxlen = max((len(c) for c in collapsed_of.values()), default=0)
    pw = [1] * (maxlen + 2)
    for i in range(1, len(pw)):
        pw[i] = pw[i - 1] * _BASE % _MOD
    out = []
    for s, c in collapsed_of.items():
        h, n = prefix_hashes(c), len(c)
        for i in range(n):
            j = i + 2 if 0 < i < n - 1 and c[i - 1] == c[i + 1] else i + 1
            tail = (h[n] - h[j] * pw[n - j]) % _MOD
            key = ((h[i] * pw[n - j] + tail) % _MOD, n - (j - i))
            for q in by_hash.get(key, ()):
                if collapsed_of[q] == c[:i] + c[j:]:
                    out.append((s, q, i, j))
    return out


def one_substitution_pairs(seqs: list[str]) -> list[tuple[str, str, int]]:
    """(a, b, i) for same-length sequence pairs that differ only at position i.
    Such a pair shares its first half or its second half exactly."""
    out, seen = [], set()
    for half in (0, 1):
        buckets = collections.defaultdict(list)
        for s in seqs:
            m = len(s) // 2
            buckets[len(s), s[:m] if half == 0 else s[m:]].append(s)
        for group in buckets.values():
            for x in range(len(group)):
                for y in range(x + 1, len(group)):
                    a, b = group[x], group[y]
                    if (a, b) in seen:
                        continue
                    diff = -1
                    for i, (u, v) in enumerate(zip(a, b)):
                        if u != v:
                            if diff >= 0:
                                diff = -2
                                break
                            diff = i
                    if diff >= 0:
                        seen.add((a, b))
                        out.append((a, b, diff))
    return out


DAMAGE = {("C", "T"), ("G", "A")}


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
    ap.add_argument("--damage-max-frac", type=float, default=0.005,
                    help="damage candidates: child / (child + parent) at most this")
    ap.add_argument("--damage-min-siblings", type=int, default=3,
                    help="damage candidates of one parent needed to call damage")
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
    rl, col = {}, {}
    for s in seqs:
        rl[s] = runs(s)
        col[s] = collapsed(rl[s])
        groups[col[s]].append(s)

    def substitution_like(edits_pos):
        """Edits of opposite sign within --pair-window bases."""
        return any(d1 * d2 < 0 and abs(p1 - p2) <= a.pair_window
                   for x, (p1, d1) in enumerate(edits_pos) for (p2, d2) in edits_pos[x + 1:])

    # Per child, every qualifying parent: (parent, edited bases, edits). An edit
    # is (parent run, child run) for a run-length change, or
    # ("ins"|"del", position in the parent, bases) for a new run.
    options = collections.defaultdict(list)
    for members in groups.values():
        if len(members) < 2:
            continue
        for s in members:
            for p in members:
                if total[p] <= total[s]:
                    continue
                idx = [j for j, (pr, sr) in enumerate(zip(rl[p], rl[s])) if pr[2] != sr[2]]
                if substitution_like([(rl[p][j][1], rl[s][j][2] - rl[p][j][2]) for j in idx]):
                    continue
                nd = sum(abs(rl[s][j][2] - rl[p][j][2]) for j in idx)
                if nd <= a.max_edits:
                    options[s].append((p, nd, [(rl[p][j], rl[s][j]) for j in idx]))
    for lng, sht, i, j in one_insertion_pairs(col):
        split = j == i + 2  # the new run sits inside a run of the shorter sequence
        if split and rl[lng][i - 1][2] + rl[lng][i + 1][2] != rl[sht][i - 1][2]:
            continue  # ...which also changed length: where is ambiguous
        # Where the new run sits, in the shorter sequence's coordinates.
        if split:
            at_sht = rl[sht][i - 1][1] + rl[lng][i - 1][2]
        else:
            at_sht = rl[sht][i][1] if i < len(rl[sht]) else len(sht)
        new = rl[lng][i]
        # Run-length changes on the runs the two share, as (short run, long run).
        shared = []
        for ks in range(len(rl[sht])):
            if split and ks == i - 1:
                continue
            kl = ks if ks < i else ks + (2 if split else 1)
            if rl[lng][kl][2] != rl[sht][ks][2]:
                shared.append((rl[sht][ks], rl[lng][kl]))
        for child, parent in ((lng, sht), (sht, lng)):
            if total[parent] <= total[child]:
                continue
            if child == lng:  # the child has an inserted run
                edits = [(pr, cr) for pr, cr in shared] + [("ins", at_sht, new[0] * new[2])]
                pos = [(pr[1], cr[2] - pr[2]) for pr, cr in shared] + [(at_sht, new[2])]
            else:  # the child lacks one of the parent's runs
                edits = [(cr, pr) for pr, cr in shared] + [("del", new[1], new[0] * new[2])]
                pos = [(cr[1], pr[2] - cr[2]) for pr, cr in shared] + [(new[1], -new[2])]
            nd = sum(abs(d) for _, d in pos)
            if nd <= a.max_edits and not substitution_like(pos):
                options[child].append((parent, nd, edits))
    cands = []
    for s, opts in options.items():
        p, nd, diffs = max(opts, key=lambda o: (total[o[0]], -o[1]))
        cands.append((s, p, nd, diffs))

    # Context of each edit, keyed on the parent sequence.
    k = a.context
    def edits_of(p, diffs):
        out = []
        for d in diffs:
            if d[0] in ("ins", "del"):  # a new run: length 0 in the parent's terms
                kind, at, bases = d
                end = at if kind == "ins" else at + len(bases)
                sign = "+" if kind == "ins" else "-"
                out.append((f"{sign}{len(bases)}", bases, 0, p[max(0, at - k):at], p[end:end + k]))
                continue
            (b, start, n), (_, _, m) = d
            out.append((f"{m - n:+d}", b, n, p[max(0, start - k):start], p[start + n:start + n + k]))
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
    # Damage: one C>T or G>A substitution from a more abundant ASV, at a low
    # fraction of it. ASVs already judged as run-length or indel candidates keep
    # that judgement.
    dmg = {}
    for x, y, i in one_substitution_pairs(seqs):
        for child, parent in ((x, y), (y, x)):
            if child in info or total[parent] <= total[child]:
                continue
            if (parent[i], child[i]) not in DAMAGE:
                continue
            tr = tracking(counts[child], counts[parent])
            if tr["frac"] > a.damage_max_frac:
                continue
            if child not in dmg or total[parent] > total[dmg[child][0]]:
                dmg[child] = (parent, i, tr)
    siblings = collections.Counter(p for p, _, _ in dmg.values())
    for s, (p, i, tr) in dmg.items():
        independent = multi and (tr["orphans"] > 0 or (tr["dispersion"] == tr["dispersion"] and tr["dispersion"] > 10))
        if independent:
            label = "real-candidate"
        elif siblings[p] >= a.damage_min_siblings:
            label = "damage"
        else:
            label = "unresolved"
        ev = f"{p[i]}>{s[i]}@{i} {p[max(0, i - k):i]}|{p[i + 1:i + 1 + k]} siblings={siblings[p]}"
        rows.append(dict(label=label, reads=total[s], parent_reads=total[p], fold=f"{total[p] / total[s]:.1f}",
                         frac=f"{tr['frac']:.4g}", samples=tr["samples"],
                         dispersion=f"{tr['dispersion']:.3g}" if multi else "NA",
                         orphans=tr["orphans"] if multi else "NA", edits=1, evidence=ev,
                         parent=p, sequence=s))

    rows.sort(key=lambda r: -r["reads"])
    with a.out.open("w", newline="") as fh:
        fields = ["label", "reads", "parent_reads", "fold", "frac", "samples", "dispersion", "orphans",
                  "edits", "evidence", "parent", "sequence"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    by = collections.Counter(r["label"] for r in rows)
    print(f"{len(seqs)} ASVs, {len(rows)} candidates: " +
          ", ".join(f"{k} {n} ({sum(r['reads'] for r in rows if r['label'] == k)} reads)" for k, n in by.most_common()) +
          ("" if multi else "; one sample, so no tracking evidence"), file=sys.stderr)

    if a.action != "flag":
        drop = {r["sequence"]: r["parent"] for r in rows if r["label"] in ("artifact", "damage")}
        # Merge into the first kept ancestor, so a dropped parent's reads are not lost.
        for c in drop:
            q = drop[c]
            while q in drop:
                q = drop[q]
            drop[c] = q
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
        print(f"{a.action}: {len(drop)} artifact and damage ASVs -> {a.out_seqtab}", file=sys.stderr)


if __name__ == "__main__":
    main()
