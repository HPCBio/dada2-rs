#!/usr/bin/env python3
"""indel_ab_compare.py — compare the arms of run_indel_ab.sh (issue #291).

Without a truth set, an ASV that only the indel-aware λ calls is judged by how
its abundance relates to its parent's across samples:

  - an error (PCR slippage or sequencing) is drawn from its parent's molecules,
    so it tracks the parent at a roughly constant small fraction: given the
    sample's reads of the pair, its count is close to binomial with one
    fraction for every sample;
  - a real variant varies independently of the parent: the fraction swings
    between samples, it is missing where the parent is deep, or present where
    the parent is absent.

For each new ASV the parent is the nearest default-arm ASV (dada2-rs
reference-eval). Counts of both come from the indel arm's chimera-free table,
the table where the reads have moved. Per new ASV:

  frac          pooled fraction new / (new + parent)
  dispersion    Pearson chi-square / (k - 1) for one fraction across the k
                samples holding either; ~1 if it tracks, >> 1 if not
  orphans       samples with the new ASV (>= 2 reads) but no parent reads
  missing       samples where the pooled fraction predicts >= 5 reads and
                none were seen

Label (provisional, for triage): `independent` if it has orphans or
dispersion > 10; `tracks` if dispersion <= 3 with no orphans; `sparse` below
10 reads or 3 samples with the pair; otherwise `unclear`. PCR stochasticity
alone makes errors somewhat overdispersed, so read dispersion as a scale, not
a test.

Usage:
  indel_ab_compare.py --binary target/release/dada2-rs <out-dir> [--threads 8]
Writes <out-dir>/compare/: new_asvs.tsv, lost_asvs.tsv, and the FASTA and
reference-eval output used for parents.
"""

import argparse
import collections
import csv
import json
import math
import os
import re
import subprocess
import sys


def table(path):
    d = json.load(open(path))
    seqs = d["sequences"]
    counts = {s: [row[i] for row in d["counts"]] for i, s in enumerate(seqs)}
    return d["samples"], counts


def write_fasta(path, seqs, counts):
    with open(path, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">a{i};size={sum(counts[s])}\n{s}\n")


def align(a, b, band=32):
    """Banded global alignment (match 5, mismatch -4, gap -8), as rows of
    (a_char, b_char) with '-' for gaps."""
    n, m = len(a), len(b)
    band = max(band, abs(n - m) + 4)
    neg = -10**9
    score = [dict() for _ in range(n + 1)]
    trace = [dict() for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(max(0, i - band), min(m, i + band) + 1):
            if i == 0 and j == 0:
                score[0][0] = 0
                continue
            best, t = neg, None
            if i > 0 and j > 0 and (j - 1) in score[i - 1]:
                v = score[i - 1][j - 1] + (5 if a[i - 1] == b[j - 1] else -4)
                if v > best:
                    best, t = v, "d"
            if i > 0 and j in score[i - 1]:
                v = score[i - 1][j] - 8
                if v > best:
                    best, t = v, "u"
            if j > 0 and (j - 1) in score[i]:
                v = score[i][j - 1] - 8
                if v > best:
                    best, t = v, "l"
            score[i][j], trace[i][j] = best, t
    i, j, out = n, m, []
    while i > 0 or j > 0:
        t = trace[i][j]
        if t == "d":
            out.append((a[i - 1], b[j - 1])); i -= 1; j -= 1
        elif t == "u":
            out.append((a[i - 1], "-")); i -= 1
        else:
            out.append(("-", b[j - 1])); j -= 1
    return out[::-1]


def run_len(s, i):
    """Length of the homopolymer run containing s[i]."""
    if i < 0 or i >= len(s):
        return 0
    lo = i
    while lo > 0 and s[lo - 1] == s[i]:
        lo -= 1
    hi = i
    while hi + 1 < len(s) and s[hi + 1] == s[i]:
        hi += 1
    return hi - lo + 1


def describe(parent, new):
    """Edits of `new` against `parent`: substitutions, and each indel as
    kind + length, with `@run<n>` when it changes a homopolymer's length: n is
    the largest run involved (the parent's run for a deletion, the lengthened
    run for an insertion), the rule the indel-aware λ uses."""
    al = align(parent, new)
    subs, indels = 0, []
    k, ppos = 0, 0  # ppos: index in parent of the next parent base
    while k < len(al):
        p, q = al[k]
        if p != "-" and q != "-":
            subs += p != q
            ppos += 1
            k += 1
            continue
        kind = "del" if q == "-" else "ins"
        bases = []
        while k < len(al) and (al[k][1] == "-" if kind == "del" else al[k][0] == "-"):
            bases.append(al[k][0] if kind == "del" else al[k][1])
            k += 1
        start = ppos
        if kind == "del":
            ppos += len(bases)
        hp = 0
        if len(set(bases)) == 1:
            b = bases[0]
            if kind == "del":
                hp = run_len(parent, start)
            else:  # an inserted run lengthens a flanking parent run of its base
                for f in (start - 1, start):
                    if 0 <= f < len(parent) and parent[f] == b:
                        hp = run_len(parent, f) + len(bases)
                        break
        indels.append(f"{kind}{len(bases)}{f'@run{hp}' if hp >= 2 else ''}")
    return subs, indels


def tracking(v, p):
    pairs = [(x, y) for x, y in zip(v, p) if x + y > 0]
    nv = sum(x for x, _ in pairs)
    tot = sum(x + y for x, y in pairs)
    frac = nv / tot if tot else 0.0
    chi2 = 0.0
    for x, y in pairs:
        e = (x + y) * frac
        var = (x + y) * frac * (1 - frac)
        if var > 0:
            chi2 += (x - e) ** 2 / var
    k = len(pairs)
    disp = chi2 / (k - 1) if k > 1 else float("nan")
    orphans = sum(1 for x, y in pairs if x >= 2 and y == 0)
    missing = sum(1 for x, y in pairs if x == 0 and (x + y) * frac >= 5)
    cooc = sum(1 for x, y in pairs if x > 0 and y > 0)
    prev = sum(1 for x, _ in pairs if x > 0)
    if nv < 10 or k < 3:
        label = "sparse"
    elif orphans > 0 or disp > 10:
        label = "independent"
    elif disp <= 3:
        label = "tracks"
    else:
        label = "unclear"
    return dict(frac=frac, dispersion=disp, orphans=orphans, missing=missing,
                cooccur=cooc, prevalence=prev, label=label)


def timings(out):
    rows = []
    for arm in ("learn", "default", "indel"):
        d = os.path.join(out, arm)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if not f.startswith("time_"):
                continue
            text = open(os.path.join(d, f)).read()
            wall = rss = None
            m = re.search(r"Elapsed \(wall clock\) time.*: (\S+)", text)
            if m:  # GNU: [h:]mm:ss.ss
                parts = [float(x) for x in m.group(1).split(":")]
                wall = sum(x * 60 ** i for i, x in enumerate(reversed(parts)))
            m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
            if m:
                rss = int(m.group(1)) / 1024 ** 2
            m = re.search(r"([\d.]+) real", text)
            if m:  # BSD
                wall = float(m.group(1))
            m = re.search(r"(\d+)\s+maximum resident set size", text)
            if m:
                rss = int(m.group(1)) / 1024 ** 3
            rows.append((arm, f[5:-4], wall, rss))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", help="run_indel_ab.sh out-dir")
    ap.add_argument("--binary", default="dada2-rs")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    cmp_dir = os.path.join(a.out, "compare")
    os.makedirs(cmp_dir, exist_ok=True)

    samples, dflt = table(os.path.join(a.out, "default", "nochim.json"))
    samples_i, ind = table(os.path.join(a.out, "indel", "nochim.json"))
    if samples != samples_i:
        sys.exit("arms have different samples")
    _, ind_pre = table(os.path.join(a.out, "indel", "seqtab.json"))
    _, dflt_pre = table(os.path.join(a.out, "default", "seqtab.json"))

    new = sorted(set(ind) - set(dflt), key=lambda s: -sum(ind[s]))
    lost = sorted(set(dflt) - set(ind), key=lambda s: -sum(dflt[s]))
    shared = set(ind) & set(dflt)
    reads = lambda t, ss: sum(sum(t[s]) for s in ss)
    print(f"samples {len(samples)}")
    print(f"pre-chimera ASVs: default {len(dflt_pre)}, indel {len(ind_pre)}")
    print(f"chimera-free ASVs: default {len(dflt)}, indel {len(ind)}; shared {len(shared)}, "
          f"new {len(new)} ({reads(ind, new)} reads), lost {len(lost)} ({reads(dflt, lost)} reads)")
    new_pre = set(ind_pre) - set(dflt_pre)
    killed = [s for s in new_pre if s not in ind]
    print(f"indel-only ASVs before chimera removal {len(new_pre)}; removed as bimeras {len(killed)}")

    for arm, d, wall, rss in timings(a.out):
        print(f"time {arm:8s} {d:4s} wall {wall if wall is None else f'{wall:.0f}s'}  "
              f"peak {rss if rss is None else f'{rss:.2f} GB'}")

    if not new:
        return
    ref_fa = os.path.join(cmp_dir, "default_asvs.fa")
    new_fa = os.path.join(cmp_dir, "new_asvs.fa")
    dseqs = sorted(dflt, key=lambda s: -sum(dflt[s]))
    write_fasta(ref_fa, dseqs, dflt)
    write_fasta(new_fa, new, ind)
    per_asv = os.path.join(cmp_dir, "new_vs_default.per_asv.tsv")
    subprocess.run([a.binary, "reference-eval", "--asvs", new_fa, "--reference", ref_fa,
                    "--near-diffs", "100000", "--threads", str(a.threads), "--per-asv", per_asv,
                    "-o", os.path.join(cmp_dir, "new_vs_default.json")],
                   check=True, stderr=subprocess.DEVNULL)
    nearest = {}
    with open(per_asv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            nearest[r["sequence"]] = int(r["best_ref"].split(";")[0][1:])

    rows, labels, kinds = [], collections.Counter(), collections.Counter()
    for s in new:
        parent = dseqs[nearest[s]]
        subs, indels = describe(parent, s)
        t = tracking(ind[s], ind.get(parent, [0] * len(samples)))
        labels[t["label"]] += 1
        kind = "indel-only" if subs == 0 and indels else ("subs-only" if not indels else "mixed")
        kinds[kind] += 1
        rows.append(dict(reads=sum(ind[s]), parent_reads=sum(ind.get(parent, [0])),
                         kind=kind, subs=subs, indels=",".join(indels) or "-",
                         **{k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in t.items()},
                         parent=parent, sequence=s))
    print("new ASVs by edit kind vs parent: " + ", ".join(f"{k} {n}" for k, n in kinds.most_common()))
    print("new ASVs by tracking label: " + ", ".join(f"{k} {n}" for k, n in labels.most_common()))
    hp = collections.Counter()
    for r in rows:
        for e in r["indels"].split(","):
            if e != "-":
                hp[re.sub(r"\d+@", "@", e) if "@run" in e else re.sub(r"\d+$", "", e) + "-nonhp"] += 1
    print("indels in new ASVs: " + ", ".join(f"{k} {n}" for k, n in sorted(hp.items())))
    with open(os.path.join(cmp_dir, "new_asvs.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(cmp_dir, "lost_asvs.tsv"), "w") as fh:
        fh.write("reads\tprevalence\tsequence\n")
        for s in lost:
            fh.write(f"{sum(dflt[s])}\t{sum(1 for x in dflt[s] if x)}\t{s}\n")
    print(f"wrote {cmp_dir}/new_asvs.tsv and lost_asvs.tsv")


if __name__ == "__main__":
    main()
