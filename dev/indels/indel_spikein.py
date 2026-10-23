#!/usr/bin/env python3
"""indel_spikein.py — spike known indel variants into a real run (issue #291).

Takes a denoised sample and writes synthetic copies of it in which a fraction of
one parent ASV's reads carry a fixed insertion or deletion. Edited reads are
real reads, errors included, so the run's error profile is kept; reads are
MOVED from parent to variant, so depth is unchanged. Each variant's sequence and
read count are then known exactly.

Read -> ASV assignment is exact: `derep --show-map` (read -> unique) chained
with the dada `map` (unique -> cluster). Run both on the same FASTQ.

Design grid (one cell per variant):
  plain  type {ins, del} x size x position {5p, mid, 3p} x fraction
         (sites chosen so the edit is not a homopolymer-length change)
  hp     type {ins, del}, size 1, the run of length >= --hp-min-run nearest
         mid-amplicon, x --hp-fractions: the controls that must NOT split
Each synthetic sample uses each parent at most once, so variants of different
parents do not interact; cells are packed into as many samples as needed.

An edited read must contain the parent's sequence around the site exactly
(--anchor nt each side, plus the deleted bases), so a read with an error there
is not eligible; when too few are, the cell moves what it can and says so.
Inserted bases take the lower quality of the two flanking bases.

Outputs, in --out-dir:
  <prefix>_sNN.fastq.gz    synthetic sample NN
  <prefix>_sNN.truth.fasta parent ASVs + that sample's variants, for reference-eval
  truth.tsv                one row per variant

Usage:
  dada2-rs derep prep/S.fastq.gz --show-map -o S.derep.json
  dada2-rs dada prep/S.fastq.gz --error-model err.json ... -o S.dada.json
  indel_spikein.py --fastq prep/S.fastq.gz --derep S.derep.json --dada S.dada.json \\
      --parents S.seqtab_nochim.json --out-dir spikes/ --seed 1
"""

import argparse
import csv
import gzip
import json
import os
import random
import sys

def read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline().rstrip("\n")
            if not h:
                return
            s = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            yield h, s, q


def load_parents(path):
    """Sequences from a seqtab JSON or a FASTA."""
    if path.endswith(".json"):
        return set(json.load(open(path))["sequences"])
    seqs, cur = set(), []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur:
                seqs.add("".join(cur))
            cur = []
        elif line:
            cur.append(line.upper())
    if cur:
        seqs.add("".join(cur))
    return seqs


def runs(seq):
    """(start, length, base) of every homopolymer run."""
    out, start = [], 0
    for i in range(1, len(seq) + 1):
        if i == len(seq) or seq[i] != seq[start]:
            out.append((start, i - start, seq[start]))
            start = i
    return out


def is_hp_change(parent, s, kind, bases):
    """True when the edit only lengthens or shortens a run of one base."""
    if len(set(bases)) != 1:
        return False
    end = s + len(bases) if kind == "del" else s
    flanks = (parent[s - 1] if s > 0 else None, parent[end] if end < len(parent) else None)
    return bases[0] in flanks


def random_insert(rng, k, left, right):
    """Random k-mer that does not extend either flank and has no run > 3."""
    while True:
        x = "".join(rng.choice("ACGT") for _ in range(k))
        if x[0] == left or x[-1] == right:
            continue
        if any(r[1] > 3 for r in runs(x)):
            continue
        return x


def choose_site(rng, parent, kind, k, position, anchor):
    """Site s (0-based) and inserted bases for a plain (non-homopolymer) edit."""
    n = len(parent)
    target = {"5p": 30, "mid": n // 2, "3p": n - 30 - k}[position]
    lo, hi = anchor, n - anchor - k
    for d in range(0, n):
        for s in (target - d, target + d):
            if not lo <= s <= hi:
                continue
            if kind == "del":
                seg = parent[s : s + k]
                if parent[s - 1] == parent[s] or parent[s + k - 1] == parent[s + k]:
                    continue
                if is_hp_change(parent, s, "del", seg):
                    continue
                return s, seg
            if parent[s - 1] == parent[s]:
                continue
            return s, random_insert(rng, k, parent[s - 1], parent[s])
    raise ValueError(f"no {kind} site of size {k} near {position}")


def choose_hp_site(parent, kind, min_run, anchor):
    """The run of length >= min_run nearest mid-amplicon; edit it by one base."""
    mid = len(parent) // 2
    cand = [
        (abs(st + ln // 2 - mid), st, ln, b)
        for st, ln, b in runs(parent)
        if ln >= min_run and st >= anchor and st + ln + anchor <= len(parent)
    ]
    if not cand:
        raise ValueError(f"no homopolymer run >= {min_run}")
    _, st, ln, b = min(cand)
    return st, b, ln


def apply_edit(parent, s, kind, bases):
    if kind == "del":
        return parent[:s] + parent[s + len(bases) :]
    return parent[:s] + bases + parent[s:]


def edit_read(seq, qual, parent, s, kind, bases, anchor):
    """Apply the edit to one read via an exact anchor window; None if absent."""
    k = len(bases) if kind == "del" else 0
    window = parent[s - anchor : s + k + anchor]
    i = seq.find(window)
    if i < 0 or seq.find(window, i + 1) >= 0:
        return None
    cut = i + anchor
    if kind == "del":
        return seq[:cut] + seq[cut + k :], qual[:cut] + qual[cut + k :]
    qi = min(qual[cut - 1], qual[cut])
    return seq[:cut] + bases + seq[cut:], qual[:cut] + qi * len(bases) + qual[cut:]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--derep", required=True, help="derep JSON with --show-map")
    ap.add_argument("--dada", required=True, help="dada JSON from the same FASTQ")
    ap.add_argument("--parents", help="seqtab JSON or FASTA of eligible parents (e.g. non-chimeric)")
    ap.add_argument("--min-parent-reads", type=int, default=500)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="spike")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sizes", default="1,3,10,50,100")
    ap.add_argument("--positions", default="5p,mid,3p")
    ap.add_argument("--fractions", default="0.005,0.01,0.02,0.05,0.1")
    ap.add_argument("--hp-fractions", default="0.001,0.005,0.01,0.02")
    ap.add_argument("--hp-min-run", type=int, default=4)
    ap.add_argument("--anchor", type=int, default=20)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    derep = json.load(open(a.derep))
    dada = json.load(open(a.dada))
    reads = list(read_fastq(a.fastq))
    if len(reads) != len(derep["map"]):
        sys.exit(f"{a.fastq}: {len(reads)} reads but derep map has {len(derep['map'])}")
    if len(derep["uniques"]) != len(dada["map"]):
        sys.exit("derep and dada disagree on the unique count; run both on the same FASTQ")
    asv_seq = [x["sequence"] for x in dada["asvs"]]
    read_asv = [dada["map"][u] for u in derep["map"]]

    members = {}
    for i, c in enumerate(read_asv):
        if c is not None:
            members.setdefault(c, []).append(i)
    allowed = load_parents(a.parents) if a.parents else set(asv_seq)
    parents = sorted(
        c for c, m in members.items() if asv_seq[c] in allowed and len(m) >= a.min_parent_reads
    )
    if not parents:
        sys.exit("no eligible parents")

    cells = []
    for kind in ("ins", "del"):
        for k in map(int, a.sizes.split(",")):
            for pos in a.positions.split(","):
                for f in map(float, a.fractions.split(",")):
                    cells.append(("plain", kind, k, pos, f))
        for f in map(float, a.hp_fractions.split(",")):
            cells.append(("hp", kind, 1, "mid", f))
    rng.shuffle(cells)

    # Pack cells into samples, each parent at most once per sample.
    samples = []
    for cell in cells:
        for smp in samples:
            free = [p for p in parents if p not in smp]
            if free:
                smp[rng.choice(free)] = cell
                break
        else:
            samples.append({rng.choice(parents): cell})

    os.makedirs(a.out_dir, exist_ok=True)
    truth = []
    for si, smp in enumerate(samples, 1):
        name = f"{a.prefix}_s{si:02d}"
        edited = {}
        variants = []
        for parent_c, (context, kind, k, pos, f) in sorted(smp.items()):
            parent = asv_seq[parent_c]
            if context == "hp":
                st, b, ln = choose_hp_site(parent, kind, a.hp_min_run, a.anchor)
                s, bases, hp_len = (st, b, ln)
            else:
                s, bases = choose_site(rng, parent, kind, k, pos, a.anchor)
                hp_len = 0
            variant = apply_edit(parent, s, kind, bases)
            if variant in asv_seq:
                sys.exit(f"{name}: variant of ASV {parent_c} equals an existing ASV")
            pool = members[parent_c]
            want = max(1, round(f * len(pool)))
            ok = [i for i in pool if edit_read(reads[i][1], reads[i][2], parent, s, kind, bases, a.anchor)]
            chosen = rng.sample(ok, min(want, len(ok)))
            for i in chosen:
                edited[i] = edit_read(reads[i][1], reads[i][2], parent, s, kind, bases, a.anchor)
            vid = f"{name}_v{len(variants) + 1:02d}"
            variants.append((vid, variant))
            truth.append({
                "sample": name, "variant": vid, "parent_asv": parent_c,
                "parent_reads": len(pool), "context": context, "type": kind,
                "size": len(bases), "position": pos, "site": s, "hp_len": hp_len,
                "fraction": f, "reads_wanted": want, "reads_eligible": len(ok),
                "reads_moved": len(chosen), "capped": int(len(chosen) < want),
                "bases": bases, "sequence": variant,
            })
        with gzip.open(os.path.join(a.out_dir, f"{name}.fastq.gz"), "wt", compresslevel=4) as out:
            for i, (h, s_, q) in enumerate(reads):
                if i in edited:
                    s_, q = edited[i]
                    h += " spike"
                out.write(f"{h}\n{s_}\n+\n{q}\n")
        with open(os.path.join(a.out_dir, f"{name}.truth.fasta"), "w") as fa:
            for c in parents:
                fa.write(f">asv{c + 1}\n{asv_seq[c]}\n")
            for vid, v in variants:
                fa.write(f">{vid}\n{v}\n")
        print(f"{name}: {len(smp)} variants, {len(edited)} reads edited", file=sys.stderr)

    with open(os.path.join(a.out_dir, "truth.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(truth[0]), delimiter="\t")
        w.writeheader()
        w.writerows(truth)
    capped = sum(t["capped"] for t in truth)
    print(f"{len(truth)} variants in {len(samples)} samples; {capped} capped by eligible reads",
          file=sys.stderr)


if __name__ == "__main__":
    main()
