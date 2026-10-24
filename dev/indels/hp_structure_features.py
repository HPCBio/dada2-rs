#!/usr/bin/env python3
"""hp_structure_features.py — structure covariates for each homopolymer run (issue #291).

Joins hp_variant_counts.py tables with three descriptions of whether each run
sits in base-paired structure, to test structure as a slippage covariate:

  rfam_paired, rfam_junction   consensus structure from an Infernal CM
                               (cmalign; e.g. Rfam RF00177 for bacterial SSU):
                               fraction of the run's bases paired, with a
                               complementary partner in this sequence; whether
                               the run spans paired and unpaired positions
  plf37, plf72                 mean pairing probability over the run from
                               RNAplfold with DNA parameters, at 37 and 72 C
                               (local folding; needs no consensus model)
  ir_stem                      longest perfect inverted-repeat stem (>= 4 nt,
                               loop 3-30) covering at least 2 run bases

The Rfam columns are optional (--cm); the other two work for any amplicon.

Usage:
  hp_structure_features.py --dataset NAME SEQTAB HP_TSV MIN_READS [--dataset ...] \\
      [--cm RF00177.cm] [--dna-par dna_mathews2004.par] --workdir DIR --out features.tsv
"""

import argparse
import collections
import csv
import glob
import json
import os
import subprocess

COMP = str.maketrans("ACGT", "TGCA")
WC = {("G", "C"), ("C", "G"), ("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")}


def parents_of(seqtab, min_reads):
    t = json.load(open(seqtab))
    reads = {s: sum(r[i] for r in t["counts"]) for i, s in enumerate(t["sequences"])}
    return [s for s, n in sorted(reads.items(), key=lambda x: -x[1]) if n >= min_reads]


def rfam_pairing(fasta, cm, workdir):
    """Per sequence, per position: 1 if paired in the consensus with a
    complementary partner in that sequence, else 0."""
    sto = os.path.join(workdir, "parents.sto")
    subprocess.run(["cmalign", "--cpu", "8", "--outformat", "Pfam", "-o", sto, cm, fasta],
                   check=True, stdout=subprocess.DEVNULL)
    rows, ss = {}, ""
    for line in open(sto):
        if line.startswith("#=GC SS_cons"):
            ss = line.split()[-1]
        elif line.strip() and not line.startswith(("#", "//")):
            k, v = line.split()
            rows[k] = v
    pair, stacks = {}, collections.defaultdict(list)
    opn = {"<": ">", "(": ")", "[": "]", "{": "}"}
    cls = {v: k for k, v in opn.items()}
    for i, ch in enumerate(ss):
        if ch in opn or (ch.isalpha() and ch.isupper()):
            stacks[ch].append(i)
        elif ch in cls:
            j = stacks[cls[ch]].pop()
            pair[i], pair[j] = j, i
        elif ch.isalpha() and ch.islower() and stacks[ch.upper()]:
            j = stacks[ch.upper()].pop()
            pair[i], pair[j] = j, i
    out = {}
    for k, aln in rows.items():
        flags = []
        for c, ch in enumerate(aln):
            if ch in ".-":
                continue
            p = pair.get(c)
            flags.append(int(p is not None and (ch.upper(), aln[p].upper()) in WC))
        out[k] = flags
    return out


def plfold(fasta, par, temp, workdir):
    """Per sequence, per position: probability of being paired (1 - unpaired)."""
    d = os.path.join(workdir, f"plfold{temp}")
    os.makedirs(d, exist_ok=True)
    with open(fasta) as fh:
        subprocess.run(["RNAplfold", "-W", "150", "-L", "100", "-u", "1", "-T", str(temp), "-P", par],
                       stdin=fh, cwd=d, check=True, stdout=subprocess.DEVNULL)
    out = {}
    for f in glob.glob(os.path.join(d, "*_lunp")):
        k = os.path.basename(f)[:-len("_lunp")]
        probs = []
        for line in open(f):
            if line.lstrip().startswith("#"):
                continue
            x = line.split()
            probs.append(1.0 - float(x[1]) if x[1] != "NA" else 0.0)
        out[k] = probs
    return out


def ir_stem(seq, pos, n, kmin=4, kmax=12, loop=(3, 30)):
    best = 0
    for k in range(kmin, kmax + 1):
        for i in range(max(0, pos - k + 2), min(pos + n - 2, len(seq) - k) + 1):
            stem = seq[i:i + k]
            rc = stem.translate(COMP)[::-1]
            lo, hi = i + k + loop[0], i + k + loop[1] + k
            down = seq.find(rc, lo, min(len(seq), hi))
            ulo, uhi = max(0, i - loop[1] - k), max(0, i - loop[0])
            up = seq.find(rc, ulo, uhi)
            if down >= 0 or up >= 0:
                best = k
                break
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", nargs=4, action="append", required=True,
                    metavar=("NAME", "SEQTAB", "HP_TSV", "MIN_READS"))
    ap.add_argument("--cm")
    ap.add_argument("--dna-par", default="/opt/homebrew/share/ViennaRNA/dna_mathews2004.par")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    seqs, tables = {}, []
    for name, seqtab, hp, mr in a.dataset:
        ps = parents_of(seqtab, int(mr))
        for i, s in enumerate(ps):
            seqs[f"{name}_p{i + 1}"] = s
        tables.append((name, hp))
    fasta = os.path.join(a.workdir, "parents.fa")
    with open(fasta, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")
    rfam = rfam_pairing(fasta, a.cm, a.workdir) if a.cm else {}
    pl = {t: plfold(fasta, a.dna_par, t, a.workdir) for t in (37, 72)}

    rows = []
    for name, hp in tables:
        for r in csv.DictReader(open(hp), delimiter="\t"):
            k = f"{name}_{r['parent']}"
            s = seqs[k]
            pos, n = int(r["run_pos"]), int(r["run_len"])
            f = dict(dataset=name, **r)
            if rfam:
                fl = rfam[k][pos:pos + n]
                f["rfam_paired"] = f"{sum(fl) / n:.3f}"
                f["rfam_junction"] = int(0 < sum(fl) < n)
            for t in (37, 72):
                p = pl[t].get(k)
                f[f"plf{t}"] = f"{sum(p[pos:pos + n]) / n:.3f}" if p else "nan"
            f["ir_stem"] = ir_stem(s, pos, n)
            rows.append(f)
    with open(a.out, "w", newline="") as fh:
        fields = list(dict.fromkeys(k for r in rows for k in r))
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(seqs)} parents, {len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
