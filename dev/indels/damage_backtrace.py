#!/usr/bin/env python3
"""damage_backtrace.py — fate of every one-substitution child of a truth allele in a pooled run (issue #291).

Usage:
  dada2-rs dada-pooled ... --pooled-record record.json
  damage_backtrace.py record.json truth.fa err.json out.tsv reads.fastq.gz [...] \\
      [--slice fwd:260 | rev:210] [--overlap N] [--omega-a 1e-40]

Parents are uniques identical to a truth sequence; children are same-length
uniques one substitution from one. Each child's fate comes from the record's
global map: absorbed into its parent's ASV, split into its own, absorbed
elsewhere, or failed. Expected reads = err_out[type][site Q] x the parent
unique's reads, with site Q the rounded mean over the child's reads at that
position (from the FASTQs, which must be the reads the record was built from).
The abundance p-value is dada's Poisson tail conditioned on at least one read;
a child is predicted to split when it falls below --omega-a.

--slice fwd:L uses truth[:L], rev:L uses revcomp(truth)[:L] (Illumina reads
truncated after primer removal). --overlap N marks children whose position
falls in read coordinates >= N (the read-pair overlap).
"""
import argparse, collections, gzip, json, math, statistics

NT = "ACGT"
RC = str.maketrans("ACGT", "TGCA")
GROUP = {"C>T": "deamination", "G>A": "deamination", "G>T": "oxidation", "C>A": "oxidation",
         "A>G": "other transition", "T>C": "other transition"}
GROUPS = ("deamination", "oxidation", "other transition", "other transversion")


def read_fasta(p):
    out, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split(";")[0]
        elif n:
            out[n] = out.get(n, "") + l
    return out


def log10_pval(n, lam):
    """log10 P(X >= n | X >= 1), X ~ Poisson(lam), as dada's abundance p-value."""
    if n <= lam:
        return 0.0
    lp = lambda k: k * math.log(lam) - lam - math.lgamma(k + 1)
    first = lp(n)
    tail, k = 0.0, n
    while True:
        t = math.exp(lp(k) - first)
        tail += t
        if t < 1e-17 * tail:
            break
        k += 1
    return (first + math.log(tail) - math.log(-math.expm1(-lam))) / math.log(10)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("record"); ap.add_argument("truth"); ap.add_argument("err"); ap.add_argument("out")
    ap.add_argument("fastq", nargs="+")
    ap.add_argument("--slice")
    ap.add_argument("--overlap", type=int)
    ap.add_argument("--omega-a", type=float, default=1e-40)
    a = ap.parse_args()

    rec = json.load(gzip.open(a.record, "rt") if a.record.endswith(".gz") else open(a.record))
    useqs = [u["sequence"] for u in rec["uniques"]]
    ucount = [u["count"] for u in rec["uniques"]]
    gmap = rec["map"]
    asv_seq = [x["sequence"] for x in rec["asvs"]]
    err = json.load(open(a.err))["err_out"]

    truth = read_fasta(a.truth)
    if a.slice:
        d, L = a.slice.split(":"); L = int(L)
        truth = {n: (s if d == "fwd" else s.translate(RC)[::-1])[:L] for n, s in truth.items()}
    names = collections.defaultdict(list)
    for n, s in truth.items():
        names[s].append(n)
    idx = {s: i for i, s in enumerate(useqs)}
    parents = {s: idx[s] for s in names if s in idx}
    print(f"{len(names)} distinct truth sequences, {len(parents)} present as uniques", flush=True)

    # one-substitution children by half-bucketing against the parents
    halves = collections.defaultdict(list)
    for s in parents:
        m = len(s) // 2
        halves[(len(s), 0, s[:m])].append(s)
        halves[(len(s), 1, s[m:])].append(s)
    child = {}
    for ui, s in enumerate(useqs):
        if s in parents:
            continue
        m = len(s) // 2
        for key in ((len(s), 0, s[:m]), (len(s), 1, s[m:])):
            for p in halves.get(key, ()):
                diff = [i for i in range(len(s)) if s[i] != p[i]]
                if len(diff) == 1:
                    pi = parents[p]
                    if ui not in child or ucount[pi] > ucount[child[ui][0]]:
                        child[ui] = (pi, diff[0])
    print(f"{len(child)} one-substitution children", flush=True)

    # site qualities from the reads
    want = {useqs[ui]: pos for ui, (_, pos) in child.items()}
    qsum, qn = collections.Counter(), collections.Counter()
    for f in a.fastq:
        with gzip.open(f, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
                pos = want.get(s)
                if pos is not None:
                    qsum[s] += ord(q[pos]) - 33; qn[s] += 1

    log_omega = math.log10(a.omega_a)
    rows = []
    for ui, (pi, pos) in child.items():
        s, p = useqs[ui], useqs[pi]
        t = f"{p[pos]}>{s[pos]}"
        cm, pm = gmap[ui], gmap[pi]
        if cm is None:
            fate = "failed"
        elif cm == pm:
            fate = "absorbed"
        elif asv_seq[cm] == s:
            fate = "split"
        else:
            fate = "elsewhere"
        q = round(qsum[s] / qn[s]) if qn[s] else None
        row = 4 * NT.index(p[pos]) + NT.index(s[pos])
        e = err[row][min(q, len(err[row]) - 1)] if q is not None else float("nan")
        exp = e * ucount[pi]
        lp = log10_pval(ucount[ui], exp) if exp > 0 else float("nan")
        rows.append(dict(parent="+".join(names[p]), pos=pos, type=t, group=GROUP.get(t, "other transversion"),
                         reads=ucount[ui], parent_reads=ucount[pi],
                         frac=ucount[ui] / (ucount[ui] + ucount[pi]), q=q, err=e, expected=exp,
                         obs_exp=ucount[ui] / exp if exp > 0 else float("inf"),
                         log10_p=lp, predicted="split" if lp < log_omega else "join", fate=fate,
                         overlap=(pos >= a.overlap) if a.overlap is not None else ""))
    with open(a.out, "w") as fh:
        k = list(rows[0])
        fh.write("\t".join(k) + "\n")
        for r in sorted(rows, key=lambda r: -r["reads"]):
            fh.write("\t".join(f"{r[x]:.4g}" if isinstance(r[x], float) else str(r[x]) for x in k) + "\n")

    print(f"\n{'group':18s} {'children':>8s} {'reads':>7s} {'split':>6s} {'split%':>7s} {'medQ':>5s} "
          f"{'med err':>9s} {'med frac':>9s} {'med o/e':>8s}")
    for g in GROUPS:
        rs = [r for r in rows if r["group"] == g]
        if not rs:
            continue
        sp = sum(r["fate"] == "split" for r in rs)
        med = lambda f: statistics.median(f(r) for r in rs)
        print(f"{g:18s} {len(rs):8d} {sum(r['reads'] for r in rs):7d} {sp:6d} {100 * sp / len(rs):6.1f}% "
              f"{med(lambda r: r['q'] or 0):5.0f} {med(lambda r: r['err']):9.2e} {med(lambda r: r['frac']):9.2e} "
              f"{med(lambda r: r['obs_exp']):8.2f}")
    print("\nsplit fraction by observed / expected reads (children with >= 2 reads):")
    bins = [(0, 1), (1, 10), (10, 100), (100, 1e9)]
    print(f"{'group':18s} " + " ".join(f"{f'{lo:g}-{hi:g}' if hi < 1e9 else f'>{lo:g}':>14s}" for lo, hi in bins))
    for g in GROUPS:
        cells = []
        for lo, hi in bins:
            rs = [r for r in rows if r["group"] == g and r["reads"] >= 2 and lo <= r["obs_exp"] < hi]
            sp = sum(r["fate"] == "split" for r in rs)
            cells.append(f"{sp}/{len(rs)}")
        print(f"{g:18s} " + " ".join(f"{c:>14s}" for c in cells))
    print("\nfates:", dict(collections.Counter(r["fate"] for r in rows)))
    split = [r for r in rows if r["fate"] == "split"]
    pred = [r for r in rows if r["predicted"] == "split"]
    hit = sum(r["predicted"] == "split" for r in split)
    print(f"predicted splits (p < {a.omega_a:g}): {hit} of {len(split)} observed splits; "
          f"{len(pred) - hit} predicted but not split")
    if a.overlap is not None:
        for ov in (True, False):
            rs = [r for r in rows if r["overlap"] == ov and r["reads"] >= 2]
            print(f"overlap={ov}: {len(rs)} children >=2 reads, split {sum(r['fate'] == 'split' for r in rs)}")


if __name__ == "__main__":
    main()
