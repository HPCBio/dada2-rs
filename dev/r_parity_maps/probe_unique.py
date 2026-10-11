#!/usr/bin/env python3
"""Probe one moved unique from compare_maps.py (#284).

Usage: [CENTRES=a,b] [ERR=err.json] probe_unique.py <maps/r> <moved-seq> <fastq>...
FASTQs are the pool's inputs, in its input order. CENTRES are md5 prefixes from
moved.tsv; ERR is the run's model (an R import for the R side).

1. Every sample containing the sequence, from R's maps, with R's centre.
2. Both centres (any md5 prefixes passed via CENTRES=a,b), and where the
   sequence differs from each.
3. Pooled mean quality per position, two ways: exact integer sums / N (dada2-rs,
   containers.rs from_qual_sums), and R's path (derepFastq mean per sample, then
   combineDereps2 sum of mean*n in input order, / N). Both rounded half away
   from zero, as C round() and Rust f64::round do. Flags positions where the
   rounded values differ, for the moved sequence and both centres.
4. The moved read's quality at each centre mismatch, and with ERR=<err.json>
   its lambda to each centre, multiplied left to right from 1.0 as R's
   compute_lambda and ours do; printed with raw IEEE-754 bits. dada uses the
   unique's rounded mean quality: per-read lambda is dada's only for a
   singleton, as in #284.
"""
import glob, gzip, hashlib, json, math, os, struct, sys

md5 = lambda s: hashlib.md5(s.encode()).hexdigest()[:8]
rnd = lambda x: math.floor(x + 0.5)  # half away from zero for x >= 0

rdir, moved, fqs = sys.argv[1], sys.argv[2], sys.argv[3:]
want_centres = set(filter(None, os.environ.get("CENTRES", "").split(",")))

# 1 + 2: R's maps
centre_seq, centre_reads, present = {}, {}, []
for p in sorted(glob.glob(os.path.join(rdir, "*.r.uniques.tsv"))):
    s = os.path.basename(p).replace(".r.uniques.tsv", "")
    with open(p) as fh:
        next(fh)
        for line in fh:
            seq, ab, cen = line.rstrip("\n").split("\t")
            h = md5(cen)
            centre_reads[h] = centre_reads.get(h, 0) + int(ab)
            if h in want_centres:
                centre_seq[h] = cen
            if seq == moved:
                present.append((s, int(ab), h))
print(f"moved {md5(moved)} len {len(moved)}; in {len(present)} sample(s), pool abundance {sum(a for _, a, _ in present)}")
for s, a, h in present:
    print(f"  {s}: {a} read(s), R centre {h}")
for h in sorted(want_centres):
    c = centre_seq.get(h)
    if c is None:
        print(f"centre {h}: not found in R's maps"); continue
    diffs = [(i, c[i], moved[i]) for i in range(min(len(c), len(moved))) if c[i] != moved[i]]
    print(f"centre {h}: len {len(c)}, {centre_reads[h]} reads in pool (R), "
          f"{len(diffs)} mismatches to moved: {diffs[:12]}")

# 3: qualities, per sample, for the moved sequence and the centres
targets = {md5(moved): moved, **centre_seq}
by_seq = {seq: h for h, seq in targets.items()}
per_sample = {h: [] for h in targets}  # (n, [int sums]) per sample, input order
moved_quals = []  # each read of the moved sequence
for fq in fqs:
    acc = {}
    with gzip.open(fq, "rt") as fh:
        while True:
            hdr = fh.readline()
            if not hdr:
                break
            seq = fh.readline().rstrip(); fh.readline(); q = fh.readline().rstrip()
            h = by_seq.get(seq)
            if h is None:
                continue
            if seq == moved:
                moved_quals.append([ord(c) - 33 for c in q])
            n, sums = acc.get(h, (0, [0] * len(seq)))
            acc[h] = (n + 1, [a + ord(c) - 33 for a, c in zip(sums, q)])
    for h, v in acc.items():
        per_sample[h].append(v)

for h, rows in per_sample.items():
    if not rows:
        print(f"{h}: no reads found in the FASTQs given"); continue
    N = sum(n for n, _ in rows)
    L = len(rows[0][1])
    flips, halves = [], 0
    for i in range(L):
        S = sum(s[i] for _, s in rows)
        exact = S / N
        r_acc = 0.0
        for n, s in rows:
            r_acc = r_acc + (s[i] / n) * n  # derepFastq mean, then sweep(..., "*")
        r_mean = r_acc / N
        halves += (2 * S) % N == 0 and (2 * S // N) % 2 == 1
        if rnd(exact) != rnd(r_mean):
            flips.append((i, exact, r_mean))
    print(f"{h}: {N} reads in {len(rows)} sample(s); {halves} positions with mean exactly k.5; "
          f"{len(flips)} rounding flips {[(i, repr(e), repr(r)) for i, e, r in flips[:8]]}")

# 4: the moved read's qualities at the mismatches, and lambda to each centre
NT = {"A": 0, "C": 1, "G": 2, "T": 3}
bits = lambda x: struct.unpack("<Q", struct.pack("<d", x))[0]
err = json.load(open(os.environ["ERR"]))["err_out"] if os.environ.get("ERR") else None
print(f"moved read(s): {len(moved_quals)}")
for k, q in enumerate(moved_quals):
    for h in sorted(centre_seq):
        c = centre_seq[h]
        mism = [i for i in range(len(c)) if c[i] != moved[i]]
        line = f"  read {k} vs {h}: " + ", ".join(f"pos {i} {c[i]}>{moved[i]} Q{q[i]}" for i in mism)
        if err is not None:
            lam = 1.0
            for i in range(len(moved)):
                row = NT[c[i]] * 4 + NT[moved[i]] if c[i] != moved[i] else NT[moved[i]] * 5
                lam = lam * err[row][q[i]]
            line += f"; lambda {lam!r} bits {bits(lam):016x}"
        print(line)
