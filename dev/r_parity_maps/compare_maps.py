#!/usr/bin/env python3
"""Compare per-sample unique -> cluster maps and transitions, R vs dada2-rs (#277).

Reads r_dada_maps.R's <s>.r.uniques.tsv / <s>.r.trans.tsv and rs_dada_maps.sh's
<s>.rs.derep.json / <s>.rs.dada.json, for every sample present on both sides.

Uniques are joined by sequence, so the two derep orders need not agree. A unique
is "moved" when its cluster centre differs. Transition cells are compared as
counts. With one error model on both sides, any difference is from screening,
alignment or the divisive loop.

Usage: compare_maps.py <r-dir> <rs-dir> [-o moved.tsv]
Exit status 0 when every sample is identical, 1 otherwise.
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import sys

TRANS = [f"{a}2{b}" for a in "ACGT" for b in "ACGT"]


def short(seq):
    return hashlib.md5(seq.encode()).hexdigest()[:8] if seq else "NA"


def read_r(r_dir, s):
    with open(os.path.join(r_dir, f"{s}.r.uniques.tsv")) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    uniq = {r["sequence"]: (int(r["abundance"]), None if r["center"] == "NA" else r["center"])
            for r in rows}
    path = os.path.join(r_dir, f"{s}.r.trans.tsv")
    if not os.path.exists(path):  # r_pooled_maps.R writes none
        return uniq, None
    trans = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r.pop("trans")
            for q, n in r.items():
                if int(n):
                    trans[(t, int(q))] = int(n)
    return uniq, trans


def read_rs(rs_dir, s):
    derep = json.load(open(os.path.join(rs_dir, f"{s}.rs.derep.json")))["uniques"]
    dada = json.load(open(os.path.join(rs_dir, f"{s}.rs.dada.json")))
    asvs = [a["sequence"] for a in dada["asvs"]]
    if len(dada["map"]) != len(derep):
        sys.exit(f"{s}: dada map has {len(dada['map'])} entries, derep has {len(derep)}")
    # map + derep must rebuild every ASV's reads. A derep in another order than
    # the one dada ran on would scramble the map and still join by sequence.
    rebuilt = [0] * len(asvs)
    for u, m in zip(derep, dada["map"]):
        if m is not None:
            rebuilt[m] += u["count"]
    if rebuilt != [a["abundance"] for a in dada["asvs"]]:
        sys.exit(f"{s}: map + derep do not reproduce the ASV abundances; wrong derep?")
    uniq = {u["sequence"]: (u["count"], None if m is None else asvs[m])
            for u, m in zip(derep, dada["map"])}
    aux = dada.get("aux")
    if aux is None:  # dada-pooled: no per-sample transitions
        return uniq, None
    ncol = aux["transitions_ncol"]
    trans = {}
    for i, n in enumerate(aux["transitions"]):
        if n:
            trans[(TRANS[i // ncol], i % ncol)] = n
    return uniq, trans


def model_check(label, model_path, per_sample):
    """Sum per-sample trans and diff against a learned model's trans.

    Run under the model's err_in with OMEGA_C = 0, per-sample dada replays the
    final learning pass, so the sum must equal the model's trans exactly. If it
    does not, the run is not replaying that pass and its maps say nothing about
    learning. Returns 1 on a mismatch.
    """
    if any(t is None for t in per_sample):
        sys.exit(f"{label} model check: a sample has no per-sample trans")
    total = {}
    for t in per_sample:
        for k, n in t.items():
            total[k] = total.get(k, 0) + n
    m = json.load(open(model_path))["trans"]
    model = {(TRANS[i], q): n for i, row in enumerate(m) for q, n in enumerate(row) if n}
    cells = sorted((k, model.get(k, 0), total.get(k, 0)) for k in set(model) | set(total)
                   if model.get(k, 0) != total.get(k, 0))
    print(f"{label} model check: {len(per_sample)} sample(s), sum {sum(total.values()):,} vs "
          f"model {sum(model.values()):,}; {len(cells)} cell(s) differ", file=sys.stderr)
    for (t, q), want, got in cells[:12]:
        print(f"    {t} Q{q}: model {want}, samples {got}", file=sys.stderr)
    return 1 if cells else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("r_dir")
    ap.add_argument("rs_dir")
    ap.add_argument("-o", "--out", help="TSV of moved uniques and differing trans cells")
    ap.add_argument("--r-model", help="err.json whose trans the R side's per-sample trans must sum to")
    ap.add_argument("--rs-model", help="err.json whose trans the dada2-rs side's per-sample trans must sum to")
    a = ap.parse_args()

    r_samples = {os.path.basename(p)[: -len(".r.uniques.tsv")]
                 for p in glob.glob(os.path.join(a.r_dir, "*.r.uniques.tsv"))}
    rs_samples = {os.path.basename(p)[: -len(".rs.dada.json")]
                  for p in glob.glob(os.path.join(a.rs_dir, "*.rs.dada.json"))}
    samples = sorted(r_samples & rs_samples)
    for s in sorted(r_samples ^ rs_samples):
        print(f"WARNING: {s} is on one side only; skipped", file=sys.stderr)
    if not samples:
        sys.exit("no sample on both sides")

    out_rows = []
    n_diff = 0
    print("sample\tuniques\treads\tcentres_r\tcentres_rs\tmoved_uniques\tmoved_reads\ttrans_cells\ttrans_L1")
    for s in samples:
        ru, rt = read_r(a.r_dir, s)
        su, st = read_rs(a.rs_dir, s)
        if set(ru) != set(su) or any(ru[k][0] != su[k][0] for k in ru):
            # Derep must agree before maps mean anything.
            print(f"{s}: DEREP DIFFERS ({len(ru)} vs {len(su)} uniques)", file=sys.stderr)
            n_diff += 1
            continue
        moved = [(seq, ru[seq][0], ru[seq][1], su[seq][1]) for seq in ru if ru[seq][1] != su[seq][1]]
        if rt is None or st is None:
            rt = st = {}  # one side has no per-sample transitions: compare maps only
        cells = sorted((k, rt.get(k, 0), st.get(k, 0)) for k in set(rt) | set(st)
                       if rt.get(k, 0) != st.get(k, 0))
        l1 = sum(abs(r - x) for _, r, x in cells)
        cr = {c for _, c in ru.values() if c}
        cs = {c for _, c in su.values() if c}
        print(f"{s}\t{len(ru)}\t{sum(v[0] for v in ru.values())}\t{len(cr)}\t{len(cs)}\t"
              f"{len(moved)}\t{sum(m[1] for m in moved)}\t{len(cells)}\t{l1}")
        if moved or cells or cr != cs:
            n_diff += 1
        for seq, ab, c_r, c_s in moved:
            out_rows.append([s, "unique", short(seq), ab, short(c_r), short(c_s), seq])
        for (t, q), r, x in cells:
            out_rows.append([s, "trans", t, q, r, x, ""])

    # The model checks cover every sample on that side, not just the shared ones:
    # the model was learned from all of them.
    if a.r_model:
        n_diff += model_check("R", a.r_model, [read_r(a.r_dir, s)[1] for s in sorted(r_samples)])
    if a.rs_model:
        n_diff += model_check("rs", a.rs_model, [read_rs(a.rs_dir, s)[1] for s in sorted(rs_samples)])

    if a.out:
        with open(a.out, "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            # For "unique" rows: id, abundance, R centre, rs centre, sequence.
            # For "trans" rows: transition, quality, R count, rs count.
            w.writerow(["sample", "kind", "key", "abundance_or_q", "r", "rs", "sequence"])
            w.writerows(out_rows)
    print(f"{len(samples)} sample(s); {n_diff} differ", file=sys.stderr)
    print("==> IDENTICAL" if n_diff == 0 else "==> DIFFERS")
    sys.exit(0 if n_diff == 0 else 1)


if __name__ == "__main__":
    main()
