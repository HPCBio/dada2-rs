#!/usr/bin/env python3
"""indel_lambda_dryrun.py — offline test of an indel-aware λ (issue #291, Phase 3).

Re-scores the uniques in a `--indel-events` table as if λ carried indel terms,
without re-running dada:

    λ' = λ × Π p(event),  p(event) = p_open × r^(len-1)

  insertions  p_open = per-base insertion rate at the inserted bases' lowest
              quality (`ins_q`)
  deletions   p_open = per-base deletion rate at the centre run length
              (`del_hp`; mixed-base deletions use run length 1)
  r           per type, from the mean length of the rate run's error events,
              assuming lengths fall off geometrically

Rates come from a separate default run (its `--indel-opportunities` and
`--indel-events`), as `--indel-rates` would take them. The guard counts an
event as 1, as today, when the largest homopolymer run it involves is >= --guard:
L for a deletion (L -> L-1), L+1 for an insertion lengthening a run of L, 1
otherwise. Guard 2 thus spares every run-length change. Terminal gaps are
excluded unless --terminal. Counted events add to `hamming`. --flank-q scales
deletion rates by flanking quality; --pair-window scores a nearby insertion and
deletion of equal length as substitutions (the aligner's drawing of them). With
--hp-ceiling,
guarded events are counted too, at that fixed probability instead of 1: a
ceiling above the worst site-specific homopolymer error rate, so recurring slips
stay absorbed but a homopolymer-length variant far above it can split.

With --truth, each bud in a spike sample is the variant itself, its family (a
variant read with one more error, in the parent's cluster, with an event at the
spike site; dada would fold it into the variant once that buds), family in
another cluster (events at a spike site of the same sample, but placed with a
sibling ASV), or an error.

A unique "buds" when dada's abundance test would split it from its cluster:
`p_A × nraw < OMEGA_A` with p_A = P(X >= reads | X >= 1, Poisson(λ' ×
cluster_reads)), nraw the sample's unique count, singletons never tested.
This is the first-bud test only: dada re-partitions after each bud, so this
ranks designs, it does not replace running them.

The events come from a full-NW re-alignment, while `lambda` comes from dada's
own comparison, which may have taken the gapless shortcut (#277/#281). Where the
two disagree (no substitutions in the events, yet lambda < 1e-3: e.g. a run
boundary substitution AAAAGG/AAAGGG drawn as a deletion plus an insertion), the
unique is skipped and counted: in dada the indel terms and λ come from one Sub.

Usage:
  indel_lambda_dryrun.py --rates DIR --events events.tsv --dada dada.json|DIR \\
      [--truth truth.tsv --derep-dir DIR] [--guard 0,1,2,3] [--terminal]
"""

import argparse
import collections
import csv
import gzip
import json
import math
import os

OMEGA_A = 1e-40


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def pooled_rates(rows, table, min_events=20):
    """Per-key rate, widening a window over neighbouring keys until it holds
    min_events events, so sparse keys borrow from their neighbours."""
    opp, ev = collections.Counter(), collections.Counter()
    for r in rows:
        if r["table"] == table:
            k = int(r["key"])
            opp[k] += int(r["opportunities"])
            ev[k] += int(r["events"])
    keys = sorted(opp)
    lo_k, hi_k = keys[0], keys[-1]
    out = {}
    for k in range(lo_k, hi_k + 1):
        w = 0
        while True:
            ks = range(max(lo_k, k - w), min(hi_k, k + w) + 1)
            e, o = sum(ev[x] for x in ks), sum(opp[x] for x in ks)
            if e >= min_events or (k - w <= lo_k and k + w >= hi_k):
                out[k] = e / o if o else 0.0
                break
            w += 1
    return out, lo_k, hi_k


def length_factor(events, kind):
    """r for a geometric length model, from read-weighted error events."""
    n = tot = 0
    for e in events:
        if e["type"] == kind and e["terminal"] == "0" and e["correct"] == "1":
            w = int(e["reads"])
            n += w
            tot += w * int(e["len"])
    mean = tot / n if n else 1.0
    return max(1.0 - 1.0 / mean, 1e-6)


def log_upper_tail(k, lam):
    """log P(X >= k), X ~ Poisson(lam), for lam < k (else None: not small)."""
    if lam <= 0:
        return -math.inf
    if lam >= k:
        return None
    lead = -lam + k * math.log(lam) - math.lgamma(k + 1)
    s, term, j = 1.0, 1.0, 1
    while term > 1e-17 * s and j < 10000:
        term *= lam / (k + j)
        s += term
        j += 1
    return lead + math.log(s)


def buds(reads, lam, cluster_reads, nraw):
    if reads <= 1:
        return False
    e = lam * cluster_reads
    lu = log_upper_tail(reads, e)
    if lu is None:
        return False
    log_norm = math.log(-math.expm1(-e)) if e > 0 else -math.inf
    return lu - log_norm + math.log(nraw) < math.log(OMEGA_A)


def load_dada(dada):
    """Per sample: unique count (nraw) and ASV sequences in cluster order."""
    paths = [dada] if dada.endswith((".json", ".json.gz")) else [
        os.path.join(dada, f) for f in sorted(os.listdir(dada)) if f.endswith((".json", ".json.gz"))]
    nraw, asvs = {}, {}
    for p in paths:
        op = gzip.open if p.endswith(".gz") else open
        d = json.load(op(p, "rt"))
        nraw[d["sample"]] = len(d["map"])
        asvs[d["sample"]] = [x["sequence"] for x in d["asvs"]]
    return nraw, asvs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rates", required=True, help="dir with opp.tsv and events.tsv of a default run")
    ap.add_argument("--events", required=True, help="--indel-events table to re-score")
    ap.add_argument("--dada", required=True, help="dada JSON, or a dir of them, for nraw")
    ap.add_argument("--truth", help="indel_spikein.py truth.tsv: report spike recall")
    ap.add_argument("--derep-dir", help="derep JSONs of the spike samples (unique order)")
    ap.add_argument("--parent-dada", help="dada JSON the spike parents were drawn from")
    ap.add_argument("--guard", default="0,2,3,4", help="run-length thresholds; 0 = no guard")
    ap.add_argument("--hp-ceiling", default="", help="comma list of probabilities for guarded events")
    ap.add_argument("--terminal", action="store_true", help="also count terminal gaps")
    ap.add_argument("--write-recovered", help="TSV of each variant's outcome per guard")
    ap.add_argument("--known", help="FASTA of true sequences: a bud matching one is a real allele")
    ap.add_argument("--write-buds", help="TSV of every bud per setting")
    ap.add_argument("--derep", help="derep JSON of an unspiked sample, for --known")
    ap.add_argument("--flank-q", action="store_true",
                    help="scale deletion rates by the deletion enrichment at the lower flanking quality")
    ap.add_argument("--end-zone", type=int, default=0,
                    help="ignore internal indels within this many nt of the centre's 3' end (read end)")
    ap.add_argument("--pair-window", type=int, default=0,
                    help="score an insertion and a deletion of equal length within this many nt as substitutions")
    a = ap.parse_args()

    rate_rows = read_tsv(os.path.join(a.rates, "opp.tsv"))
    rate_events = read_tsv(os.path.join(a.rates, "events.tsv"))
    p_ins, qlo, qhi = pooled_rates(rate_rows, "ins_q")
    p_del, llo, lhi = pooled_rates(rate_rows, "del_hp")
    r_ins, r_del = length_factor(rate_events, "I"), length_factor(rate_events, "D")
    print(f"rates: r_ins={r_ins:.4f} r_del={r_del:.4f}  "
          f"p_ins Q20={p_ins.get(20, 0):.2e} Q60={p_ins.get(60, 0):.2e} Q90={p_ins.get(90, 0):.2e}  "
          f"p_del L1={p_del.get(1, 0):.2e} L3={p_del.get(3, 0):.2e} L6={p_del.get(min(6, lhi), 0):.2e}")

    p_sub, sqlo, sqhi = pooled_rates(rate_rows, "sub_q")

    # Deletion enrichment by the lower flanking quality, relative to how common
    # bases of that quality are (the rates tables have no flanking dimension).
    edges = [0, 10, 20, 30, 35, 40, 60, 80, 93, 94]
    qbin = lambda q: max(i for i in range(len(edges) - 1) if q >= edges[i])
    base_q = collections.Counter()
    for r in rate_rows:
        if r["table"] == "ins_q":
            base_q[qbin(int(r["key"]))] += int(r["opportunities"])
    del_q = collections.Counter()
    for e in rate_events:
        if e["type"] == "D" and e["terminal"] == "0" and e["correct"] == "1" and e["qual_left"] != "NA" and e["qual_right"] != "NA":
            del_q[qbin(min(int(e["qual_left"]), int(e["qual_right"])))] += int(e["reads"])
    tb, td = sum(base_q.values()), sum(del_q.values())
    flank = {b: (del_q[b] / td) / (base_q[b] / tb) if base_q[b] and td else 1.0 for b in base_q}
    if a.flank_q:
        print("deletion enrichment by flanking-Q bin:", "  ".join(
            f"Q{edges[b]}-{edges[b+1]-1}:{flank[b]:.2g}" for b in sorted(flank)))

    def pair_up(events):
        """Split off insertion/deletion pairs of equal length within --pair-window."""
        if not a.pair_window:
            return events, []
        rest, pairs, used = [], [], set()
        ins = [e for e in events if e["type"] == "I"]
        dels = [e for e in events if e["type"] == "D"]
        for i in ins:
            cand = [d for d in dels if id(d) not in used and d["len"] == i["len"]
                    and abs(int(d["center_pos"]) - int(i["center_pos"])) <= a.pair_window]
            if cand:
                d = min(cand, key=lambda d: abs(int(d["center_pos"]) - int(i["center_pos"])))
                used.add(id(d)); used.add(id(i)); pairs.append((i, d))
        rest = [e for e in events if id(e) not in used]
        return rest, pairs

    def p_pair(pair):
        i, _ = pair
        q = int(i["qual_min"]) if i["qual_min"] != "NA" else sqhi
        return p_sub[min(max(q, sqlo), sqhi)] ** int(i["len"])

    nraw, asvs = load_dada(a.dada)
    by_unique = collections.defaultdict(list)
    for e in read_tsv(a.events):
        by_unique[(e["sample"], e["unique"])].append(e)

    def p_event(e):
        k = int(e["len"])
        if e["type"] == "I":
            q = int(e["qual_min"]) if e["qual_min"] != "NA" else qhi
            return p_ins[min(max(q, qlo), qhi)] * r_ins ** (k - 1)
        L = int(e["hp_len"]) or 1
        p = p_del[min(max(L, llo), lhi)] * r_del ** (k - 1)
        if a.flank_q and e["qual_left"] != "NA" and e["qual_right"] != "NA":
            p *= flank.get(qbin(min(int(e["qual_left"]), int(e["qual_right"]))), 1.0)
        return min(p, 1.0)

    def involved_run(e):
        L = int(e["hp_len"])
        if e["type"] == "D":
            return max(L, 1)
        return L + 1 if L else 1

    variants, spiked, sites, by_variant_unique = [], {}, collections.defaultdict(list), {}
    if a.truth:
        seq_to_unique = {}
        for f in sorted(os.listdir(a.derep_dir)):
            d = json.load(open(os.path.join(a.derep_dir, f)))
            for i, u in enumerate(d["uniques"]):
                seq_to_unique[(d["sample"], u["sequence"])] = str(i)
        parent_seq = [x["sequence"] for x in json.load(open(a.parent_dada))["asvs"]]
        for t in read_tsv(a.truth):
            t["unique"] = seq_to_unique.get((t["sample"], t["sequence"]))
            variants.append(t)
            by_variant_unique[(t["sample"], t["unique"])] = t
            spiked[(t["sample"], parent_seq[int(t["parent_asv"])])] = t
            sites[t["sample"]].append((int(t["site"]) - 5, int(t["site"]) + int(t["size"]) + 5))

    def classify(sample, uid, es):
        # The variant itself, in whichever cluster it was placed.
        if (sample, uid) in by_variant_unique:
            return "variant"
        t = spiked.get((sample, asvs[sample][int(es[0]["cluster"])]))
        if t is not None:
            lo, hi = int(t["site"]) - 5, int(t["site"]) + int(t["size"]) + 5
            if any(lo <= int(e["center_pos"]) <= hi for e in es):
                return "family"
        if any(lo <= int(e["center_pos"]) <= hi for e in es for lo, hi in sites[sample]):
            return "family-other"
        if known_seqs and unique_seq.get((sample, uid)) in known_seqs:
            return "real"
        return "error"

    known_seqs = set()
    if a.known:
        cur = []
        for line in open(a.known):
            line = line.strip()
            if line.startswith(">"):
                if cur:
                    known_seqs.add("".join(cur))
                cur = []
            else:
                cur.append(line)
        if cur:
            known_seqs.add("".join(cur))
        dereps = [a.derep] if a.derep else [
            os.path.join(a.derep_dir, f) for f in sorted(os.listdir(a.derep_dir))]
        unique_seq = {}
        for path in dereps:
            d = json.load(open(path))
            unique_seq.update({(d["sample"], str(i)): u["sequence"] for i, u in enumerate(d["uniques"])})
    disagree = {k for k, es in by_unique.items() if es[0]["nsubs"] == "0" and float(es[0]["lambda"]) < 1e-3}
    print(f"skipped {len(disagree)} uniques whose events and lambda come from different alignments")
    outcomes, bud_rows = [], []
    ceilings = [None] + [float(c) for c in a.hp_ceiling.split(",") if c]
    for guard, ceiling in [(g, c) for g in map(int, a.guard.split(",")) for c in (ceilings if g else [None])]:
        budded = []
        for (sample, uid), es in by_unique.items():
            if (sample, uid) in disagree:
                continue
            e0 = es[0]
            internal = [e for e in es if a.terminal or e["terminal"] == "0"]
            if a.end_zone:
                internal = [e for e in internal if int(e["center_len"]) - int(e["center_pos"]) > a.end_zone]
            internal, pairs = pair_up(internal)
            guarded = [e for e in internal if guard and involved_run(e) >= guard]
            counted = [e for e in internal if e not in guarded]
            factor = math.prod(p_event(e) for e in counted) * math.prod(p_pair(pp) for pp in pairs)
            counted = counted + [pp[0] for pp in pairs]
            if ceiling is not None and guarded:
                factor *= ceiling ** len(guarded)
                counted = counted + guarded
            if not counted:
                continue
            lam = float(e0["lambda"]) * factor
            if buds(int(e0["reads"]), lam, int(e0["cluster_reads"]), nraw[sample]):
                budded.append(((sample, uid), es, counted, classify(sample, uid, es)))
        cls = collections.Counter(b[3] for b in budded)
        false = [b for b in budded if b[3] == "error"]
        label = f"guard run>={guard}" if guard else "no guard"
        if ceiling is not None:
            label += f", ceiling {ceiling:g}"
        for (s_, u), es, c, k in budded:
            bud_rows.append((label, s_, u, es[0]["cluster"], es[0]["reads"], es[0]["cluster_reads"],
                             es[0]["nsubs"], k, ";".join(f"{e['type']}{e['len']}@{e['center_pos']}:{e['bases'][:12]}"
                                                         f"(hp{e['hp_len']})" for e in c)))
        print(f"\n== {label}{' (+terminal)' if a.terminal else ''}: {len(budded)} buds: "
              f"{cls['variant']} variants, {cls['family']} family, "
              f"{cls['family-other']} family in another cluster, {cls['real']} real alleles, "
              f"{cls['error']} errors")
        for (s_, u), es, c, k in budded:
            if k == "real":
                desc = ", ".join(f"{e['type']}{e['len']}@{e['center_pos']}:{e['bases'][:8]}(hp{e['hp_len']},q{e['qual_min']})" for e in c)
                print(f"    real {s_} u{u} reads={es[0]['reads']} cluster_reads={es[0]['cluster_reads']}  {desc}")
        if false:
            top = sorted(false, key=lambda b: -int(b[1][0]["reads"]))[:5]
            for (s_, u), es, c, _ in top:
                desc = ", ".join(f"{e['type']}{e['len']}@{e['center_pos']}:{e['bases'][:8]}"
                                 f"(hp{e['hp_len']},q{e['qual_min']})" for e in c)
                print(f"    error {s_} u{u} reads={es[0]['reads']} cluster_reads="
                      f"{es[0]['cluster_reads']} nsubs={es[0]['nsubs']}  {desc}")
        if variants:
            hit = {k for k, _, _, c in budded if c == "variant"}
            rows = collections.defaultdict(lambda: [0, 0])
            for t in variants:
                key = "hp-control" if t["context"] == "hp" else f"size {int(t['size']):>3}"
                rows[key][1] += 1
                # Recovered today (an ASV already) or would bud now.
                today = t["sequence"] in asvs[t["sample"]]
                now = (t["sample"], t["unique"]) in hit
                if today or now:
                    rows[key][0] += 1
                outcomes.append((f"{guard}:{ceiling}", t["variant"], int(today), int(today or now)))
            print("  spikes recovered (today or now):",
                  "  ".join(f"{k}:{v[0]}/{v[1]}" for k, v in sorted(rows.items())))
            moved = collections.defaultdict(lambda: [0, 0])
            for t in variants:
                if t["context"] == "plain" and int(t["size"]) <= 10:
                    n = int(t["reads_moved"])
                    b = "1" if n <= 1 else "2-4" if n <= 4 else "5-19" if n <= 19 else "20-99" if n <= 99 else ">=100"
                    moved[b][1] += 1
                    if t["sequence"] in asvs[t["sample"]] or (t["sample"], t["unique"]) in hit:
                        moved[b][0] += 1
            print("  size<=10 by reads moved:", "  ".join(
                f"{b}:{moved[b][0]}/{moved[b][1]}" for b in ["1", "2-4", "5-19", "20-99", ">=100"] if b in moved))
    if a.write_buds:
        with open(a.write_buds, "w") as fh:
            fh.write("setting\tsample\tunique\tcluster\treads\tcluster_reads\tnsubs\tclass\tcounted_events\n")
            for r in bud_rows:
                fh.write("\t".join(map(str, r)) + "\n")
    if a.write_recovered:
        with open(a.write_recovered, "w") as fh:
            fh.write("guard\tvariant\trecovered_today\trecovered\n")
            for o in outcomes:
                fh.write("\t".join(map(str, o)) + "\n")


if __name__ == "__main__":
    main()
