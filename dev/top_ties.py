#!/usr/bin/env python3
"""Flag derep inputs whose top uniques tie on count: the only case #239 changes.

A tie there is the one input on which `assign_center`'s tie rule matters. None
found means #239 cannot change that run's output, so an A/B is unnecessary.

Usage: top_ties.py [--pooled] derep.json [...]
  default   one line per sample (per-sample dada, pseudo, learn-errors)
  --pooled  also sums counts across the given dereps, as dada-pooled's merge does
"""
import json, sys
from collections import Counter

args = sys.argv[1:]
pooled = "--pooled" in args
paths = [a for a in args if a != "--pooled"]

def top(counts):
    c = sorted(counts, reverse=True)
    return c[0], (c[1] if len(c) > 1 else None), sum(1 for x in c if x == c[0])

pool, ties = Counter(), 0
for p in paths:
    u = json.load(open(p))["uniques"]
    t1, t2, n = top(x["count"] for x in u)
    ties += n > 1
    print(f"{p.split('/')[-1]:<32} n_uniq={len(u):>7} top={t1:>7} second={t2!s:>7}{'  TIE' if n > 1 else ''}")
    if pooled:
        for x in u:
            pool[x["sequence"]] += x["count"]
print(f"per-sample ties: {ties} of {len(paths)}")
if pooled:
    t1, t2, n = top(pool.values())
    print(f"POOLED n_uniq={len(pool)} top={t1} second={t2} tied_at_top={n}{'  TIE' if n > 1 else ''}")
