#!/usr/bin/env python3
"""indel_fit_proto.py — compare candidate indel-rate fits out of sample (issue #291).

Fits each candidate on one set of `--indel-opportunities` tables and scores it
on another by Poisson deviance per held-out event (lower is better), so a fit
that chases noise in sparse cells is penalised. Also reports calibration
(predicted / observed events by quality band): under-prediction is the
dangerous direction, since it splits errors off as variants.

Tables are summed over the chosen samples (`--train-samples`, `--test-samples`;
default all rows). For a `learn-errors --indel-history` table, select the final
iteration, e.g. `--train-samples iter3`.

Candidates
  insertions (ins_q_novel, ins_q_run; rate by inserted-base quality):
    pool        per-cell rate, widening a window until it holds 20 events (current)
    loess       weighted local-linear fit of log10((events+1)/opportunities) on Q,
                span 0.75, weights = opportunities (DADA2's loessErrfun recipe);
                the top quality is taken from its own rate when --cap-q is set
    glm_lin     Poisson GLM, log rate linear in Q
    glm_spline  Poisson GLM, natural spline in Q (df 4)
  deletions (del_fq_L1..L5; rate by run length class and lower flanking quality):
    product     rate(run length) x enrichment(flanking Q), from the two marginals (current)
    additive    Poisson GLM: run length + spline(flanking Q)
    interaction Poisson GLM: run length x spline(flanking Q)
    pool_joint  per-cell pooling within each run length
  base identity (del_hp_X, ins_hp_X; rate by run length):
    L, L+base, L*base  Poisson GLMs

Usage:
  indel_fit_proto.py --train A.tsv --test B.tsv [--train-samples iter3] [--test-samples iter3] [--cap-q 93]
"""

import argparse
import collections
import csv
import math
import warnings

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

warnings.filterwarnings("ignore")


def load(path, samples):
    keep = set(samples.split(",")) if samples else None
    rows = collections.defaultdict(lambda: [0, 0])
    with open(path) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            if keep and r["sample"] not in keep:
                continue
            c = rows[(r["table"], int(r["key"]))]
            c[0] += int(r["opportunities"])
            c[1] += int(r["events"])
    if not rows:
        raise SystemExit(f"{path}: no rows for samples {samples}")
    df = pd.DataFrame([(t, k, o, e) for (t, k), (o, e) in rows.items()], columns=["table", "key", "opp", "ev"])
    return df


def deviance(obs, pred):
    """Poisson deviance summed over cells; pred are expected counts."""
    pred = np.maximum(pred, 1e-300)
    term = np.where(obs > 0, obs * np.log(obs / pred), 0.0) - (obs - pred)
    return 2 * term.sum()


# --- insertion fits: return a rate function key -> rate --------------------------------

def fit_pool(tr):
    tr = tr.sort_values("key")
    keys, o, e = tr.key.to_numpy(), tr.opp.to_numpy(), tr.ev.to_numpy()
    rates = {}
    for i, k in enumerate(keys):
        w = 0
        while True:
            lo, hi = max(0, i - w), min(len(keys) - 1, i + w)
            if e[lo:hi + 1].sum() >= 20 or (lo == 0 and hi == len(keys) - 1):
                rates[k] = e[lo:hi + 1].sum() / max(o[lo:hi + 1].sum(), 1)
                break
            w += 1
    ks = np.array(sorted(rates))
    return lambda q: np.array([rates[ks[np.abs(ks - x).argmin()]] for x in q])


def wloess(x, y, w, at, span=0.75):
    """Weighted local-linear regression with tricube kernel (R loess degree 1)."""
    n = len(x)
    k = max(int(math.ceil(span * n)), 2)
    out = []
    for x0 in at:
        d = np.abs(x - x0)
        h = np.sort(d)[k - 1] or 1e-9
        u = np.clip(d / (h * 1.0000001), 0, 1)
        kw = (1 - u ** 3) ** 3 * w
        X = np.column_stack([np.ones(n), x - x0])
        W = np.diag(kw)
        try:
            beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
            out.append(beta[0])
        except np.linalg.LinAlgError:
            out.append(np.average(y, weights=kw) if kw.sum() > 0 else y.mean())
    return np.array(out)


def fit_loess(tr, cap_q):
    tr = tr[tr.opp > 0]
    body = tr[tr.key != cap_q] if cap_q is not None else tr
    x, o, e = body.key.to_numpy(float), body.opp.to_numpy(float), body.ev.to_numpy(float)
    y = np.log10((e + 1) / o)
    cap_rate = None
    if cap_q is not None and (tr.key == cap_q).any():
        c = tr[tr.key == cap_q].iloc[0]
        cap_rate = (c.ev + 1) / c.opp

    def f(q):
        q = np.asarray(q, float)
        r = 10 ** wloess(x, y, o, np.clip(q, x.min(), x.max()))
        if cap_rate is not None:
            r = np.where(q == cap_q, cap_rate, r)
        return r
    return f


def fit_glm(tr, formula_rhs, cap_q=None):
    tr = tr[tr.opp > 0].copy()
    tr["cap"] = (tr.key == cap_q).astype(float) if cap_q is not None else 0.0
    rhs = formula_rhs + (" + cap" if cap_q is not None and tr.cap.sum() > 0 else "")
    y, X = patsy.dmatrices(f"ev ~ {rhs}", tr, return_type="dataframe")
    m = sm.GLM(y, X, family=sm.families.Poisson(), offset=np.log(tr.opp)).fit()
    design = X.design_info

    def f(df):
        df = df.copy()
        df["cap"] = (df.key == cap_q).astype(float) if cap_q is not None else 0.0
        (Xn,) = patsy.build_design_matrices([design], df, return_type="dataframe")
        return np.exp(Xn @ m.params).to_numpy()
    return f, len(m.params)


def score_ins(train, test, cap_q):
    print("== insertions: held-out deviance per held-out event (lower is better), calibration pred/obs by Q band")
    bands = [(0, 19), (20, 29), (30, 39), (40, 59), (60, 79), (80, 92), (93, 93)]
    for table in ("ins_q_novel", "ins_q_run"):
        tr, te = train[train.table == table], test[test.table == table]
        te = te[te.opp > 0]
        kmin, kmax = tr.key.min(), tr.key.max()
        te = te[(te.key >= kmin) & (te.key <= kmax)]
        n_ev = te.ev.sum()
        fits = {}
        fits["pool"] = (fit_pool(tr), "cells")
        fits["loess"] = (fit_loess(tr, cap_q), "smooth")
        g, k = fit_glm(tr, "key", cap_q)
        fits["glm_lin"] = (lambda q, g=g: g(pd.DataFrame({"key": q})), k)
        g, k = fit_glm(tr, f"cr(key, df=4)", cap_q)
        fits["glm_spline"] = (lambda q, g=g: g(pd.DataFrame({"key": q})), k)
        print(f"   {table} (test events {n_ev}):")
        for name, (f, k) in fits.items():
            pred = te.opp.to_numpy() * f(te.key.to_numpy())
            d = deviance(te.ev.to_numpy(), pred) / max(n_ev, 1)
            cal = []
            for lo, hi in bands:
                sel = (te.key >= lo) & (te.key <= hi)
                if te.ev[sel].sum() >= 5:
                    cal.append(f"Q{lo}-{hi}:{pred[sel.to_numpy()].sum() / te.ev[sel].sum():.2f}")
            print(f"      {name:11s} params={k!s:6s} dev/event={d:7.3f}   " + " ".join(cal))


def joint_del(df):
    x = df[df.table.str.startswith("del_fq_L")].copy()
    x["L"] = x.table.str[-1].astype(int)
    return x[x.opp > 0]


def score_del(train, test, cap_q):
    print("== deletions (run length x flanking Q): held-out deviance per held-out event, calibration pred/obs by run length")
    tr, te = joint_del(train), joint_del(test)
    te = te[(te.key >= tr.key.min()) & (te.key <= tr.key.max())]
    n_ev = te.ev.sum()
    preds = {}
    # product of marginals (current model)
    rl = tr.groupby("L")[["opp", "ev"]].sum()
    rate_L = rl.ev / rl.opp
    fq = tr.groupby("key")[["opp", "ev"]].sum().sort_index()
    mean = fq.ev.sum() / fq.opp.sum()
    enrich = fit_pool(fq.reset_index().assign(table="x"))
    preds["product"] = (te.opp.to_numpy() * rate_L.reindex(te.L).to_numpy() * enrich(te.key.to_numpy()) / mean, "marginals")
    for name, rhs in (("additive", "C(L) + cr(key, df=4)"), ("interaction", "C(L) * cr(key, df=3)")):
        g, k = fit_glm(tr, rhs, cap_q)
        preds[name] = (te.opp.to_numpy() * g(te), k)
    pj = np.empty(len(te))
    for L in te.L.unique():
        f = fit_pool(tr[tr.L == L])
        sel = (te.L == L).to_numpy()
        pj[sel] = f(te.key.to_numpy()[sel])
    preds["pool_joint"] = (te.opp.to_numpy() * pj, "cells")
    print(f"   test events {n_ev}")
    for name, (pred, k) in preds.items():
        d = deviance(te.ev.to_numpy(), pred) / max(n_ev, 1)
        cal = " ".join(f"L{L}:{pred[(te.L == L).to_numpy()].sum() / max(te.ev[te.L == L].sum(), 1):.2f}"
                       for L in sorted(te.L.unique()) if te.ev[te.L == L].sum() >= 5)
        print(f"      {name:11s} params={k!s:9s} dev/event={d:7.3f}   {cal}")


def score_base(train, test):
    print("== base identity (rate by run length): held-out deviance per held-out event")
    for kind in ("del_hp", "ins_hp"):
        def long(df):
            x = df[df.table.str.match(f"{kind}_[ACGT]$")].copy()
            x["base"] = x.table.str[-1]
            x["L"] = x.key.clip(upper=5)
            x = x.groupby(["base", "L"], as_index=False)[["opp", "ev"]].sum()
            return x[x.opp > 0]
        tr, te = long(train), long(test)
        n_ev = te.ev.sum()
        out = []
        for name, rhs in (("L", "C(L)"), ("L+base", "C(L) + C(base)"), ("L*base", "C(L) * C(base)")):
            g, k = fit_glm(tr.assign(key=tr.L), rhs)
            pred = te.opp.to_numpy() * g(te.assign(key=te.L))
            out.append(f"{name}: {deviance(te.ev.to_numpy(), pred) / max(n_ev, 1):.3f} ({k} params)")
        print(f"   {kind} (test events {n_ev}): " + "   ".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--train-samples")
    ap.add_argument("--test-samples")
    ap.add_argument("--cap-q", type=int, help="capped top quality fitted on its own (PacBio: 93)")
    a = ap.parse_args()
    train, test = load(a.train, a.train_samples), load(a.test, a.test_samples)
    score_ins(train, test, a.cap_q)
    score_del(train, test, a.cap_q)
    score_base(train, test)


if __name__ == "__main__":
    main()
