# The LOESS error model: a silent floor, and the probe that found it

**On binned-quality data, the default error function returned an error model
pinned at the `1e-7` floor in every cell — with no warning, no error, and output
that was structurally indistinguishable from a real model.** Fixed in
[#97](https://github.com/HPCBio/dada2-rs/pull/97). It bit exactly the regime
sequencing is moving toward: NovaSeq bins quality to 4 distinct values, and the
Illumina i100 data we have seen carries 3.

Three results, and the order they happened in matters:

1. **The bug was found by a probe aimed at something else.** We were evaluating
   whether to adopt an external LOESS crate. The comparison harness turned up a
   defect in *our* code, not theirs.
2. **Two of that evaluation's headline claims were wrong, and both were the
   instrument's fault.** "The crate fails at n=3" was a guard in our own throwaway
   probe; "it is catastrophic on binned input" was a missing API being scored as
   an accuracy failure. Corrected in public, and the harness is now committed so
   the numbers do not have to be re-derived from memory.
3. **The recommendation did not change** — do not adopt — but for different
   reasons than first given, and the crate is better than the first write-up said.

## The bug: a tricube weight of exactly zero

`fit_local_at` set `max_dist` to the distance of the farthest *included*
neighbour. That point therefore gets `u == 1.0`, a tricube weight of exactly
`0.0`, and is then dropped by the `w > 0.0` filter. **A neighbourhood of
`n_local` points yields at most `n_local - 1` usable observations.**

With the DADA2 defaults of `span = 0.75` and `degree = 2` (so `p = 3`
coefficients):

| populated quality columns | `n_local` | surviving | result |
|---|---|---|---|
| 5 | 3 | 2 | `2 < p` → `None` |
| 6 | 4 | 3 | `3 == p` → fits |

`None` became `NaN` in `extrapolate_flat`, and `NaN` became `min_error_rate` in
`loess_errfun`. Every cell, all 480 of them, at `1e-7`.

Two things kept it hidden:

- **The output looked fine.** A uniform matrix is a valid matrix. Nothing
  downstream could tell it from a learned model, and no message was emitted.
- **The obvious path avoided it.** `binned_qual_errfun` interpolates between
  anchors and never calls LOESS, so only users running the *default* errfun on
  binned input were affected — a plausible thing to do, and silent when done.

R's `loess(rlogp ~ q, weights = tot)` returns finite, sensible values at all 40
columns on the identical input, so this was ours alone.

## The fix

Two changes, deliberately scoped:

- **Fit the highest degree the surviving neighbourhood can identify**, bounded
  by both the positively-weighted observation count and the number of *distinct*
  abscissae — a degree-`d` polynomial needs `d+1` distinct x values however many
  points it has — retrying at lower degree if the solve is numerically singular.
  This engages **only** in the degenerate regime, so dense fits stay
  byte-identical and the [concordance guardrail](../benchmarking.md#5-concordance-validation-tooling)
  is unaffected. It degrades all the way down to a single populated column,
  which fits a local constant: quality-independent, but an honest estimate of
  the observed rate.
- **`loess_errfun` and `pacbio_errfun` now return `Result`**, matching
  `binned_qual_errfun`, and report the populated-column count rather than
  emitting a structurally-valid but unusable matrix. Because the degree fallback
  handles everything down to one column, the guard only fires on a transition
  matrix with *no* observations at all — which also used to floor silently.

The general point: **a numerical routine that cannot fit should say so.**
Returning a plausible-looking constant is worse than returning an error, because
it consumes the one signal a user has that something went wrong.

## Should we adopt `loess-rs`? No — and here is the honest table

Evaluated at 0.9.0 against R 4.6.0 `stats::loess` and dada2 1.40.0, on the F3D0
MiSeqSOP transition matrix (12 transitions × 40 quality bins, 480 cells).
Everything is measured as `|log10(rate_arm / rate_R)|` **after** the `10^` and
the `[1e-7, 0.25]` clamp — the number that actually reaches inference, not the
raw fit residual.

| arm | surface | median | p95 | max |
|---|---|---|---|---|
| `loess-rs` | `direct` | **2.3e-15** | 7.1e-15 | 2.5e-14 |
| ours | `direct` | 4.4e-14 | 1.0e-12 | 7.7e-12 |
| `loess-rs` | `interpolate` | 5.1e-15 | 2.3e-01 | **2.4e-01** |
| ours | `interpolate` | 1.0e-12 | 7.3e-04 | 1.2e-03 |

- **Their `custom_weights` is exactly R's `weights=`** — confirmed in source, and
  empirically cleaner than our own port on the `direct` surface. That feature
  exists because we asked for it upstream.
- **Their interpolation surface is not R's `ehg128`.** It matches to machine
  epsilon across the interior (q=15–36) and diverges at the boundary vertices by
  up to 0.24 log10, about 1.7× in rate. Ours is ~200× closer there — **and the
  interpolate surface is what the `r-dada2` preset uses**, so on the path we care
  about most, ours wins.

In short: their better arm is one we do not use, and the worse arm is the one we
do. The rest of the decision is unglamorous and decisive — perf is irrelevant at
n=41 points and 12 fits per round; `nalgebra` + `wide` + `num-traits` against our
455 zero-dependency lines matters given the crates.io surgery the WFA git
dependency already forces ([#63](https://github.com/HPCBio/dada2-rs/issues/63));
0.2.2 → 0.9.0 in six months is API churn our concordance guardrail makes
expensive; and there is **no public predict-at-new-x**, so it cannot produce a
value at an unpopulated quality column at all.

**Decision: parked, not adopted.** What would reopen it: `ehg128` parity at the
boundaries, a published R-reference suite with stated tolerances, and a 1.0.

## Two withdrawn claims, and what they cost

The first write-up of that evaluation contained two errors, both flattering to
us and both the instrument's fault rather than the crate's:

1. **"`loess-rs` fails at n_valid = 3."** It does not; it matches R to ~4e-16
   there. The failure was an `if idx.len() < 4 { return NaN }` guard in our own
   throwaway probe, misread as the crate declining the fit. So `loess-rs` was
   better than us across the *whole* sparse range before #97, not merely at some
   of it.
2. **"It is catastrophic on binned input."** That was a **coverage gap scored as
   an accuracy failure**: with no public predict-at-new-x it cannot emit a value
   at an unpopulated column, and the harness was writing those cells out as the
   `1e-7` floor. On the production-shaped binned case it matches R to 3.4e-16 on
   the 192 cells it *can* cover, and is silent on the other 288. Still a genuine
   blocker for us — `loess_errfun` must emit a rate for every column — but a
   missing-API criticism, not an accuracy one.

The rule both share: **when your instrument reports that someone else's tool
failed, suspect the instrument first.** The same discipline that applies to a
null result applies to a flattering one. Rebuilding the probe as a committed
harness — `dev/loess-oracle`, `./run.sh <learn_errors.json>` — is what exposed
both, and is why this is not a story that has to be re-derived next time.

### One result that goes against us

The corrected sparsity sweep found a place where `loess-rs` is genuinely better,
on the surface we actually ship. Max `|log10|` vs R, interpolate, anchors only:

| populated columns | ours | `loess-rs` |
|---|---|---|
| 3 | 2.5e-01 | 1.5e-01 |
| 4 | **7.3e-01** | 2.8e-04 |
| 5 | **5.4e-01** | 2.8e-03 |
| 6 | **7.1e-01** | 1.9e-02 |
| 8 | 2.1e-01 | 4.5e-02 |

Our kd-tree vertex partition has very little to work with at 4–6 anchors. No
product impact today — nothing ships a sparse-input interpolate path, and R
itself is unreliable in this regime — but it is a real gap on the `r-dada2`
surface and it belongs in any serious binned-quality work.

## Still open: robustness, weighting and the boundary

Tracked in [#96](https://github.com/HPCBio/dada2-rs/issues/96). These are not
cosmetic — they move rates in the region that decides calls:

| arm | median `|log10|` | max | cells moving >2× |
|---|---|---|---|
| bisquare + MAD, 4 iterations | 3.5e-02 | 4.1e-01 | 8/480 |
| boundary = reflect | 4.4e-02 | 3.2e-01 | 13/480 |
| boundary = extend | 2.3e-02 | 2.5e-01 | 0/480 |

Robustness iterations systematically **raise** error rates at q=31–37 (+0.07 to
+0.09 log10 mean signed) — the high-count region that dominates inference — so
adopting them would make calling more conservative. Whether that is an
improvement is an empirical question, and `lrs_robust_mad` is itself unstable at
small n (5.0 log10 divergence at 7 anchors, 4/84 cells floored): **robustness
iterations need enough points to estimate a scale from**, which is precisely what
binned data does not have.

### Prior art: R got here first

benjjneb/dada2#1307 collects community variants, of which `loessErrfun_mod4` is
reported most robust in practice: `degree = 1`, `span = 0.95`,
`weights = log10(tot)`, plus a post-hoc floor at the max-Q column. Four deltas,
each aimed at something we have also hit:

- **`degree = 1`** solves the degeneracy of #95 globally, where our fallback
  reduces degree only where the neighbourhood cannot support a quadratic. Ours
  keeps dense fits byte-identical; mod4 perturbs *every* fit and so is not
  R-concordant by construction. Both philosophies are worth an A/B.
- **`weights = log10(tot)`** attacks leverage by compressing the weight range
  rather than downweighting outliers — on F3D0, `tot` spans 75–94,451 (~1250×)
  and `log10(tot)` compresses that to 1.9–5.0 (~2.6×). A blunter instrument than
  bisquare + MAD, at the same target.
- **The max-Q floor** targets real oscillation, but it is not monotonicity, and
  it hardcodes the `X40` column — so it does not transfer to PacBio Q93 or any
  other maximum. A portable version would reference the last populated column.

mod4's robustness is reported from practice, not from a controlled evaluation
against a truth set. So the arm list is: stock, bisquare + MAD, `log10(tot)`
weights, mod4's full parameter set, and a portable floor — **scored on ASV-level
set identity, not on how the curves look.**

## What this dictates

- **Binned quality is the stress case for the error model, not an edge case.**
  Every defect on this page was invisible on dense MiSeq data and obvious on 3–4
  distinct quality values. New error-model work gets a sparse arm from the start.
- **Keep the smoother in-house.** 455 zero-dependency lines that are closer to
  R on the surface we ship beat a pre-1.0 dependency that cannot predict at an
  unpopulated column. Revisit only against the three conditions above.
- **A fit that cannot be made must be reported, never floored.** The `Result`
  return is the durable part of #97; the degree fallback merely makes the error
  rare.
- **Do not tune the smoother on curve shape.** Every arm in #96 changes rates in
  the high-count region; the only acceptable score is ASV-level churn.
- One user-facing note worth carrying upstream: R's `Error rates could not be
  estimated (this is usually because of very few reads)` misattributes this
  failure. It is `tryCatch` catching an all-`NA` prediction, and on binned data
  the cause is too few *distinct quality values*, not too few reads.
