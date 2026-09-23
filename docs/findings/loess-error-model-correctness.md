# The LOESS error model: fidelity, and a silent floor

**The error model is an amplifier.** From the first weeks of the port it was
clear that small deviations in the fitted error rates move the ASVs and the
counts that come out the other end — so "close enough" is not a property the
smoother is allowed to have. Everything on this page follows from that.

Two arcs, opposite in character:

1. **Fidelity (#4, #14).** Chasing the error model to R took three distinct
   fixes and ruled out two plausible hypotheses. It ended with native LOESS
   **bit-equivalent to R's `loess(surface = "direct")` to machine precision**
   (1.97e-13), and with every remaining difference from R DADA2 attributed to
   one thing: the kd-tree interpolation surface R defaults to.
2. **A silent floor (#95, #97).** On binned-quality data the default errfun
   returned an error model pinned at the `1e-7` floor in **every one of 480
   cells** — no warning, output structurally indistinguishable from a real
   model. R handled the same input. It bit exactly the regime sequencing is
   moving toward: NovaSeq bins quality to 4 distinct values, the Illumina i100
   data we have seen carries 3.

The two are connected, and not pleasantly: **the boundary convention that fixed
the first arc is the channel through which the second one failed silently.**

Along the way, three claims had to be withdrawn, and all three were the
instrument's fault rather than the thing being measured. That pattern is the
third result on this page.

## Why fidelity here is not pedantry

Error rates feed `pval.rs` pointwise, so a shifted rate moves abundance
p-values, which moves partition decisions. The effect is small but it is
systematically *there*:

- On F3D0, the three #4 fixes moved steady-state clusters from **194 to 132** —
  R's answer is 132 — and total transition counts from 1,579,172 to 1,700,083
  against R's 1,700,323 (−0.014%).
- On a 363-sample run, three ways of producing the error model ranked
  consistently in the final output: R `learnErrors` ≻ our `learn-errors` with an
  external R loess ≻ our native loess. The worst case was 240 sequences out of
  3.6M and ≤1 ASV per sample — close to noise, but it *ordered* the same way
  every time, which is the signature of a real effect rather than jitter.

So the target is not aesthetic. It is that the choice of smoother should not be
visible in the ASV table.

## Arc one: chasing R (#4, #14)

This is where the project changed what it was optimising for. The original port
aimed at **workflow completeness** — FASTQ in, sequence table out, ASVs and
counts that looked reasonable — with correctness against R not yet the emphasis.
Once that pipeline existed end to end, attention turned to R concordance, and
the first place a real deviation surfaced was the LOESS fit. Small differences
in the fitted rates were producing different ASVs. The habit of scoring changes
on ASV-level concordance rather than on how close a curve looks, which the rest
of these pages take for granted, starts here.

### The low-Q extrapolation fix

Our `loess_predict` returned a polynomial extrapolation at every x. R's
`predict.loess` returns `NA` outside the fitted-data range, and `loessErrfun`
then flat-fills from the nearest fitted value. On F3D0, where nothing is
observed below Q=12, we extrapolated A2C at Q=0 to **0.0083 against R's
0.0637**.

Fixing it — return `None` outside `[x_min, x_max]`, let `extrapolate_flat`
handle the boundary — dropped mean absolute `err_mat` difference per row from
**1–30% to 0.05–0.15%**.

Two other fixes landed alongside it: R's `nconsist=0` init pass, which we
lacked entirely (run `dada` with `MAX_CLUST=1`, fit the errfun on the
accumulated trans, force the diagonals to 1.0), and a Poisson fix that was what
finally moved steady state from 194 clusters to 132.

**Remember the `None`-outside-range convention.** It is correct, it is what R
does, and it is the exact mechanism that made the second arc silent.

### Two hypotheses, falsified

The residual against R was ~1e-3 to 1e-2, concentrated at low-Q edges. The
obvious explanation was conditioning: we solved the local weighted least squares
with **normal equations**, whose condition number goes as cond(X)², while R's
`loessc.c` uses QR. At low-Q edges the local x-range shrinks and the basis
`[1, q, q²]` goes nearly collinear — the story fit the evidence exactly.

It was wrong. Weighted Householder QR with a centred basis came out
**bit-identical**, iteration counts included. Branch deleted. Conditioning was
never the bottleneck, and the mechanism that "explained" the edge signature had
nothing to do with it.

What did explain it was smaller and duller: **`n_local` rounding.** R's docs say
`round(α · n)`; `loessc.c` actually computes `(int)(s · n)`, which is *floor*.
We used `ceil`. At `nv = 35`, span 0.75, that is 26 against our 27.

With that one change, native LOESS matched R `loess(surface = "direct")` at
**1.97e-13** — machine precision on real 362-sample data. It also closed four
other checklist items at a stroke (tricube boundary epsilon, tie-breaking, NaN
bookkeeping, exact-vs-approximate statistics): they would have shown up in that
comparison and did not.

### What was left, and the port that addressed it

| comparison | max abs diff |
|---|---|
| ours ↔ R `surface = "direct"` | **1.97e-13** |
| ours ↔ R `surface = "interpolate"` (R DADA2's default) | ~1e-3 to 1e-2 |

No unexplained divergence remained — all of it was the kd-tree surface R
defaults to. So `ehg124`/`ehg128` were ported (kd-tree partition subdividing at
the *median data point*, plus cubic Hermite blending) and exposed through
`LoessConfig` and `--loess-preset`: `default` keeps the historical direct
surface, `r-dada2` mirrors R DADA2's interpolate surface.

On the 362-sample data the two presets produce **visibly different error models
and an identical set of ASVs** — which is the most reassuring possible outcome
for a knob, and the empirical bound on how much this particular choice matters.

R DADA2's use of `surface = "interpolate"` appears to be R's default rather than
a deliberate choice; on an integer-Q grid, where every data point is already a
vertex, direct evaluation is arguably the more accurate of the two.

### Two ways to run R's error model inside the Rust workflow

The concordance push produced two escape hatches that outlived their diagnostic
purpose and are now supported features, documented together in [Using an
external error model](../walkthroughs/external-error-models.md):

- **`--errfun external --errfun-cmd "<command>"`** hands the fit to any external
  program. `dada2-rs` writes a transition TSV, runs `<command> <trans-tsv>
  <err-tsv>` once per self-consistency iteration, and reads the result back.
  `examples/external_errfun/` ships R references (`loess_reference.R`,
  `loess_reference_direct.R`, `pacbio_reference.R`, a modified `loess_modified.R`)
  and a pure-stdlib Python `noqual.py`. It is the way to prototype a model for a
  new chemistry, or to run a published R errfun, without rebuilding anything.
- **`scripts/learnerrors_to_dada2rs.R`** converts an R DADA2 `learnErrors()`
  `.rds` straight to a dada2-rs JSON error model for `--error-model`, so R's own
  fitted model can drive Rust inference.

Both were built to isolate the fit from everything around it — and both showed
the same thing. **With either R-sourced model, ASVs and counts hew much closer
to R's**, which is what established that the smoother was the thing to chase.
Neither is bit-exact with full R DADA2, though, and the most likely reason is
not the fit at all: **the two runs may be learning from different input
sequences.**

That hypothesis has since been measured, and it is larger than it sounds. A
model's sensitivity to *which samples were drawn* at a fixed `--nbases` budget is
comparable to its sensitivity to the budget itself — see
[`--nbases` and error-model convergence](learn-errors-nbases-convergence.md),
and [#68](https://github.com/HPCBio/dada2-rs/issues/68) for the sample-level
accumulation that causes it. An error-model comparison that does not pin the
input set is measuring two things at once.

### A correction that reshaped the presets

For a while this investigation believed R DADA2 did **not** clamp the fitted
rates, and the `r-dada2` preset was built to match: interpolate surface, no
clamp. That was wrong. R clamps off-diagonals to `[1e-7, 0.25]` in
`errorModels.R:53-56`, in R rather than C, under a comment reading `# HACKY`.

The error came from working off a *summary* of the R source that omitted those
lines, and from three "reference" scripts that were themselves incomplete ports
of `loessErrfun` — so the reference and the implementation shared a blind spot,
and the comparison could not see it. A measured 1.06e-2 hotspot at `A2C q=0..9`,
confidently attributed to the clamp, was a mismeasurement: real R DADA2 produces
0.25 there too.

Corrected, the presets differ in **exactly one thing** — Direct versus
Interpolate. Clamp bounds, `n_local` flooring and blending are shared.

## Arc two: the silent floor (#95, #97)

### A tricube weight of exactly zero

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

**This is where the two arcs meet.** `loess_predict` returns `None` outside the
fitted range because that is what R does, and adopting it is what closed the
low-Q gap in arc one. It made `None` an *expected* value with a defined
fallback — so when a genuinely unfittable neighbourhood later produced the same
`None`, the pipeline handled it quietly and correctly-looking, instead of
failing. A correct fix built the channel the later fault travelled down. Worth
remembering when adding a sentinel: every "this means fall back" path is also a
path a real failure can hide in.

Two things kept it hidden:

- **The output looked fine.** A uniform matrix is a valid matrix. Nothing
  downstream could tell it from a learned model, and no message was emitted.
- **The obvious path avoided it.** `binned_qual_errfun` interpolates between
  anchors and never calls LOESS, so only users running the *default* errfun on
  binned input were affected — a plausible thing to do, and silent when done.

R's `loess(rlogp ~ q, weights = tot)` returns finite, sensible values at all 40
columns on the identical input, so this was ours alone.

### The fix

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

## Three withdrawn claims, and the pattern they share

The `loess-rs` write-up contained two errors, both flattering to us and both the
instrument's fault rather than the crate's:

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

Add the R-clamp mistake from arc one and there are three, with one shape between
them. In each case the *instrument* was wrong — a guard in our own probe, a
missing API scored as inaccuracy, a summary of source standing in for the source
— and in each case the error flattered us or simplified the story. The rule:
**when a measurement says the other implementation is wrong, suspect the
measurement first.** The same discipline this project already applies to
a null result applies to a flattering one. Rebuilding the probe as a committed
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

- **The error model is an amplifier, so fidelity work on it is justified by
  default.** Small rate changes move partition decisions; the #4 fixes moved
  F3D0 from 194 clusters to R's 132. The bar is that the choice of smoother
  should not be visible in the ASV table — which is also the *only* acceptable
  score for a proposed change, not how the curves look.
- **Binned quality is the stress case for the error model, not an edge case.**
  Every defect on this page was invisible on dense MiSeq data and obvious on 3–4
  distinct quality values. New error-model work gets a sparse arm from the start.
- **Keep the smoother in-house.** 455 zero-dependency lines that are closer to
  R on the surface we ship beat a pre-1.0 dependency that cannot predict at an
  unpopulated column. Revisit only against the three conditions above.
- **A fit that cannot be made must be reported, never floored.** The `Result`
  return is the durable part of #97; the degree fallback merely makes the error
  rare.
- **Two supported routes exist for running R's error model** — `--errfun
  external` for the fit, `scripts/learnerrors_to_dada2rs.R` for a finished
  `learnErrors()` model. They are the answer for anyone who needs R parity today,
  and the reference arm for any future comparison.
- **Pin the input set before comparing error models.** Otherwise the fit and the
  sample draw move together, and the residual cannot be attributed.
- **`--loess-preset` is the fidelity knob**, and the two presets differ in
  exactly one thing: the fitting surface. `default` is direct, `r-dada2` is R's
  interpolate. On 362 samples they give different error models and the same
  ASVs.
- **Do not tune the smoother on curve shape.** Every arm in #96 changes rates in
  the high-count region; the only acceptable score is ASV-level churn.
- One user-facing note worth carrying upstream: R's `Error rates could not be
  estimated (this is usually because of very few reads)` misattributes this
  failure. It is `tryCatch` catching an all-`NA` prediction, and on binned data
  the cause is too few *distinct quality values*, not too few reads.
