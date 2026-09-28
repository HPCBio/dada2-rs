# LOESS: what is ported from R, and what is not

**Everything R DADA2's `loessErrfun` actually exercises is ported, and it
matches R's `stats::loess` to double-precision round-off.** On the 29-column
MiSeq grid, and on every anchors-only sparse case from 3 to 8 populated quality
columns, the largest per-cell difference from R is under 1e-14 in log10 rate.
With the training set pinned, `learn-errors` reproduces R's whole
self-consistency loop bit for bit: 0 of 656 transition cells differ in either
read direction, and `err_out` agrees to 2e-11.

What is not ported falls into three groups:

1. **Parts of `loess()` that `loessErrfun` never calls** — robustness
   iterations, fit statistics, standard errors, multi-predictor options. They
   cannot affect a DADA2 error model.
2. **Parameters that are fixed rather than exposed** — `span = 0.75` and
   `degree = 2`, hardcoded in R's `loessErrfun` too.
3. **Two degenerate corners where R's own output carries no information** —
   predictions *between* three anchors, and the `direct` surface at the exact
   midpoint of two anchors. We differ from R in both, deliberately. Binned data
   has its own errfun and does not need to go through LOESS.

This page is the ledger. For how the port got here, and why the error model is
worth this much attention, see [The LOESS error model: fidelity, and a silent
floor](loess-error-model-correctness.md).

## The errfun layer

Each R DADA2 error function has a native counterpart.

| R DADA2 | dada2-rs | status |
|---|---|---|
| `loessErrfun` | `--errfun loess` (default) | **Exact.** Same `log10((errs+1)/tot)` response, `weights = tot`, NA on `tot == 0`, flat fill outside the fitted range, `[1e-7, 0.25]` clamp, diagonals as `1 − Σ off-diagonals`. |
| `PacBioErrfun` | `--errfun pacbio` | **Ported.** LOESS on Q0–92 plus `(count + 1) / (total + 4)` at Q93, unclamped as in R. Verified at ASV level (PacBio concordance fixture, 93 = 93 ASVs); no per-fit oracle run yet. Two cosmetic gaps: R prints a message when Q93 is absent and ours falls back silently, and R stops if Q93 is present but not the last column, which our contiguous quality columns cannot produce. |
| `noqualErrfun` | `--errfun noqual` | **One deliberate difference.** Ours clamps to `[1e-7, 0.25]` and R does not. Measured inert: real rates sit three orders inside the window. No LOESS on either side. |
| `makeBinnedQualErrfun` | `--errfun binned-qual` | **Ported, not yet verified against R.** Same piecewise-linear interpolation, same stop on out-of-range data, same warning when an observed extreme misses an anchor, same clamp. The tight-tolerance comparison against R's version is open in [#120](https://github.com/HPCBio/dada2-rs/issues/120). |
| any user errfun | `--errfun external --errfun-cmd …` | **Wire-compatible.** Runs an R or Python errfun inside our self-consistency loop. R's own `loessErrfun` run this way matches R end-to-end to 4e-16. See [Using an external error model](../walkthroughs/external-error-models.md). |

`learnErrors()` output can also be converted and used directly with
[`learnerrors_to_dada2rs.R`](../scripts/pipeline-helpers.md#learnerrors_to_dada2rsr).

## The smoother: `loess()` arguments

`loessErrfun` calls `loess(rlogp ~ q, df, weights = tot)` and leaves every
other argument at its default. The table covers each argument and states
whether it applies to that call.

| `loess()` / `loess.control()` argument | R default | dada2-rs | notes |
|---|---|---|---|
| `formula` | — | one predictor (`q`) | Multi-predictor fitting is not ported; the error model has one predictor. |
| `weights` | none | **ported** | Kernel weight × prior weight, as R. |
| `na.action` | `na.omit` | **equivalent** | Non-finite responses and zero weights are dropped; with `weights = tot` these are the same rows. |
| `span` | 0.75 | **fixed at 0.75**, not exposed | `nf = min(n, floor(n·span))` exactly as R, with no floor at `degree + 1`. `span > 1` bandwidth inflation is ported but unreachable. |
| `enp.target` | — | not ported | An alternative way of setting `span`. |
| `degree` | 2 | **fixed at 2**, not exposed | The kernel takes any degree, but only the errfun calls it. Degree drops automatically where the neighbourhood cannot identify a quadratic (see below). |
| `family` / `iterations` / `iterTrace` | `gaussian`, 4 | **not ported** | Robustness iterations are used only with `family = "symmetric"`, which `loessErrfun` does not request. Whether to adopt them anyway is [#96](https://github.com/HPCBio/dada2-rs/issues/96). |
| `surface` | `interpolate` | **ported, both values** | `--loess-surface`, default `interpolate` since [#205](https://github.com/HPCBio/dada2-rs/issues/205). |
| `cell` | 0.2 | **ported** | `--loess-cell`. Leaf threshold `floor(cell · span · n)`, as R's `fc`. |
| `statistics` / `trace.hat` | `approximate` | not ported | Only affect fit summaries (`enp`, residual SE, trace of the hat matrix), not fitted values. |
| `normalize` / `parametric` / `drop.square` | — | not applicable | Multi-predictor only. |
| `predict(..., se = TRUE)` | `FALSE` | not ported | `loessErrfun` does not request standard errors. |

## The smoother: R's Fortran kernel

The things that make the `interpolate` surface agree with R to round-off rather
than to about 1e-3 all live in `loessf.f`. Each one is cited on the Rust
function that ports it (`src/loess.rs`).

| R routine | behaviour | status |
|---|---|---|
| `ehg126` | kd-tree bounding box padded by 0.5% of the data range per side | **ported** ([#213](https://github.com/HPCBio/dada2-rs/pull/213)) |
| `ehg124` | cut at the median data point; a cell is a leaf when `count ≤ fc` **or** the cut lands on its own boundary | **ported** ([#217](https://github.com/HPCBio/dada2-rs/pull/217)). Vertex sets equal R's `kd$xi` in both the dense and the 3-anchor regime, pinned by tests. |
| `ehg124` / `ehg131` | second leaf test, `diam ≤ fd` | **inert in R**: `lowesd` sets `fd = 0`. Not ported. |
| `ehg129` | tie-shifting for equal abscissae | **unreachable**: the fitted abscissae are distinct quality values. Not ported. |
| `ehg127` | local design centred at the query point; bandwidth = distance to the `nf`-th neighbour × `sqrt(max(1, span))`; tricube kernel | **ported** ([#217](https://github.com/HPCBio/dada2-rs/pull/217)). Centring alone moved agreement from 8.6e-12 to 1.6e-14. |
| `ehg127` | weighted least squares by QR (`dqrdc`/`dqrsl`), condition-number check, SVD pseudoinverse fallback | **different method, same answer**: we solve the normal equations and drop one degree on a singular solve. After centring, the two differ in how a degenerate fit is *detected*, not in the fitted values. |
| `ehg128` | cubic Hermite blending of vertex values and slopes | **ported** |
| `predict.loess` | `NA` outside the data range on `interpolate` | **ported**, then flat-filled as `loessErrfun` does |
| `predict.loess` | polynomial **extrapolation** outside the data range on `direct` | **deliberately not ported.** Ours returns no prediction and flat-fills, so `--loess-surface direct` behaves the way `loessErrfun` would if R's direct surface returned `NA` there. `loess_reference_direct.R` and the oracle both apply the same masking. |

## Measured agreement

`dev/loess-oracle` fits the same `(q, rlogp, tot)` triples with our real
`loess_predict` and with R `stats::loess`, and scores `|log10(rate / rate_R)|`
after the clamp. This is the number that reaches inference. Data: the pinned
362-sample MiSeq SOP forward model.

| case | `direct` vs R direct | `interpolate` vs R interpolate |
|---|---|---|
| full grid (29 populated columns) | 1.6e-14 | 7.5e-15 |
| anchors only, 3 to 8 columns | ≤ 8e-15 | ≤ 8e-15 |
| binned shape, 3 anchors, **scored between anchors** | see below | 2.4e-01 |

End to end, on the same 362-sample run, pooled, with both tools trained on one
pinned 117-sample, 3.01e8-base manifest
([#205](https://github.com/HPCBio/dada2-rs/issues/205)):

| comparison | result |
|---|---|
| forward / reverse `trans`, ours vs R | **0 / 656 cells differ**, both directions |
| self-consistency iterations | 6 / 7, equal |
| `err_out`, ours vs R | max 1.8e-11 / 2.1e-11 |
| our interpolate vs R's interpolate, inside our loop | median 2.4e-14, max 1.8e-11 |
| ASVs | 731 vs 730. The extra ASV comes from a [saturated-birth tie-break](r-parity-floor-and-ceiling.md) and is independent of LOESS. |
| reads | +14 in 2.8M |

## Where we differ from R, deliberately

Both differences are in regimes where R warns and its answer carries no
information. Neither applies to dense Illumina or HiFi quality grids.

**Between three anchors, `interpolate`.** With three populated columns, R
reports `pseudoinverse used` and `reciprocal condition number 0`, and its
between-anchor values come from a rank-deficient solve. We match R exactly *at*
the anchors. Between them we differ by up to 0.24 log10. Matching R here would
mean reproducing a degenerate pseudoinverse; we have not tried to. Binned data
should use `--errfun binned-qual`, which interpolates between the observed bins
directly.

**Midpoint of two anchors, `direct`.** A query point exactly halfway between
two of three anchors is equidistant from both, so the tricube zeroes every
neighbour. R warns `all weights zero`, returns 0, and a log10 rate of 0 becomes
a 100% error rate clamped to 0.25. We return no fit, and the cell takes the
`1e-7` floor: 24 of 492 cells on the binned-shaped oracle case. Neither answer
is informative; ours errs low, R's errs high. `interpolate` is unaffected.

**One populated column.** `interpolate` cannot fit and fails loudly, naming
`--errfun binned-qual`; R fails here too. `direct` fits a local constant.

**No populated columns.** An all-empty transition matrix is an error, not a
uniform `1e-7` matrix. R's equivalent is its `Error rates could not be
estimated` message.

## Not ported, and whether it should be

| item | reason | revisit when |
|---|---|---|
| Robustness iterations (`family = "symmetric"`) | Not used by `loessErrfun`; moves high-Q rates by up to 0.4 log10 | [#96](https://github.com/HPCBio/dada2-rs/issues/96): A/B scored on ASV churn, not curve shape |
| `span` / `degree` as CLI flags | Fixed in R's `loessErrfun`; exposing them invites tuning on curve shape | [#98](https://github.com/HPCBio/dada2-rs/issues/98): a native span=2 binned errfun (dada2#938) needs `span > 1`, already ported |
| Community errfuns (`loessErrfun_mod4`, span=2, …) | Not part of R DADA2 | Run today through `--errfun external`; see [binned-quality prior art](loess-error-model-correctness.md#prior-art-r-got-here-first) |
| Fit statistics, `se = TRUE` | Do not affect fitted values | a use for uncertainty in the error model |
| Multi-predictor LOESS | The error model has one predictor | — |

## Re-checking any of this

```sh
dev/loess-oracle/run.sh <learn_errors.json> out/   # per-fit vs R, ~1 minute
```

Any native `learn-errors` JSON works, since it carries the `trans` block. The
harness needs only base R. For a whole-model comparison, pin the training set
first and check that `sum(trans)` agrees on both sides before looking at
anything else. `--train-samples` pins which files are read, but `nbases` still
caps how much of them is used; leaving R at its 1e8 default once trained R on a
third of the manifest, and that training-data difference looked like a smoother
difference. `dev/concordance/write_reference.R` now stops when `nbases` would
truncate the manifest ([#214](https://github.com/HPCBio/dada2-rs/pull/214)).

## What this dictates

- **LOESS is closed as a source of R divergence.** A future error-model
  difference from R should be looked for in the training set or the
  self-consistency loop first.
- **Do not port the remaining `loess()` surface for completeness.** Nothing
  unported can reach a DADA2 error model. Port a piece only when a decided
  feature needs it, as `span > 1` was ported for [#98](https://github.com/HPCBio/dada2-rs/issues/98).
- **Do not chase R in the degenerate corners.** Few-column data belongs on
  `--errfun binned-qual`, and the open question there is agreement with
  `makeBinnedQualErrfun` ([#120](https://github.com/HPCBio/dada2-rs/issues/120)),
  not with LOESS.
- **Any deliberate divergence starts from an exact baseline.** Robustness,
  alternative weights or a different span are now measurable against a port
  that is exact, so an A/B isolates the change being tested.
