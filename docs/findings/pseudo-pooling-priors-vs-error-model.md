# Pseudo-pooling: priors, not a re-fitted error model

**Verdict:** our `dada-pseudo` — select priors from round 1, re-run each sample
with them flagged, error model unchanged — is the behaviour to keep. R DADA2
appears to additionally **re-fit the error model between its two pseudo rounds**;
emulating that is measurably worse on the 362-sample MiSeq benchmark on every
axis we can measure: **3,118 fewer reads recovered, 709 more ASVs, and 72× more
ASVs that prior flagging cannot account for**. R emulation is available behind
`--reestimate-err-between-rounds` but is **not** the default.

**The re-fit is unintended, confirmed with DADA2's author.** Re-estimating error
rates *within* a `dada()` call is intended `selfConsist` behaviour; carrying that
re-estimated model out of pseudo-pooling's first step and into its second is not —
the error model supplied by the caller is meant to hold across both steps. So our
priors-only default is not merely a defensible reading of the published
definition, it is the intended semantics, and `--reestimate-err-between-rounds`
emulates a behaviour its author does not intend.

**R's behaviour is confirmed, both by direct observation and at the table level.**
Tracing R's own `dada_uniques` shows `pool="pseudo"` denoising its second round
with an error model the caller never supplied, while `pool=FALSE` uses the supplied
model throughout ([details](#observed-directly-r-round-2-uses-a-different-model)).
Independently, enabling our emulation moves us from 31 ASVs apart from R to 3, and
**30 of R's 31 extra ASVs are precisely the ones the re-fit produces**
([details](#r-really-does-re-fit-the-model-between-rounds)). See
[issue 100](https://github.com/HPCBio/dada2-rs/issues/100).

## Background

Two things prompted this. First, a long-standing observation (CJ, in R, some years
ago): `dada(pool="pseudo")` and a hand-rolled two-round run — round 1 per-sample,
round 2 per-sample with round-1 priors filtered by the same rule — gave
*different* results, with the manual form recovering more. They should have been
identical: same inputs, same `PSEUDO_PREVALENCE = 2` / `PSEUDO_ABUNDANCE = Inf`
selection over the pre-chimera table. The motivation for splitting the run was
parallel per-sample denoising on a cluster.

Second, reading R's `dada()` offers a mechanism. Pseudo-pooling there is not two
`dada()` calls; it is two turns of the self-consistency loop:

| `R/dada.R` | Behaviour |
| --- | --- |
| `300` | `erri <- err` — this turn denoises with whatever `err` currently holds |
| `335` | prior flag is `names(drpi$uniques) %in% c(priors, pseudo_priors)` |
| `371-378` | `err <- errorEstimationFunction(cur)` — refitted from this turn's transitions |
| `391-393` | `if(!pseudo \|\| (pseudo && nconsist >= 2))` — with `pseudo`, turn 1 does not break |
| `399-401` | pseudo priors selected from `makeSequenceTable(clustering)` |

The re-fit at 371 has no `selfConsist` guard and `errorEstimationFunction`
defaults to `loessErrfun`. So round 2 would denoise with a model derived from
round 1 rather than the supplied one — and this is invisible in the return value,
which reports only round 1's `err_in` (line `434`).

## What pseudo-pooling is supposed to do

From the DADA2 documentation:

> The purpose of priors is to increase sensitivity to a restricted set of
> sequences, including singleton detection, without increasing false-positives
> from the unrestricted set of all possible amplicon sequences that must be
> considered by the naive algorithm.

This wording is the basis of the sharpest test on this page. If the error model is
held fixed, the **only** input that changes between rounds is the prior flags.
A round-2 ASV that appears in *no* sample's round-1 output and is *not* a prior
therefore has no legitimate source — it is an off-target call from the
unrestricted space, precisely what the restricted-prior design exists to avoid.

## Setup

`dev/run_pseudo_equivalence.sh` (with `REESTIMATE=1`) +
`dev/compare_pseudo_equivalence.py`. Four arms over one shared set of derep JSONs
and one shared error model, so only the plumbing under test can differ:

| Arm | What it runs |
| --- | --- |
| A | `dada-pseudo` (default) |
| B | manual: `dada` → `make-sequence-table` → `seq-table-to-fasta --prevalence 2` → `dada --prior` |
| C | *(only if the prior sets differ)* round 2 with arm A's priors, to separate selection from denoising |
| D | `dada-pseudo --reestimate-err-between-rounds` |

Data: the F1000 'MiSeqSOP' benchmark, 362 samples → 724 outputs (R1 and R2
denoised separately). A 20-sample subset was used for development and is reported
where it differs instructively.

!!! warning "`total_reads` is not a recovery measure"
    A cluster's abundance includes members that failed the `OMEGA_C` attribution
    test, so a dada output's `total_reads` always equals the input total and is
    **invariant** across configurations. An early version of this analysis
    reported "+0 reads recovered" from it — an artifact of a constant, not a
    result. The quantity that varies is reads belonging to uniques that failed
    (`map[i] is None`), whose counts come from the derep; that is what
    `--derep` enables in the comparator.

## Result 1: our two forms are equivalent

Arm A vs arm B, 724 outputs: **0 differ**. Identical prior sets (1,138
sequences), identical per-sample ASVs, abundances and unique→cluster maps —
177,456 ASVs and 6,702,610 reads exact.

This is equivalence *by construction* (`dada-pseudo` round 2 is `mark_priors` +
`dada_uniques`, exactly what `dada --prior` does), and it is now verified at
production scale. Consequence: **the cluster-parallel workflow is sound** — round 1
per sample, collect priors centrally, round 2 per sample, with no penalty relative
to the single built-in command.

!!! note "A green equivalence result can be vacuous"
    On a 2-sample fixture the arms also matched — but round 2 equalled round 1,
    so the priors were inert and the prior path was never exercised. The
    comparator now detects and reports that case. On the benchmark, priors added
    +140 to +189 ASVs per sample, so the path was genuinely tested.

## Result 2: sensitivity gain over per-sample calling

Round 1 *is* the per-sample baseline, which makes the gain directly measurable:

| | reads lost to failed uniques | rescued vs per-sample |
| --- | --- | --- |
| round 1 (per-sample) | 89,702 (51,257 uniques) | — |
| **pseudo, default** | **20,819** | **+68,883 (−76.8%)** |
| pseudo, re-estimated | 23,937 | +65,765 |

Pseudo-pooling rescues ~69k reads, about 1% of the 6.7M total. This is the
documented purpose working as described.

## Result 3: re-estimating the model is worse on every axis

Arm A vs arm D:

| | default | re-estimated |
| --- | --- | --- |
| total ASVs | 177,456 | 178,165 (**+709**) |
| failed uniques | 18,335 | 21,119 (+2,784) |
| reads lost | **20,819 (0.311%)** | 23,937 (0.357%) |
| samples differing | — | 718/724 |
| reads reassigned among shared ASVs | — | 11,864 |

The re-fit **fragments clusters and loses reads**: +709 ASVs while recovering
3,118 fewer reads, handing back 4.5% of the sensitivity gain above.

## Result 4: ASV provenance — the decisive measure

Round-2 ASVs absent from that sample's round 1, classified:

| | rescued via the prior set (intended) | in no round 1, not a prior |
| --- | --- | --- |
| default | 98,712 | **1** |
| re-estimated | 99,289 | **72** |

**72× more unexplained ASVs with re-estimation.** These are off-target calls from
the unrestricted sequence space — the false positives the restricted-prior design
is specifically meant to exclude. Combined with the read counts, the two
behaviours are not a trade-off: priors-only is better on sensitivity *and* on
specificity.

### Why the default is 1 and not 0

On the 20-sample subset this count was **0**, and it was tempting to treat
"unexplained == 0" as an invariant. At scale it is 1, and the mechanism is
understood: flagging a prior changes the partition, so a *non-prior* unique can
newly clear `OMEGA_A` because its cluster or center moved. Abundance p-values are
computed against the partition, so this second-order effect is expected in rare
cases.

**Consequence:** unexplained-ASV count is a quantity to *monitor with a bound*,
never a hard equality assertion — at 1 in ~98,713 it is noise, while 72 is
signal. An `== 0` assertion derived from the small subset would have failed on
real data.

## Observed directly: R round 2 uses a different model

`dev/probe_r_pseudo_err.R` traces `dada2:::dada_uniques` — the Rcpp entry point
every per-sample inference passes through — and records the `err` matrix handed to
it on each call, comparing against the matrix the caller supplied. No dada2 source
changes, just a `trace()`. Two shallow samples, because this is a binary property
of the code path rather than of the data.

| run | `dada_uniques` calls | error model used |
| --- | --- | --- |
| `pool=FALSE` (control) | 2 | both the supplied matrix, Δ = 0 |
| `pool="pseudo"` | 4 | round 1: supplied, Δ = 0 · **round 2: different, Δ = 4.4e-2** |

The control matters: it shows the deviation is specific to pseudo-pooling, not
something `dada()` does generally. `selfConsist` is off in both runs, so a caller
has every reason to expect the supplied matrix to be used throughout.

This is the mechanism, not an inference from it. Reproduce with:

```bash
Rscript dev/probe_r_pseudo_err.R          # uses the data/dada2 fixtures
```

## R really does re-fit the model between rounds

The decisive test does not require patching R. If R re-fits, its ASV table should
sit closer to our `--reestimate-err-between-rounds` arm than to our default. Run
on the same MiSeq data via `dev/benchmark/bench_pooled.py --pool pseudo --run-r`
with `--loess-preset r-dada2` on both arms (so the loess surface is not a second
difference), compared with `dev/compare_pseudo_vs_r.py` on the **pre-chimera**
tables:

| | ASVs | Jaccard vs R | shared | R-only | ours-only |
| --- | --- | --- | --- | --- | --- |
| R `pool="pseudo"` | 781 | — | — | — | — |
| our default | 752 | 0.9579 | 750 | **31** | 2 |
| our `--reestimate…` | 782 | **0.9962** | 780 | 1 | 2 |

Our own two arms differ by 0.0434 Jaccard, and the two distances to R differ by
0.0383 — 88% of that gap, so the test has real power rather than splitting noise.

Both controls check out in the run's own artifacts: the two arms' error models are
**byte-identical** (`errors_fwd/rev.json`), so the flag is the only difference
between them, and both were fitted with `surface: interpolate, cell: 0.2` and
R's `[1e-7, 0.25]` clamps — i.e. R's `loessErrfun` surface, not our `direct`
default. The loess confound is therefore controlled rather than merely bounded by
the residual.

The set-level detail is what makes it conclusive:

| | count |
| --- | --- |
| ASVs R has that our default lacks | 31 |
| ASVs the re-fit produces that our default lacks | 32 |
| **overlap** | **30 (96.8% of R's extras)** |
| R extras the re-fit does not produce | 1 |
| re-fit extras R lacks | 2 |

**The between-round re-fit accounts for ~97% of the pseudo-mode ASV difference
between the two implementations.** Those shared extras are low-abundance — median
6 reads, max 22, 14 of 30 at ≤5 reads — i.e. the borderline class, consistent with
the 72 prior-unexplained ASVs measured above.

The residual after enabling emulation (3 ASVs in ~783, 0.4%) bounds *everything*
else — alignment kernel, tie-breaks, loess surface — and is a useful concordance
result in its own right.

### Chimera removal absorbs most of it

The same comparison on post-chimera tables:

| | ASVs | Jaccard vs R | R-only |
| --- | --- | --- | --- |
| R | 415 | — | — |
| our default | 409 | 0.9855 | 6 |
| our `--reestimate…` | 416 | 0.9928 | 1 |

The 31-ASV gap collapses to 6: roughly **80% of the extra ASVs the re-fit produces
are removed as chimeras**. Two consequences — the end-to-end impact on a final
table is far smaller than the denoising-stage difference, and the fact that these
ASVs are predominantly chimera-like is further evidence they are spurious rather
than rescued signal.

## What this dictates

1. **Keep priors-only as the default.** Every measured axis supports it: more
   reads recovered, fewer ASVs, ~1/72nd the unexplained ASVs. It is also what the
   published definition of pseudo-pooling describes.
2. **Keep R emulation opt-in, and document what it emulates.**
   `--reestimate-err-between-rounds` exists to compare behaviours, and records
   itself in the output's `params` block (`reestimated_err_between_rounds`) because
   `--error-model` alone no longer describes what round 2 used. Flag-off output is
   byte-identical to before the flag existed. Its help text should say it
   reproduces R's *unintended* propagation, not an alternative definition of
   pseudo-pooling — otherwise it reads as a legitimate choice between two designs.
3. **The concordance-vs-design tension is resolved, not balanced.** This was
   previously framed as a judgement call: our default is better on every metric we
   can compute (reads recovered, ASV count, prior-explained provenance) and matches
   the published description, while matching R would close a known ~4% pseudo-mode
   ASV gap (~1.5% post-chimera). That framing assumed R's behaviour was a
   legitimate alternative definition. It is not — it is unintended. So the residual
   ASV gap against R in pseudo mode is **expected and correct**, and should be
   documented as such rather than treated as a defect to close. Concordance testing
   in pseudo mode should compare against R *with the re-fit suppressed*, or accept
   a known offset.
4. **The upstream question is answered; a report is now optional.** The chain was
   already complete — a mechanism in `dada.R`, the code path observed executing,
   the predicted output direction, the specific ASVs, both controls verified. The
   one open question was whether the re-fit is intended, and it is not. What
   remains is a decision about whether to file, not further evidence-gathering.
   Worth noting the behaviour is invisible from R: it is a single `dada()` call,
   and the return value reports only round 1's `err_in`.
5. **Bound, do not assert, the provenance check.** See above.
6. **Cover pseudo mode in the concordance guardrail**, which currently exercises
   per-sample and pooled.

## Why this was hard to see from R

R's entire pseudo run is a single `dada()` call, so round-1 output is never
exposed and the round-2-vs-round-1 comparison cannot be made from outside. Only a
comparison against a separate per-sample run reaches it. Our two-round
decomposition — `--priors-out` plus the manual recipe — makes that comparison
routine, which is why the provenance measure was available here at all.

## Reproducing

```bash
# input must be filtered first: non-ACGT aborts the run (issue 101)
dada2-rs filter-and-trim --fwd sample.fastq.gz --filt filt/sample.fastq.gz \
    --trunc-len 240 --max-n 0 --max-ee 2 --trunc-q 2

REESTIMATE=1 THREADS=24 bash dev/run_pseudo_equivalence.sh \
    "filt/*.fastq.gz" ./pseudo_eq_out loess 16 2
```

Writes `report.txt` (equivalence) and `report.reestimate.txt` (the R-emulation
A/B). Budget roughly 4–5× a single `dada-pseudo` run; the arms are sequential, so
peak memory stays that of one run.

When comparing pre-existing outputs by hand, note that the comparator joins dada
outputs to dereps on the output's recorded `input_file`, not on sample names —
sample names and derep filenames need not resemble each other. `--derep` is
variadic, so keep it last on the command line.
