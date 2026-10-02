# What R parity can and cannot tell you

**Verdict:** agreement with R DADA2 is our sharpest fidelity instrument, and it
has **a floor and a ceiling**. Below the floor, differences are coin flips that
two R releases would also produce — chasing them is chasing noise. Above the
ceiling sit places where R's *actual* behaviour is not its *intended* behaviour,
and we deliberately follow the intent instead. A parity number quoted without
both bounds invites two opposite mistakes: treating an irreducible tie-break as
a bug, and treating a deliberate divergence as a regression.

Concretely, on a 95-sample pooled PacBio HiFi run we match R **exactly** at
pre-chimera (2791 = 2791) and sit **21 reads apart in 2,385,908** — and that
remainder is not reducible, because **1,596 of 2,818 divisions in that run
decide on a tie-break rather than on a statistic**.

## The floor: saturated births make ordering the comparator

DADA2 spawns a new cluster from the member with the smallest abundance p-value.
In `get_pA`, a raw whose expected-read count is small enough returns a p-value
that **underflows to exactly `0.00e0`**. When it does, the statistic that is
supposed to rank candidates carries no information at all, and `b_bud`'s
tie-break takes over: lowest p, then **most reads**, then **lowest position in
the cluster's member list**.

On the pinned 95-sample PacBio run, **1,596 of 2,818 divisions (57%)** were
decided in that regime. At the abundances involved — 7 to 76 reads — the
reads term ties too, so the decision falls through to position, which is an
artifact of input ordering that neither we nor R DADA2 ever chose deliberately.

### What that costs, measured

After [issue 219](https://github.com/HPCBio/dada2-rs/issues/219) fixed the one
*systematic* ordering difference (see below), the residual against R is:

| | dada2-rs | R DADA2 |
|---|---|---|
| pre-chimera ASVs | **2791** | **2791** |
| post-chimera ASVs | 2045 | 2046 |
| reads placed | 2,385,887 | 2,385,908 |

21 reads in 2.39M, 0.0009%. The pre-chimera exclusive sets **mirror each
other** — 10 ASVs on each side, abundances 7–76 on both, with 9/12/13/17/21/24/76
appearing in *each* list, and read totals agreeing to 1 in 2,414,418. They are
the same organisms resolved to a different base: one Cx5 homopolymer indel, one
substitution, one hamming-5 pair.

That is the signature of a coin flip, not of a defect. Chimera removal was
separately confirmed exactly equivalent on these tables by cross-feeding both
directions — 2046 = 2046 on R's input, 2045 = 2045 on ours, zero differences —
so even the post-chimera gap is input-driven rather than an implementation
difference.

### How large the floor is, measured

Running the same input under different, deliberate member orders measures the
floor directly ([issue 157](https://github.com/HPCBio/dada2-rs/issues/157)).
`DADA2RS_MEMBER_ORDER` selects `insertion` (the default), `sorted` (derep order)
or `shuffle:<seed>`, with the error model pinned so only order varies; five
seeds form the null. ASVs that change between arms, per table:

| dataset | mode | ASVs that change | largest | distinct sequences (edit > 12) |
|---|---|---|---|---|
| MiSeq SOP, 362 samples | pooled / pseudo / per-sample | 0–3 | 16 reads | one 11-mismatch pair, pooled |
| PacBio HiFi, 95 samples | pooled | 16–30 | 76 reads | 8% of changes, ≤ 18 reads |
| | pseudo / per-sample | 28–39 | 26 reads | none |
| NovaSeq ITS2, 30 samples | pooled | 4–16 | 117 reads | none |
| | pseudo / per-sample | 19–38 | 97 reads | none |

Read totals barely move: at most 0.11% of reads (PacBio per-sample). The
changes are **renames**, not organisms gained or lost: an ASV named after a
1–3-edit variant of its centre (a substitution or a homopolymer indel), at
unchanged abundance. That is 82–92% of changes per-sample and pseudo; pooled
runs have more swaps 4–12 edits apart (41% on PacBio, 48% on ITS2 reverse
reads), still between equally abundant pairs. Typically two members tie at
`pA = 0` with equal reads, order picks the
centre, and the other cannot reach `OMEGA_A` against it. On ITS2 the traced
cases were doubletons, 28 or more orders of magnitude short of the threshold.

The floor grows with read length and diversity, from a handful of ASVs on
MiSeq V4 to 1–1.6% of the table on full-length PacBio. **The rs-vs-R residual
on pooled PacBio, about 7 ASVs of the same kinds, sits below the floor of
16–30**, so it is not evidence of a defect.

`sorted` is a fair stand-in for R's ordering: R and dada2-rs both start from
derep order and share the same `swap_remove` member updates. It changes 0 ASVs
on ITS2 per-sample and pseudo, 2–4 pooled, and falls inside the shuffle null on
PacBio, where large partitions scramble member order in both
implementations.

Before primer-trimming length variants were removed, pooled ITS2 also flipped
the names of ASVs of up to 692 reads, each between two length variants of one
molecule. That part is a prep artifact, not this floor: see [Primer trimming
makes length variants](primer-trimming-length-variants.md).

### Why this matters for how parity is read

**A parity claim is only as fine-grained as this floor.** On a saturated pooled
run, "we differ from R by 7 ASVs" and "we agree with R" are the same statement.
Two DADA2 releases, or the same release on reordered input, would produce
differences of the same kind and magnitude, and the table above gives that
magnitude per platform.

It also bears on how much weight exact R equivalence can carry as a
*correctness* criterion. Where the tie-break decides, R's answer is not more
correct than ours — it is the answer that R's input ordering happened to
produce. Fidelity here means reproducing R's *procedure*, including its
arbitrary parts, not discovering a truth R has and we lack.

### The part that was a real bug

None of the above excuses a *systematic* ordering difference, and there was one.
`b_bud` scans candidates with `for r in 1..` under the comment `r=0 is the
center` — R has the identical loop and comment (`cluster.cpp:284`). Neither
implementation ever *moves* the centre to position 0; the invariant is inherited
from the input, because R's `combineDereps2` ends with
`order(derepCounts, decreasing = TRUE)`.

Our per-sample derep does the same (with a lexical tie-break, since issue #4).
**Our pooled merge did not** — it built the pool in first-seen order. On this
run, cluster 0's centre sat at position 10,020 while position 0 held a
66,937-read organism that was therefore permanently unbuddable, absent from all
2,818 divisions, while R called it with 106,853 reads. Exactly 1 of 2,810
clusters had the invariant broken, and it was cluster 0.

Fixing it moved us from +2,450 reads and 2,810 pre-chimera ASVs to −21 reads and
an exact 2,791. **The distinction is the whole point of this page:** a
systematic ordering difference is a bug and must be fixed; the residual
tie-break sensitivity underneath it is a floor and must be recognised.

!!! note "Why it hid for so long"
    `learn-errors` never touches the pooled merge — it runs per-sample, as R's
    `learnErrors` does — so every error model was built through the correctly
    sorted per-sample derep. Our `trans` matrices came out bit-identical to R
    (0/656 cells) on the very runs containing this bug. **Our strongest
    fidelity evidence was produced by the one path the bug could not reach.**
    A parity result is evidence about the code path that produced it, and
    nothing else.

!!! note "The same underflow reaches a second decision"
    `OMEGA_C`'s final pass tests `p >= omega_c` on the same `calc_pA` output.
    At `p = 0.00e0` that is **false** for `dada`'s `1e-40` (reads left unplaced)
    and **true** for `learn-errors`' `0` (reads placed), so the two defaults we
    ship diverge precisely on the saturated set. There the underflow is not a
    tie to be broken but a step change, which is a different problem from the
    one above — tracked as
    [issue 228](https://github.com/HPCBio/dada2-rs/issues/228).

### A second instance: the centre tie-break

The same invariant broke a second way. `assign_center` makes the most abundant
raw the centre; when the top two tie, R's strict `>` keeps the **first**
(`cluster.cpp:378`) and ours kept the **last** (`max_by_key`). The centre then
sat off position 0 and the raw there could never bud
([issue 239](https://github.com/HPCBio/dada2-rs/issues/239)).

It needs a count tie at the top of a derep, so pooled runs are immune in practice
(the 362-sample MiSeq pool's top two are 251,498 and 211,477 reads). Per sample
it fired on **4 of 724** MiSeq inputs (three samples), and it **removed real
reads from the final table**:

| sample | tied top | merged reads, before | after | R DADA2 |
|---|---|---|---|---|
| M3D147 | 176 = 176 (reverse) | 1,875 (70 ASVs) | 2,109 (71) | 2,109 (71) |
| M2D19 | 2 = 2 (reverse) | 0 (no pairs merge) | 7 (1) | 7 (1) |
| M3D149 | 1 = 1 (both) | 0 | 1 (1) | 1 (1) |

After the fix our merged table matches R's exactly, **73 of 73 (sample, ASV)
cells**. R ran the same chain (`dada`, `mergePairs`, `makeSequenceTable`) with
the same error models, so only the tie-break differs.

The loss happens at denoising. The unbuddable position-0 raw is its own
organism, so it fails `OMEGA_C` against the other tied centre and its reads are
left unassigned; in M3D147 that dropped a 234-read organism, 11% of the sample's
merged reads. (Before [issue 204](https://github.com/HPCBio/dada2-rs/issues/204)
the same reads were instead counted toward the wrong ASV.)

No fixture or benchmark sample had tied top counts, so no parity comparison
could see this. `dev/top_ties.py` finds the inputs where it applies. None of the
three samples was in that run's `learn-errors` training set, which filled its
`--nbases` budget from earlier files, so its error models were unaffected.

## The ceiling: where we follow the intent, not the behaviour

Parity is a means, not the goal. Where R's implemented behaviour diverges from
what it is meant to do, matching it would propagate a defect.

The worked case is **pseudo-pooling**. R's `dada(pool="pseudo")` re-fits the
error model between its two rounds — confirmed by direct observation of R's own
`dada_uniques`, and at the table level. Per DADA2's author that re-fit is **not
intended**: re-estimating rates *within* a `dada()` call is `selfConsist`
behaviour, but carrying that model out of pseudo-pooling's first step into its
second is not, and the caller's error model is meant to hold across both.

So our `dada-pseudo` keeps the supplied model for both rounds. Emulating R is
available behind `--reestimate-err-between-rounds` and is **not** the default,
because it is measurably worse on the 362-sample MiSeq benchmark on every axis
we can measure (3,118 fewer reads recovered, 709 more ASVs, 72× more ASVs prior
flagging cannot account for). Full account:
[Pseudo-pooling: priors, not a re-fitted error model](pseudo-pooling-priors-vs-error-model.md).

**Both behaviours are pinned in CI**, each against an R reference built the same
way, so neither the default nor the emulation can drift unnoticed. That is the
pattern to follow whenever we depart from R: do not merely document the
divergence, pin both sides of it.

## What this dictates

- **Quote parity with its floor.** On pooled runs, report the share of divisions
  at `pA = 0.00e0` alongside any ASV difference. Without it, a reader cannot
  tell a defect from a coin flip. The verbose `Division ... pA=` lines are the
  source; a trace's `members[].pval` is the post-hoc `omega_c` value and is
  **not** the loop's `pA`.
- **Judge residual differences by shape, not count.** Mirrored exclusive sets at
  matched abundances, with read totals agreeing to ~1 in 2.4M, are a tie-break.
  A one-sided difference, or one concentrated at high abundance, is not. Pairs
  within a few edits at equal abundance are renames; compare a difference's
  size against the measured floor for its platform and mode before treating it
  as a finding.
- **Remeasure the floor on new data.** `dev/run_member_order.sh` and
  `dev/compare_member_order.py` run and summarise the arms; the floor depends on
  amplicon, platform and pooling mode, so these numbers do not transfer.
- **Do not tune toward the floor.** Any change justified by moving a handful of
  saturated-regime ASVs closer to R is unfalsifiable at that resolution. Require
  an effect that clears the floor, or a mechanism.
- **Systematic ordering differences are bugs.** The invariant `r=0 is the
  center` is inherited from input ordering and maintained by no code in either
  implementation. Any new path that builds a `raws` vector — a new pooling mode,
  a new merge — must sort descending by abundance with a deterministic
  tie-break, and should be tested for it. Whatever picks the centre must take
  the *first* maximum, or a tie at the top re-breaks it.
- **A parity result is evidence only about the path that produced it.**
  Bit-identical `trans` said nothing about `dada-pooled`, because the two do not
  share the code in question.
- **Where R's behaviour is unintended, follow the intent and pin both.** Document
  the divergence, ship the emulation as opt-in, and give each an R reference in
  CI.

## Provenance

- [Issue 219](https://github.com/HPCBio/dada2-rs/issues/219) — pooled merge not
  abundance-sorted; the measurements above
- [Issue 204](https://github.com/HPCBio/dada2-rs/issues/204) — `calc_pA`
  returned 1.0 where R's `ppois(reads-1, 0, lower.tail = FALSE)` gives 0.0 at
  zero expected reads, fixed in the same arc
- [Issue 239](https://github.com/HPCBio/dada2-rs/issues/239) — the centre
  tie-break; the per-sample measurements above
- [Issue 157](https://github.com/HPCBio/dada2-rs/issues/157) — the member-order
  experiment that measured the floor above
- [Issue 100](https://github.com/HPCBio/dada2-rs/issues/100) — the pseudo-pooling
  re-fit
- [LOESS error-model correctness](loess-error-model-correctness.md) — the other
  half of the fidelity story: where parity *is* achievable, it is achievable to
  ~1e-15
