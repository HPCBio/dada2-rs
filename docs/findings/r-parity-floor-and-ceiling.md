# What R parity can and cannot tell you

**Verdict:** agreement with R DADA2 is our sharpest fidelity instrument, and it
has **a floor and a ceiling**. Below the floor, differences are coin flips that
two R releases would also produce — chasing them is chasing noise. Above the
ceiling sit places where R's *actual* behaviour is not its *intended* behaviour,
and we deliberately follow the intent instead. A parity number quoted without
both bounds invites two opposite mistakes: treating an irreducible tie-break as
a bug, and treating a deliberate divergence as a regression.

Concretely, on the pooled 362-sample MiSeq SOP and 30-sample NovaSeq ITS2 runs,
both read directions, we now reproduce R's ASV set **exactly** and differ by
**0 to 4 single reads** across the whole sample × ASV table, and ITS2's merged
table by 5, while the measured floor on the same runs is up to 14 ASVs. Pooled
PacBio HiFi has R's ASV set exactly and differs by 15 single reads, every one
of them from [a gapless test that only x86-64 R
takes](#x86-64-r-and-the-gapless-shortcut); with that test emulated, the
95-sample table and the learned transition counts are identical to R's.
Getting there took removing four *systematic* differences, three in ordering
and one in merging, each of which had been hiding inside what looked like
floor; one of them, [the pool's tie-break](#a-third-instance-the-pools-tie-break),
was attributed to the floor on this very page until it was measured.

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

With the systematic differences removed (see below), the residual against R
`dada(pool=TRUE)`, with pinned error models and inputs in R's order, is:

| pooled table | ASVs only on one side | sample × ASV cells that differ | reads moved |
|---|---|---|---|
| MiSeq SOP 362, forward | 0 | 2 of 205,782 | 2 |
| MiSeq SOP 362, reverse | 0 | 0 of 179,919 | 0 |
| NovaSeq ITS2 30, forward | 0 | 2 of 15,852 | 2 |
| NovaSeq ITS2 30, reverse | 0 | 4 of 16,212 | 4 |
| PacBio HiFi 95 | 0 | 30 of 62,585 (0 emulating x86 R) | 30 (0) |

On 16S and ITS2 every remaining difference is a pair of cells in one sample: a
single read that one implementation assigns to one ASV and the other to a
neighbour. Per-sample totals and the ASV set agree. That is a read-assignment
tie, a level below which ASV is born.

The same holds after merging. NovaSeq ITS2's merged, pre-chimera table has R's
ASV set exactly, and its 5 differing cells, one read each, are the
per-direction ties above carried through `merge-pairs`. That needed a fourth
systematic difference fixed, this one not about order: `merge-pairs` had aligned
with the denoising scores instead of R's `mergePairs` scores, and counted the
overlap differently, so on long ITS2 variants with a 13–15-base overlap it
rejected 93 reads' worth of pairs that R merged
([issue 272](https://github.com/HPCBio/dada2-rs/issues/272)).

PacBio's 30 cells have the same shape, 15 single reads, but they are not
ties. Each is a read 1–2 bases shorter or longer than a large centre, which x86-64
R compares without gaps and we align; see [the ceiling](#x86-64-r-and-the-gapless-shortcut).
The R reference was built on x86-64. An earlier residual of 10 ASVs and 558
reads came from a stale error model, not from denoising
([issue 269](https://github.com/HPCBio/dada2-rs/issues/269)).

### How large the floor is, measured

Running the same input under different, deliberate member orders measures the
floor directly ([issue 157](https://github.com/HPCBio/dada2-rs/issues/157)).
`DADA2RS_MEMBER_ORDER` selects `insertion` (the default), `sorted` (derep order)
or `shuffle:<seed>`, with the error model pinned so only order varies; five
seeds form the null. ASVs that change against `insertion`, per table, with
pooled runs in R's pool order (`--pool-tiebreak first-seen`):

| dataset | mode | ASVs that change | largest | distinct sequences (edit > 12) |
|---|---|---|---|---|
| MiSeq SOP, 362 samples | pooled, forward / reverse | 0–2 / 0 | 16 reads | none; the 2 are one 11-mismatch pair |
| | pseudo / per-sample | ≤ 3 | 16 reads | none |
| PacBio HiFi, 95 samples | pooled | 19–24 | 76 reads | 11% of changes, ≤ 21 reads |
| | pseudo / per-sample | 29–35 / 28–36 | 32 reads | none |
| NovaSeq ITS2, 30 samples | pooled, forward / reverse | 10–14 / 10–14 | 117 reads | none / 18% of changes |
| | pseudo, forward / reverse | 24–34 / 28–39 | 78 reads | none |
| | per-sample, forward / reverse | 25–32 / 26–39 | 97 reads | none / 2% |

The ITS2 rows were re-measured with the run's actual binned-quality anchors
(2,11,25,37); the first measurement used anchors that missed them
([issue 264](https://github.com/HPCBio/dada2-rs/issues/264)). The ranges moved
little — pooled was 4–16, pseudo and per-sample 19–38 — and the largest pooled
change is still 117 reads. The PacBio rows were re-measured with a current
error model, replacing the stale model of
[issue 269](https://github.com/HPCBio/dada2-rs/issues/269): pooled 19–24 (was
17–36), largest change still 76 reads, 63% of changes within 3 edits; pseudo
and per-sample 28–36 (was 28–39), 96% within 3 edits. The MiSeq pseudo and
per-sample row is the original measurement: neither mode pools dereps, so the
pool's tie-break cannot reach it.

Read totals barely move: at most 0.12% of reads (PacBio per-sample), 0.06% on
ITS2. The changes are **renames**, not organisms gained or lost: an ASV named
after a 1–3-edit variant of its centre (a substitution or a homopolymer indel),
at unchanged abundance. That is 81–93% of changes per-sample and pseudo in the
original measurements, and 96% on PacBio re-measured; pooled
runs have more swaps further apart (half of the ITS2 reverse changes), still
between equally abundant pairs. Typically two members tie at `pA = 0` with
equal reads, order picks the centre, and the other cannot reach `OMEGA_A`
against it.

The floor grows with read length and diversity, from 0–2 ASVs on MiSeq V4 to
19–36, about 1% of the table, on full-length PacBio. **The residual against R
is now below anything the floor produces:** on 16S and ITS2 R lands exactly
where our unshuffled order does, and on PacBio so does x86-64 R with its gapless
test emulated, while single shuffles move 2 to 14 ASVs on 16S and ITS2 and
20–24 on pooled PacBio. Against R the shuffles differ by the same counts, 20–24
ASVs and 200–351 cells; our default differs by 0 ASVs and 30 cells, and by none
with the emulation.

`sorted` is a fair stand-in for R's ordering: R and dada2-rs both start from
derep order and share the same `swap_remove` member updates. On every pooled
MiSeq and ITS2 table it matches R exactly as `insertion` does, and it changes
0–2 ASVs on ITS2 per-sample and pseudo. That holds for pooled runs only under
R's pool tie-break; with the old `lexical` pool it did not.

Before primer-trimming length variants were removed, pooled ITS2 also flipped
the names of ASVs of up to 692 reads, each between two length variants of one
molecule. That part is a prep artifact, not this floor: see [Primer trimming
makes length variants](primer-trimming-length-variants.md).

### Why this matters for how parity is read

**A parity claim is only as fine-grained as this floor.** On a saturated pooled
run, a difference smaller than the floor and "we agree with R" are the same
statement. Two DADA2 releases, or the same release on reordered input, would
produce differences of the same kind and magnitude, and the table above gives
that magnitude per platform.

The converse matters as much: **a difference the size of the floor is not
thereby floor.** The pool's tie-break moved 2 ASVs on MiSeq and 8–12 on ITS2,
squarely inside the floor's range, and was a systematic difference all the same.
It showed only because the residual was measured directly against R with the
tie-break as the single variable, not compared with the floor's size.

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
**Our pooled merge did not** — it built the pool in first-seen order. On the
PacBio run, cluster 0's centre sat at position 10,020 while position 0 held a
66,937-read organism that was therefore permanently unbuddable, absent from all
2,818 divisions, while R called it with 106,853 reads. Exactly 1 of 2,810
clusters had the invariant broken, and it was cluster 0.

Fixing it moved us from +2,450 reads and 2,810 pre-chimera ASVs to −21 reads and
2,791 ASVs, the same count as R. **The distinction is the whole point of this
page:** a systematic ordering difference is a bug and must be fixed; the
residual tie-break sensitivity underneath it is a floor and must be recognised.

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

### A third instance: the pool's tie-break

The fix for [issue 219](https://github.com/HPCBio/dada2-rs/issues/219) sorted
the pool by abundance and broke ties **by sequence**, copying the per-sample
rule from `derepFastq`. R's pool is built differently: `combineDereps2` takes
`unique()` over the samples' sequences in input order, then a stable
`order(decreasing = TRUE)`, so **tied uniques keep their order of first
appearance across the inputs**
([issue 260](https://github.com/HPCBio/dada2-rs/issues/260)). With most of a
large pool's uniques tied at 1 or 2 reads, the two rules order nearly the whole
tail differently, and the order decides saturated births.

Measured against R with everything else pinned:

| pooled table | sequence order, vs R | first-seen order, vs R |
|---|---|---|
| MiSeq SOP 362, forward | 2 ASVs; 28 cells, 32 reads | 0; 2 cells, 2 reads |
| MiSeq SOP 362, reverse | identical | identical |
| NovaSeq ITS2 30, forward | 12 ASVs; 55 cells, 432 reads | 0; 2 cells, 2 reads |
| NovaSeq ITS2 30, reverse | 8 ASVs; 42 cells, 91 reads | 0; 4 cells, 4 reads |
| PacBio HiFi 95, stale model | 17 ASVs; 474 cells, 762 reads | 10 ASVs; 404 cells, 558 reads |

The two MiSeq forward ASVs are the pair this page used to cite as the textbook
floor case — two 2-read uniques eleven mismatches apart, both saturating
against cluster 0. Under R's pool order we and R make the same call.

`dada-pooled` now uses first-seen order by default; `--pool-tiebreak lexical`
keeps the old, input-order-independent rule. The cost is R's: pooled results
depend on the order the samples are given. Pass an explicit, byte-sorted list,
not a locale-sorted glob, and keep it with the results; each output records its
`pool_input_index` and the rule in `params.pool_tiebreak`.

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

### x86-64 R and the gapless shortcut

R DADA2's results depend on the CPU it runs on. `raw_align` skips the alignment
when the positional k-mer distance equals the composition distance, which is
meant to show the pair has no indel. On x86-64 R computes the positional
distance with `kord_dist_SSEi`, which compares sequences of unequal length over
the shorter one; everywhere else, `kord_dist` returns -1 for unequal lengths and
the shortcut cannot fire. The SSE version's own comment says it returns -1
too; its length check was lost in DADA2 commit `4a89b96` (2018).

When an indel lies within k bases of either end, the two distances still come
out equal, so x86-64 R lines the pair up without gaps. Every base past the indel
then counts as a substitution, λ collapses (4.3e-8 to 1.6e-17 in the traced
case), and the comparison is dropped as below `E_minmax`. The read then
stays with a small cluster when it is closer to a large one.

dada2-rs keeps the documented behaviour, which is also what R does on ARM.
`DADA2RS_GAPLESS_X86=1` reproduces x86-64 R, for parity runs against an x86-64
reference. On the pooled 95-sample PacBio run, learned and denoised, it takes
the residual from 15 reads to **0 of 62,585 cells**, and the learned `trans` from
10 differing cells to 0 of 1,504. A per-read trace with the emulation agreed
with x86-64 R's on every event (1,450 lines, birth p-values aside). The emulation
is covered by a unit test, not by a CI reference; the full-run result is
recorded in [issue 277](https://github.com/HPCBio/dada2-rs/issues/277), along
with a one-line fix proposed upstream.

## What this dictates

- **Quote parity with its floor.** On pooled runs, report the share of divisions
  at `pA = 0.00e0` alongside any ASV difference. Without it, a reader cannot
  tell a defect from a coin flip. The verbose `Division ... pA=` lines are the
  source; a trace's `members[].pval` is the post-hoc `omega_c` value and is
  **not** the loop's `pA`.
- **Measure the residual, do not infer it from the floor.** A difference inside
  the floor's range can still be systematic. Compare against R at the cell
  level (`dev/compare_seqtab_matrix.py`) with one variable changed at a time.
- **Judge residual differences by shape, not count.** Mirrored exclusive sets at
  matched abundances, with read totals agreeing to ~1 in 2.4M, are a tie-break.
  A one-sided difference, or one concentrated at high abundance, is not. Pairs
  within a few edits at equal abundance are renames.
- **Record the reference's CPU architecture.** R DADA2 gives different answers
  on x86-64 and ARM ([above](#x86-64-r-and-the-gapless-shortcut)). Compare
  against an R reference from a known architecture, and emulate x86-64 when the
  reference was built there.
- **Remeasure the floor on new data.** `dev/run_member_order.sh` and
  `dev/compare_member_order.py` run and summarise the arms; the floor depends on
  amplicon, platform and pooling mode, so these numbers do not transfer.
- **Do not tune toward the floor.** Any change justified by moving a handful of
  saturated-regime ASVs closer to R is unfalsifiable at that resolution. Require
  an effect that clears the floor, or a mechanism — the pool's tie-break had
  both: R's own code, and an effect that vanished on four of five tables.
- **Systematic ordering differences are bugs.** The invariant `r=0 is the
  center` is inherited from input ordering and maintained by no code in either
  implementation. Any new path that builds a `raws` vector must sort descending
  by abundance with **R's** tie-break for that path — lexical within one
  sample (`derepFastq`), first-seen across a pool (`combineDereps2`) — and
  should be tested for it. Whatever picks the centre must take the *first*
  maximum, or a tie at the top re-breaks it.
- **A parity result is evidence only about the path that produced it.**
  Bit-identical `trans` said nothing about `dada-pooled`, because the two do not
  share the code in question.
- **Where R's behaviour is unintended, follow the intent and pin both.** Document
  the divergence, ship the emulation as opt-in, and give each an R reference in
  CI.

## Provenance

- [Issue 219](https://github.com/HPCBio/dada2-rs/issues/219) — pooled merge not
  abundance-sorted; the PacBio measurements in that section
- [Issue 204](https://github.com/HPCBio/dada2-rs/issues/204) — `calc_pA`
  returned 1.0 where R's `ppois(reads-1, 0, lower.tail = FALSE)` gives 0.0 at
  zero expected reads, fixed in the same arc
- [Issue 239](https://github.com/HPCBio/dada2-rs/issues/239) — the centre
  tie-break; the per-sample measurements above
- [Issue 260](https://github.com/HPCBio/dada2-rs/issues/260) /
  [PR 262](https://github.com/HPCBio/dada2-rs/pull/262) — the pool's tie-break;
  the residual tables above
- [Issue 157](https://github.com/HPCBio/dada2-rs/issues/157) — the member-order
  experiment; [issue 264](https://github.com/HPCBio/dada2-rs/issues/264) — the
  ITS2 re-measure with corrected bins
- [Issue 269](https://github.com/HPCBio/dada2-rs/issues/269) — the PacBio
  residual from a stale error model
- [Issue 277](https://github.com/HPCBio/dada2-rs/issues/277) — the remaining
  15 PacBio reads, traced to x86-64 R's gapless shortcut
- [Issue 272](https://github.com/HPCBio/dada2-rs/issues/272) — `merge-pairs`
  alignment scores and overlap counting; the merged ITS2 comparison above
- [Issue 100](https://github.com/HPCBio/dada2-rs/issues/100) — the pseudo-pooling
  re-fit
- [LOESS error-model correctness](loess-error-model-correctness.md) — the other
  half of the fidelity story: where parity *is* achievable, it is achievable to
  ~1e-15
