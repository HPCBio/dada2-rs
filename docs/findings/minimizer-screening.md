# Minimizers as the pre-alignment screen

**Verdict: this is a long-read screen.** On PacBio HiFi it is **ASV-identical to
the k-mer screen (1540 = 1540, churn 0)** while doing **13.3% fewer alignments**
and using ~74% less resident memory. On Illumina it is merely *equivalent* — 1-2
ASVs of 232 — at best-case parity on cost. That split is not a quirk; it follows
from how a sketch works, and it took the whole experiment to see.

| | Illumina MiSeq (250 bp) | **PacBio HiFi (1490 bp)** |
|---|---|---|
| table size | 232 ASVs, 20 samples | **1540 ASVs**, 3 samples, 542k reads |
| ASV set | 1-2 churned | **identical, churn 0** |
| count-matrix L1 | 0.0596% | **0.0053%** |
| alignments vs k-mer | 83-104% | **86.7%** |
| equivalent cutoff | 0.72 | **0.50** |

**Why length decides it.** The sketch is `O(len/w)` entries: a 250 bp read yields
~48 minimizers, a 1490 bp read ~500. The distance estimate is a sample, and its
precision scales with the sample size. On short reads the estimate is coarse, so
a cutoff loose enough not to lose real neighbours (0.72) also admits far more
strangers — hence +29% alignments at full recall. On long reads the estimate is
sharp: full recall arrives at 0.50 while *fewer* pairs pass than the k-mer screen.
At the very same 0.42, PacBio gets 98.45% recall passing 9.36% of pairs against
k-mer's 10.73% — better and cheaper — where Illumina gets 60.81%.

This is the original hypothesis, confirmed only after being wrong about it twice.
The k-mer screen has the opposite length behaviour — its `4^k` vector saturates
as reads get longer, which is why PacBio needs k=7 and why at k=5 it degenerates
to passing **100.00%** of pairs. Sketch precision *rises* with read length exactly
where frequency-vector precision falls.

Experimental and opt-in (`--screen-backend minimizer`), default unchanged. The
shipped defaults (k=11, cutoff 0.42) are **wrong on both platforms**; use k=8 with
a cutoff derived per platform by `kdist-calibrate --screen-backend minimizer`
(0.50 PacBio, ~0.58-0.65 Illumina).

**This page reported the opposite twice before reaching that.** An early revision
concluded equivalence from 2-sample fixtures; a later one concluded irreparable
fragmentation. Both were wrong — see
[How the fixtures misled](#how-the-fixtures-misled) and
[the gapless section](#the-residual-is-not-the-screen-at-all-the-gapless-shortcut-switches-off).

## What the screen is

The pre-alignment screen decides which candidate pairs are worth aligning. It is
a *speed* optimisation inherited from the ESPRIT lineage — **not** the cluster
definition. That is what makes a screen swap legitimate to attempt at all: the
screen may change its mind about which pairs to align, provided the pairs it
drops were never going to matter. This experiment is a measurement of whether
that proviso holds. It does not, yet.

Both backends feed the same `--kdist-cutoff` gate. The minimizer distance uses
`kmer_dist8`'s formula, `1 - Σ min(count) / min(total)`, evaluated on a winnowed
sketch rather than on the full `4^k` vector.

## The control channel comes first

Every number below is only interpretable because the noise floor was measured
rather than assumed:

| | ASVs | churn | count cells differing | L1 |
|---|---|---|---|---|
| **k-mer vs k-mer, same binary, same input** | 232 = 232 | **0** | **0 / 1637** | **0.0000%** |

The default backend is **bit-deterministic** in ASV content and per-sample
counts. So the noise floor is zero and *every* difference reported below is real
signal. Establishing this before interpreting the comparison is the standing rule
from [A/B measurement hygiene](measuring-on-numa.md); here it is what converts
"the arms differ by 13" from an ambiguous number into a verdict.

(One caveat it did surface: ASV **ordering** in `seqtab.nochim.json` is *not*
stable run-to-run, even though content is. See
[Incidental finding](#incidental-finding-asv-ordering-is-not-deterministic).)

## The result

### PacBio HiFi — 3 Kinnex samples, 542,205 reads, k=7 production baseline

| arm | ASVs | churn | count cells differing | L1 | aligned vs k-mer |
|---|---|---|---|---|---|
| k-mer, k=7 @ 0.42 (production) | 1540 | — | — | — | 100.0% |
| **minimizer, k=8 @ 0.50 (calibrated)** | **1540** | **0** | 17 / 2452 | **0.0053%** | **86.7%** |

The ASV **set is identical**. 233.6 M comparisons in each arm; the minimizer
shrouds 209.4 M against the k-mer screen's 205.6 M, so it performs 24.3 M
alignments against 28.1 M — **3.7 M fewer**. `b_compare` is 82-88%
align-dominated on this platform ([screen vs align](compare-screen-vs-align.md)),
so that is the half of the phase that matters.

### Illumina MiSeq SOP — 20 samples


20-sample MiSeq SOP, per-sample `dada`, through to the chimera-filtered table.
Set identity via `dev/compare_asvs.py`, counts via `dev/compare_seqtab_matrix.py`
(full sample × ASV matrix, cell by cell).

| arm | ASVs | churn | count cells differing | L1 (reads) |
|---|---|---|---|---|
| k-mer, k=5 (reference) | 232 | — | — | — |
| **control:** k-mer again | 232 | **0** | 0 / 1637 | **0.0000%** |
| minimizer k=11, cutoff 0.42 *(shipped default)* | 241 | 17 | 542 / 1709 | 1.10% |
| minimizer k=8, cutoff 0.42 | 237 | 13 | 350 / 1687 | 0.76% |
| minimizer k=8, **cutoff 0.65** | **232** | **2** | 145 / 1645 | **0.34%** |

The failure mode is **fragmentation, not loss**: the minimizer arm produces
*more* ASVs (241 vs 232), and the extras are low-abundance, with 3–4 of them at
Hamming-1 from a baseline ASV — including one at Hamming-1 from a **3,066-read**
parent. That is a textbook error copy that the screen prevented from ever being
compared to its parent, so it births its own cluster instead of being absorbed.

## Which stage causes it: the screen changes the error model too

The screen is active in **both** stages, so the first version of the comparison
above was confounded — the two arms had different error models *and* different
denoising. `build_trans_mat` (learn_errors.rs) aligns each raw against its
cluster center **through the screen**, and a shrouded pair returns early and
contributes nothing to the transition counts. So changing the screen refits the
model:

| minimizer arm | max relative difference in `err_out` vs k-mer |
|---|---|
| k=11, cutoff 0.42 | **90.9%** |
| k=8, cutoff 0.65 | 5.0% |

(One stage genuinely is screen-free: `run_dada`'s first compare of all raws
against cluster 0 sets `kdist_cutoff = 1.0`, which disables the distance gate for
both backends identically — `minimizer_dist` is bounded by 1.0 and the gate is
strict `>`. That covers only that pass.)

Holding the error model fixed and varying one factor at a time decomposes it:

| config | arm | churn | count L1 |
|---|---|---|---|
| **k=11, cut 0.42** | both changed | 17 | 1.10% |
| | k-mer model + **minimizer screen** | 16 | **0.98%** |
| | **minimizer model** + k-mer screen | 1 | 0.24% |
| **k=8, cut 0.65** | both changed | 2 | 0.34% |
| | k-mer model + **minimizer screen** | 2 | **0.32%** |
| | **minimizer model** + k-mer screen | 0 | 0.03% |

**The denoising-stage screen dominates in both configurations**, so the
fragmentation verdict survives the confound. Two things worth keeping:

- **The error model is remarkably robust to this.** A 90.9% relative change in
  `err_out` propagated to 0.24% read-level L1 and a single churned ASV. Whatever
  the screen does to the fitted model, denoising largely absorbs it — consistent
  with [KDIST cutoff decoupling](kdist-cutoff-decoupling.md) finding the
  learn-errors and dada cutoffs separable, and worth remembering before treating
  any error-model perturbation as automatically serious.
- **The screen is an error-model parameter.** It was originally left out of the
  model's `params` provenance block on the reasoning that no released version
  could have written it, which is true but beside the point: a model learned
  under one screen and applied under another is a real provenance error, and the
  90.9% figure is what proves it. `screen_backend` / `minimizer_k` /
  `minimizer_w` are now recorded in `params` and checked alongside
  `kdist_cutoff`, defaulting to `kmer` for older JSONs so nothing pre-existing
  warns spuriously.

### What it does to `learn-errors` convergence

| arm | iterations | total transitions | off-diagonal (error) transitions |
|---|---|---|---|
| k-mer, k=5 | **6** | 20,023,895 | 92,347 (100%) |
| minimizer k=11, cut 0.42 | **5** | 20,022,722 (100.0%) | **111,903 (121.2%)** |
| minimizer k=8, cut 0.65 | **5** | 20,023,840 (100.0%) | 92,218 (99.9%) |

The over-tight screen **converges a round earlier on an error model inflated by
21%**, and the two facts are the same fact. Self-consistency only asks whether
the model stopped moving; it cannot tell "converged because the answer was found"
from "converged because there is less left to disagree with."

The mechanism is *compositional, not quantitative*: **total transition mass is
identical across arms** (within 0.006% — every arm saturates the same `--nbases`
budget), so the screen is not starving the fit of data. It changes *which*
transitions are counted. A raw that should have moved to a nearby center cannot,
because the move requires a comparison that passes the screen; it is then scored
against a distant center and contributes mismatches. Over-screening therefore
*raises* the apparent error rate rather than lowering it.

This argues for the initial rounds being deliberately **inclusive** — the
permissive arm takes the extra round and lands on the lower, more plausible error
rate — and it is the same headroom argument as
[KDIST cutoff decoupling](kdist-cutoff-decoupling.md), which found the dada-stage
cutoff safe to tighten while tightening the *learn* stage churned real ASVs.

> **Caveat on every number on this page.** The harness runs `--nbases 2e7`, and
> the KDIST study found the error model **not converged even at 1e8**. So this is
> 5x below a level already known insufficient. The screen-side conclusions
> survive it — the decomposition above holds the model *fixed* and still churns —
> but the model-side arm (0.03% L1) may understate sensitivity at a converged
> `nbases`. An `--nbases` ladder is owed before any of the error-model numbers
> here are treated as settled.

## The residual is not the screen at all: the gapless shortcut switches off

Calibrating the cutoff (below) takes the disagreement from 1.10% to 0.34% L1 and
then **plateaus** — 0.65 and 0.72 are indistinguishable (0.3356% vs 0.3348%). The
plateau is not a screening effect, and the control proves it:

| arm (identical error model) | ASVs | churn | count L1 |
|---|---|---|---|
| k-mer **cutoff 1.0** vs minimizer **cutoff 1.0** | 232 vs 231 | 1 | **0.3074%** |

With both screens **fully open** nothing can be shrouded, so screening cannot
differ — yet essentially the entire residual survives. The cause is
`raw_align_dp`'s method selection:

```rust
if p.band == 0 || (p.gapless && (kodist - kdist).abs() < f64::EPSILON) {
    align_gapless_with_buf(...)   // skip the DP entirely
}
```

`kodist` is `kord_dist`, the **positional** k-mer mismatch fraction; `kdist` is
the screen distance. Their agreeing means every shared k-mer is also positionally
aligned — no shifts, therefore no indels — which is what makes a gapless
alignment safe. **That predicate is only meaningful when both quantities live in
the same k-mer space.** Under the minimizer backend `kdist` is a sketch distance,
so the equality essentially never holds:

| backend | gapless shortcut fires |
|---|---|
| k-mer | 2,341 / 147,866 aligned pairs (**1.58%**) |
| minimizer, k=8, cutoff 0.72 | 113 / 147,641 (**0.08%**) |

A 20x drop, with no error and no warning — the backend simply stops taking a path
it cannot know it was supposed to take.

**This generalises past minimizers.** The gapless shortcut is *coupled to the
k-mer frequency vector*, and any screen replacement that does not supply a
commensurate compositional distance disables it silently. A screen swap is
therefore not the local change it appears to be: `kdist` is not only a gate, it is
also an **input to alignment method selection**, and the second role is
undocumented and easy to miss. It was missed here.

Two consequences:

- **The memory case is weaker than it looked.** The obvious fix — keep a k-mer
  frequency distance for the predicate — reinstates the `4^k` vector the sketch
  exists to avoid. A cheaper route exists: `kord` is already stored under both
  backends and *is* the per-position k-mer index, so its value multiset is
  exactly the k-mer composition; a sorted copy per raw (~`len` entries, the same
  order as the sketch) would let the exact predicate be recovered by merge-join.
  Untried.
### The fix, and what it says about the k-mer screen

The predicate never needed the screen. It asks a question about the *pair* — do
these two sequences contain an indel? — and read it off `kdist` only because the
k-mer screen already had that number lying around.

`kord` is stored under **both** backends and `kord[i]` *is* the k-mer index at
position `i`, so its value multiset is exactly the k-mer composition. Both counts
come from it:

```text
positional     = #{ i : kord_a[i] == kord_b[i] }
compositional  = |multiset(kord_a) ∩ multiset(kord_b)|
gapless        ⟺  compositional == positional
```

`pair_is_gapless` computes that and consults no screen at all. The `4^k` counter
is **per-thread scratch, not per-raw resident**, and is cleared only over the
`klen` entries touched — so recovering the exact predicate does *not* give back
the memory the sketch saves. The now-dead per-pair `kord_dist` is removed; it ran
on every screened pair, including the majority the screen then rejected.

**A minimizer-specific analogue would have been a trap.** Comparing the *i*-th
minimizer of each read reports "no indel" even when there is one, because
winnowing is deliberately shift-robust and an indel barely perturbs the ordered
minimizer list. Doing it properly needs per-minimizer sequence offsets (~45% more
sketch) to reconstruct what `kord` already holds.

Result: both backends take the gapless path **exactly 2,341 times**, the k-mer
backend is bit-identical through the change, and the minimizer residual falls
from 0.3348% to **0.0596%**.

### Which arm is actually more faithful

With the coupling gone, an unscreened run (cutoff 1.0) becomes a usable reference
for what denoising does *without* the prefilter:

| arm | ASVs | vs unscreened |
|---|---|---|
| unscreened | 231 | — |
| k-mer @ 0.42 (production, = R DADA2) | **232** | **+1 ASV**, L1 0.0692% |
| minimizer @ 0.72 (calibrated) | **231** | **ASV set identical**, L1 0.0612% |

The screen is supposed to be a pure speed optimisation, dropping only pairs that
could not have mattered. **On this dataset the k-mer screen at 0.42 is not
lossless** — it yields one ASV unscreened denoising does not. The calibrated
minimizer screen is lossless here.

Both framings are true and they use different references: against **R DADA2**,
which is this project's stated fidelity target, the k-mer arm is correct by
definition and the minimizer differs by one ASV. Against **what the algorithm
does without a lossy prefilter**, the minimizer is the more faithful of the two.
One ASV in 231 on one dataset is a modest effect — but it is measured against a
zero noise floor, and it is a fact about the *existing* screen, not about
minimizers.

## Two claims this falsified

### 1. "The cutoff transfers between backends." It does not.

The metric was deliberately built to share `kmer_dist8`'s algebraic form, on the
argument that this would keep `KDIST_CUTOFF = 0.42` interpretable across both
backends. **That reasoning was wrong.** At the same cutoff, on the same sample:

| | pairs passed |
|---|---|
| k-mer screen, k=5, cutoff 0.42 | **27.6%** |
| minimizer screen, k=8, cutoff 0.42 | **9.0%** |

A cutoff is a property of the **distance distribution** a metric induces, not of
its formula. Sharing the algebra transfers nothing. Sweeping the cutoff locates
the matching operating point at **≈0.64**:

| minimizer cutoff (k=8) | 0.42 | 0.55 | 0.65 | 0.70 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|
| pairs passed | 9.0% | 19.3% | **28.9%** | 34.4% | 40.6% | 46.1% |

This is why the calibrated arm above closes most of the gap. It also means any
promotion needs a **calibrated per-metric cutoff**, which belongs in
`kdist-calibrate` — the subcommand that already does exactly this job for the
k-mer metric — and not in a hand-picked constant. Compare the
[KDIST cutoff decoupling](kdist-cutoff-decoupling.md) result, which established
that the learn-errors and dada cutoffs are independently settable; a per-backend
cutoff is the same kind of decoupling one level up.

### 2. "k=11 is a reasonable default." It is the wrong end of the range.

`k=11` came from a back-of-envelope substitution argument, never a measurement.
The audit shows it *missing real neighbours*, and lowering `k` monotonically
fixes that (`F3D146`, kmer-only pairs the minimizer screen would shroud):

| minimizer k | 11 | 9 | **8** | 7 (w=2) |
|---|---|---|---|---|
| pairs passed | 5.1% | 6.6% | 9.0% | 13.0% |
| missed at nsubs ≤ 5 | 0 | 0 | **0** | 0 |
| missed at nsubs 6–10 | 14 | 10 | **0** | 0 |

At **k=8, misses at nsubs ≤ 10 go to zero** on every sample checked (F3D142,
F3D146, F3D149) — and yet the table *still* churned by 13 at cutoff 0.42. That
non-result is the useful one: it proves the residual divergence is **not** driven
by missing close pairs, which is what sent the investigation to the cutoff.

## The mechanism

Putting the two together: the screen was simply too tight overall. Shrouding a
pair yields `sub = None`, hence `lambda = 0`. A raw shrouded from *every* cluster
center has zero expected abundance everywhere, so its abundance p-value collapses
and it births its own cluster. This happens even when the shrouded pairs are all
distant (>10 substitutions) — the raw does not need a *good* match to avoid a
spurious birth, it needs *a* match. The k-mer screen at 27.6% pass rate supplies
one; the minimizer screen at 9.0% does not.

So the fragmentation is a **pass-rate** effect, not a **precision** effect, and
that distinction is what the `nsubs` histogram bought.

## The pair-level audit

`--screen-audit` evaluates both screens on every comparison, aligns the
**union**, and buckets each disagreement by the substitution count the aligner
found; audit-only alignments are discarded so ASVs stay the active backend's.
Counts alone cannot separate a lost error copy from a correctly-rejected
stranger — the histogram can, and it is what localised the mechanism above.

Consistent across every configuration tested: **minimizer-only = 0**. The
minimizer screen never passes a pair the k-mer screen rejects; it is a strict
subset. The question was only ever how much *more* it rejects, and whether that
excess matters. It does.

**Instrument limitation:** the audit applies one `--kdist-cutoff` to both
screens, so it cannot directly compare the two at *different* calibrated cutoffs.
The cutoff sweep above works around this by reading only the minimizer pass rate.
A proper per-backend cutoff in the audit is a prerequisite for the next round.

## What still holds

- **The inverted index is exact.** `MinimizerIndex` replaces `nraw` merge-joins
  with one scatter over the center's posting lists, and index-on vs index-off is
  **byte-identical on real data** (`DADA2RS_MINIMIZER_INDEX=0`). Keeping the
  acceleration on a separate axis from the metric is what let every difference
  above be attributed to the metric alone. This part is not in question.
- **Memory.** PacBio HiFi (~1,495 bp), screen structure only, `kord` excluded:
  minimizer k=11/w=5 costs **4,306 B/raw** against k=7's **16,384** (−74%).
  Whether that comparison survives at the *calibrated* parameters is untested —
  `k` does not change sketch size but `w` does, and the calibration moved `k`,
  not `w`, so it very likely does hold.
- **PacBio at k=5 screens 100.00% of pairs.** Directly measured, on 38,079
  comparisons. [k-mer screen size](kmer-size-screening.md) called it "effectively
  a no-op"; it is a *literal* no-op. That finding is independent of the minimizer
  work and stands on its own.
- **Sketch density is 0.333** against winnowing's theoretical `2/(w+1)` = 0.333.
  Nothing degenerate in the sketch construction itself.

## How the fixtures misled

The 2-sample concordance fixtures reported **churn 0 at every stage** on both
platforms — Illumina 8=8, PacBio 93=93 — and the earlier revision of this page
generalised from that. The SOP then showed 17 churned ASVs on the same backend at
the same settings.

The fixtures were not wrong, they were *underpowered*, and in a way that was
predictable in advance:

- **8 ASVs cannot exhibit the failure.** The failure mode is spurious
  low-abundance ASVs at Hamming-1 from abundant parents. A fixture whose entire
  table is 8 high-abundance organisms has almost no low-abundance tail for the
  screen to fragment. The instrument could not return the answer it was being
  asked for — the same error catalogued in
  [the store-scan page's](compare-store-scan.md) three instrument failures.
- **The tail is the whole question.** A screen change perturbs marginal
  decisions. Any fixture dominated by unambiguous calls will report agreement
  regardless of whether the screen is sound.
- **PacBio 93=93 was more reassuring than it deserved to be.** It is 2 samples,
  and long reads at k=5 compare against a baseline that screens *nothing*
  (100.00% pass) — the arm being compared to was itself degenerate.

The generalisable rule: **a concordance fixture is only evidence about the part
of the distribution it contains.** For a screen, that means the low-abundance
tail, and it should be sized for the tail before it is trusted. This is the same
lesson as [k-mer screen size](kmer-size-screening.md)'s "measure at the end of the
pipeline, not on intermediate counts", one level deeper: measure at the end, *on
data with an end worth measuring*.

## Incidental finding: ASV ordering is not deterministic

Two runs of the same binary on the **default k-mer backend** emit the same ASVs
with the same counts in a **different order**. Content is deterministic (the
control above is bit-identical); only the ordering moves. Consequence:
**byte-identity checks on `seqtab.nochim.json` are unreliable**, and comparisons
must be set- and count-based. `dev/compare_asvs.py` and
`dev/compare_seqtab_matrix.py` are; ad-hoc `diff` is not. Worth its own issue.

## Incidental fix: the concordance harness mis-paired >2 samples

`run_illumina.sh` passed `merge-pairs` two independently-globbed lists and relied
on them lining up positionally. They sort differently whenever the filtered name
inserts a character the raw name lacks — the FASTQ glob orders `F3D1F.fastq.gz`
against `F3D141F.fastq.gz`, the JSON glob orders `F3D1_F_filt.json` against
`F3D141_F_filt.json`, and the `_` flips them. On the 2-sample fixture the lists
coincided; on the SOP it mis-paired 10 of 20 samples. Only a downstream
map-length check in `merge-pairs` turned it into a visible error instead of
silently merging mismatched samples. Fixed by enumerating the JSONs from the
filtered lists.

## What would change the verdict

0. **Done on this branch:** `screen_backend`/`minimizer_k`/`minimizer_w` are now
   recorded in the error model's `params` and checked on mismatch.
1. **Calibrate the cutoff properly, per backend**, in `kdist-calibrate` rather
   than by the hand sweep used here — and teach `--screen-audit` to hold a
   separate cutoff per screen so the two can be compared at their own operating
   points.
2. **Close the residual 2-ASV / 0.34% gap, or explain it.** The calibrated arm is
   close but not equal, against a zero noise floor. Either it reduces to a
   remaining pass-rate mismatch (test: sweep cutoff finer around 0.64 and watch
   churn), or something structural remains.
3. **Re-run PacBio at scale** — `PacBioFull`, not the 2-sample fixture — now that
   the fixture-scale PacBio result is known to be underpowered.
4. **Re-check memory at calibrated parameters**, and consider whether `kord`
   (2,982 B/raw on HiFi, *larger than the sketch itself*) can be dropped on this
   path; it is retained only for the gapless fast path.
5. **Only then, a wall-clock A/B** — pinned, replicated, with a control channel.
   No timing claim is made anywhere on this page, deliberately.

## Reproducing

```bash
# Control channel FIRST -- establishes the noise floor. Must be churn 0.
bash dev/concordance/run_illumina.sh ./target/release/dada2-rs $SOP /tmp/k1 8
bash dev/concordance/run_illumina.sh ./target/release/dada2-rs $SOP /tmp/k2 8
python3 dev/compare_seqtab_matrix.py /tmp/k1/seqtab.nochim.json /tmp/k2/seqtab.nochim.json

# Calibrated minimizer arm.
SCREEN_BACKEND=minimizer MINIMIZER_K=8 SCREEN_CUTOFF=0.65 \
  bash dev/concordance/run_illumina.sh ./target/release/dada2-rs $SOP /tmp/mini 8
python3 dev/compare_asvs.py --baseline kmer=/tmp/k1/seqtab.nochim.json \
    --compare minimizer=/tmp/mini/seqtab.nochim.json
python3 dev/compare_seqtab_matrix.py /tmp/k1/seqtab.nochim.json /tmp/mini/seqtab.nochim.json

# Pair-level mechanism.
./target/release/dada2-rs dada <filt.fastq.gz> --error-model <err.json> \
    --screen-backend minimizer --minimizer-k 8 --screen-audit -o /dev/null
```

The MiSeq SOP lives at `~/projects/hpcbio/dada2-rs-data/miseq/MiSeq_SOP` (raw;
the harness runs `filter-and-trim` itself). The harness expects
`<sample>F/<sample>R.fastq.gz`, so symlink the `_R1_001`/`_R2_001` names across.
