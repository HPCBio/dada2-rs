# Minimizers as the pre-alignment screen

**An experimental alternative to k-mer screening that shows promise where a lot
of screening happens.** What to expect depends on the platform, and the gap is
wide enough to state up front:

* **PacBio HiFi — exactly identical, and 18.7% faster.** At a calibrated cutoff
  of 0.45 or above, 95 samples reproduced the k-mer screen's ASV calls, counts
  and reads bit for bit, in every sample individually — while doing **22% fewer
  alignments** and finishing in **59.29s against 72.94s** (3 reps, control 0.2%,
  disjoint ranges). The k-mer screen is simply over-provisioned on this
  platform. This is the one workload where the backend costs nothing at all.
* **Illumina — mostly concordant, not identical.** ASV sets and counts agree
  closely on every dataset tested, with real differences of typically a fraction
  of a percent of reads and a small number of low-abundance ASVs. No cutoff was
  ever exactly concordant.

So it is a drop-in replacement on HiFi and not on Illumina. Why the platforms
differ is not established; the working hypothesis is that HiFi's error profile
leaves near neighbours far enough apart that both screens resolve them
identically, while Illumina's dense 1-2 mismatch clouds sit exactly where the two
screens' distance distributions diverge. That is untested.

Where it pays off is where the screen dominates — diverse pools, in which most
pairs are dissimilar, so the screen runs on every comparison while the aligner
runs on almost none. And the k-mer screen's share is not a constant; it grows with
the working set, because the frequency vector is `4^k` bytes per raw and its cost
is bound by whether those vectors fit cache:

| workload | build | k-mer pass rate | **k-mer screen ns/comp** | **k-mer screen share** | minimizer ns/comp |
|---|---|---|---|---|---|
| MiSeq SOP 16S, one sample | dev | 26.8% | 44 | **0.9%** | ~20 |
| PacBio HiFi, per-sample | dev | 9.9% | 660 | 4.0% | 20 |
| **PacBio HiFi, per-sample, 95 samples** | native | **14.07%** | **2460** | **12.2%** | **32** |
| MiSeq 16S, pooled (#127) | native | 25.0% | 841 | 16.7% | — |
| NovaSeq soil 16S, per-sample | native | 5.52% | 454 | **29.8%** | 30 |
| NovaSeq ITS2, per-sample | native | 1.74% | 308 | **43.9%** | 30 |
| **NovaSeq ITS2, pooled** | native | **0.70%** | **1305** | **76.5%** | **33** |

> **Two build configurations appear here.** Rows marked *native* are
> `RUSTFLAGS="-C target-cpu=native" --profile release-native` on the cluster;
> *dev* rows are a plain `--release` developer build on Apple Silicon. That
> matters for the k-mer column specifically: `kmers.rs`'s `sum-of-min` sweep is
> written to auto-vectorise under `target-cpu=native`, so a dev build understates
> what the k-mer screen can do, while the minimizer's branchy, scatter-write path
> gains little from vectorisation. The *native* rows are therefore the fair
> comparison, and they are the ones every screen-share and speed claim on this
> page rests on. The dev rows are included for the working-set trend only, where
> the effect is an order of magnitude larger than any codegen difference — and if
> anything a native dev build would push the 44 ns figure *lower*, widening the
> range rather than narrowing it.

**The k-mer screen goes 44 -> 1305 ns/comp (30x) across these; the minimizer goes
20 -> 33 (1.7x).** The frequency vector is `4^k` bytes per raw and its cost is
bound by whether those vectors fit cache; the sketch is ~10x smaller and the
inverted index localises access, so pool size barely registers. On pooled ITS2 the
k-mer screen is **76.5% of `b_compare`** — the single largest cost in the pipeline
— and the minimizer replaces it with 9.3%.

What that is worth, at *matched alignment work* so the screen is the only variable:

| dataset | mode | cutoff | speedup at matched alignments | accuracy cost |
|---|---|---|---|---|
| MiSeq SOP 16S | per-sample | 0.70 | ~0% (screen is 0.9%) | churn 0 |
| **PacBio HiFi, 95 samples** | per-sample | **0.45** | **18.7%** (control 0.2%) | **none — exactly identical** |
| PacBio HiFi, 95 samples | per-sample | 0.42 | **34.0%** | 3 of 2409 pooled ASVs, −16 reads |
| PacBio HiFi, 95 samples | per-sample | 0.48 | 3.8% | none — exactly identical |
| NovaSeq soil 16S | per-sample | 0.65 | **16.6%** | churn 277/17398, L1 0.842% |
| NovaSeq ITS2 | per-sample | 0.62 | **15%** | churn 7/3808, L1 0.365% |
| **NovaSeq ITS2** | **pooled** | 0.62-0.63 | **28.3%** (3 reps, control 0.5%; 20% on an exclusive node) | churn 18/3028, L1 0.335-0.769% |

**The right cutoff is the one that matches the k-mer screen's PASS RATE**, and it
must be derived per dataset because that rate spans 0.70%-26.8% across the five
workloads above. Matched-pass predicted the sweep optimum in every case:

| dataset | k-mer pass | minimizer cutoff matching it | sweep optimum |
|---|---|---|---|
| ITS2 pooled | 0.70% | 0.62 | 0.62 |
| ITS2 per-sample | 1.93% | 0.62 | 0.62 |
| soil 16S per-sample | 5.52% | 0.65 | 0.65 |
| **pooled 16S** | **4.54%** | **0.65** | **0.65** |
| PacBio HiFi, 95 samples | 14.07% | 0.50 | **0.45** (one step tighter) |
| MiSeq SOP | 27.5% | **0.64** | 0.70 — but see below |

**Four of these are real tests; one is not, and the MiSeq SOP row was wrong.**
An earlier revision gave MiSeq SOP a matched cutoff of ~0.80, which contradicts
this page's own pass-rate table further down — 0.65 already passes 28.9% there,
above the 27.5% target — so the matching value is ~0.64, and
`kdist-calibrate --derive-cutoff` computes exactly 0.64 on three SOP samples.
Against a swept optimum of 0.70 that is a 6-point miss.

It is also a workload where **nothing is at stake**: the screen is 0.9% of
`b_compare` on the SOP and churn is 0 across a wide range, so "the optimum" there
means "a value that churned nothing", not a performance optimum. The SOP cannot
adjudicate this rule in either direction.

**PacBio is the exception that does matter**: matched-pass
picked 0.50, and 0.45 is equally exact for 22% fewer alignments — so the rule
landed inside the right region but at its expensive end. On a platform with a wide
exactly-concordant plateau, matching the pass rate is a safe starting point and
not the answer; only the ASV table finds the cheap edge of the plateau. See also
[calibrating on read retention](#calibrating-on-read-retention-the-cutoff-is-064),
which converges on ~0.64 across every screen-dominated dataset and is the smoother
signal of the two.

That is what `analyze_kdist_curves.py` reports, and it is a better target than
100% recall — which picked 0.72 on Illumina and 0.50 on PacBio, one of them
20 points too loose.

Experimental and opt-in (`--screen-backend minimizer`), default unchanged. The
shipped defaults (k=11, cutoff 0.42) are **wrong everywhere**; use **k=8** with a
cutoff swept per dataset — **0.62** on both ITS2 modes, **0.45** on PacBio HiFi,
0.62-0.70 on 16S. The per-backend default this branch ships (0.63) is calibrated
for Illumina and is 0.18 too loose for HiFi, where it costs ~60% more alignments
than 0.45 for identical output.

**This page reached that verdict after being wrong three times**: it first
concluded equivalence from 2-sample fixtures, then irreparable fragmentation, then
"a long-read screen" from PacBio alone. The last was closest but still wrong about
the mechanism — long reads matter because they enlarge the working set, not
because of read length as such, which is why a 250 bp high-diversity pool reaches
a *higher* screen share than 1.5 kb HiFi. See
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

### PacBio HiFi — 95 samples, shared error model, full cutoff grid

The authoritative PacBio result. 95 pre-filtered samples, one `learn-errors`
model shared by every arm, `dada` per-sample at 48 threads, 3 timed replicates
interleaved across arms, timing control channel **0.2%**.

Per-sample ASV calls, against the k-mer screen at its production 0.42:

| minimizer cutoff | ASV calls | churn | abundance L1 | alignments vs k-mer |
|---|---|---|---|---|
| **control:** k-mer again | 21314 | **0** | **0.0000%** | 100.0% |
| 0.40 | 21314 | 76 | 0.0016% | **49.4%** |
| 0.42 | 21314 | 30 | 0.0003% | **59.2%** |
| **0.45** | **21314** | **0** | **0.0000%** | **77.6%** |
| 0.48 | 21314 | 0 | 0.0000% | 96.2% |
| 0.50 | 21314 | 0 | 0.0000% | 109.7% |
| 0.52 | 21314 | 0 | 0.0000% | 121.5% |
| 0.55 | 21314 | 0 | 0.0000% | 137.0% |
| 0.60 | 21314 | 0 | 0.0000% | 159.4% |

**From 0.45 through 0.60 the two screens are exactly identical** -- same ASV
calls, same counts, same 2,434,599 reads, in all 95 samples individually as well
as in total. This is a stronger concordance than any Illumina workload produced,
where no cutoff was ever exactly concordant, and it is the single most important
difference between the platforms on this page.

It also means the cutoff can be chosen on cost alone inside that window, and
**0.45 does 22% fewer alignments than the k-mer screen for bit-identical
output**. The k-mer screen at 0.42 is over-provisioned on HiFi -- the same shape
as the [band-size](band-size-platform-defaults.md) finding, arrived at
independently.

Wall clock, 3 replicates interleaved across arms, **control channel 0.2%**:

| arm | median | min-max | vs k-mer | accuracy cost |
|---|---|---|---|---|
| k-mer (baseline) | 72.94s | 72.27-73.44 | — | — |
| **control:** k-mer again | 73.08s | 72.94-73.56 | +0.2% | none |
| minimizer @ 0.42 | **48.17s** | 48.05-48.79 | **−34.0%** | 3 of 2409 pooled ASVs, −16 reads |
| **minimizer @ 0.45** | **59.29s** | 58.89-60.08 | **−18.7%** | **none — exactly identical** |
| minimizer @ 0.48 | 70.15s | 69.60-71.48 | −3.8% | none — exactly identical |
| minimizer @ 0.52 | 82.39s | — | +13.0% | none — exactly identical |

**−18.7% for bit-identical output** is the strongest result on this page. Every
other speedup here is bought with churn; this one is free, and its replicate
range (58.89-60.08) is disjoint from the baseline's (72.27-73.44).

The 0.52 row is the other half of the lesson: three grid steps looser costs
**+13.0%**, because past the concordance plateau the screen only adds
alignments. The window is asymmetric -- tightening to the plateau's cheap edge
gains 18.7%, loosening past it loses 13%, and both produce the same ASV table.

> An earlier revision of this page reported this platform from **3** Kinnex
> samples whose two arms each ran their own `learn-errors`, and drew a "13.3%
> fewer alignments" claim from that pair. Different error models make alignment
> counts incomparable, and the claim was withdrawn. Everything above holds one
> model fixed across every arm, which is now enforced structurally by
> `dev/run_screen_sweep_pacbio.sh`.

### What the PacBio churn below 0.45 actually is

Worth characterising, because it is a different failure mode from the Illumina
one and the obvious diagnosis is wrong. At cutoff 0.42, 30 of 21,314 per-sample
ASV calls differ. Every sample's ASV **count** and read total are identical, so a
set diff alone cannot say what changed. `compare_asvs.py --per-sample` pairs the
churn up inside each sample:

| at cutoff 0.42 | count | reads |
|---|---|---|
| within-sample substitutions | 14 | conserved (13 exact, 1 off by one) |
| spurious greedy pairing (H=908, not a substitution) | 1 | — |
| genuine pooled deaths | 2 | −22 |
| genuine pooled birth | 1 | +6 |

Pooled, the inventory moves from 2409 to **2408** ASVs — 3 of 2409 churn, net −1
ASV and −16 reads out of 2.43 M. The 14 substitutions are **cross-sample
reattributions**: each gained sequence is at Hamming **0** to a baseline ASV in a
*different* sample, so the sequence was already known and only the sample it is
credited to moved. Pair distances are median 3–4, not 1, so this is not
fragmentation into a near neighbour; it is a different representative chosen from
a small cloud of 3–4-mismatch variants. All 15 are `birth_type = Abundance`.

**It is not a marginal-significance effect, and OMEGA_A is not the lever.** The
natural reading -- borderline births flipping either side of the threshold -- is
falsified by the birth statistics. Reproducing one flip on a single sample
(`SRR28725000`, H=3, abundance 24 -> 24, both arms 281 calls):

| | `birth_e` | `birth_fold` | vs `OMEGA_A = 1e-40` |
|---|---|---|---|
| lost under k-mer | 1.67e-232 | 2.39e232 | ~190 orders past |
| gained under minimizer | 9.43e-222 | 4.24e221 | ~180 orders past |

Both candidates clear the threshold by ~180-190 orders of magnitude. The least
significant birth anywhere in that sample is `pval = 1.24e-42`, a factor of 80
inside `OMEGA_A`, and the two arms' birth landscapes are **identical** --
same 280 abundance births, same 209 underflows, same least-significant ASVs.
Only the identity of one center differs.

So the mechanism is the **candidate set, not the threshold**. At 0.42 the
minimizer performs 59.2% of the k-mer screen's alignments; a pair it declines
changes which uniques end up in a cluster, and therefore which unique is eligible
to be its birth center. Both candidates carry the same 24 reads, which is why the
counts are conserved exactly. Tuning `OMEGA_A` would not touch this.

> **`birth_pval` is unusable for this.** It underflows to exactly 0 for **74%**
> of ASVs in the sample (209 of 281), including both sides of the flip. Use
> `birth_e` or `birth_fold`. `--failed-uniques` is also the wrong instrument --
> its columns are `sequence / sample / reads` with no birth statistics, and the
> two arms' files are byte-identical.

### The length-variability rationale was wrong, not merely weak

I proposed length-variable amplicons (ITS2) as a promising short-read case on the
grounds that `kmer_dist8` normalises by `min(len) - k + 1` while a sketch
normalises by `min(total_a, total_b)` — "what each read actually contains". **The
distinction does not exist.** A length-L read has exactly `L - k + 1` k-mers, so

```text
min(len_a, len_b) - k + 1  ==  min(#kmers_a, #kmers_b)
```

which is the *same* min-based normalisation. Both screens are length-adaptive in
identical fashion. The real difference between them is **what is counted** — every
k-mer versus a winnowed sample — not how it is scaled. This was derivable in one
line and I did not derive it before offering it as motivation.

The measurements agree. PacBio HiFi reads are strongly variable-length and show no
directional bias between the screens:

| dataset | % reads off modal length | net / gross |
|---|---|---|
| MiSeq SOP (`truncLen` applied) | 0.00% | 0.150 |
| **ITS2 (cutadapt-trimmed)** | **6.08%** (range 1.4-50.4%) | — |
| PacBio HiFi | **85.30%** | 0.172 |

`net/gross` is indistinguishable between a fixed-length and an 85%-variable
dataset: reassignment in both, bias in neither. And ITS2 turns out **less**
length-variable by this measure than the dataset that already showed no effect.

An ITS2 run is still worth doing — as a third dataset, and to measure the
screen/align split on a high-diversity pool, where a low pass rate raises the
screen's share of runtime. It is not worth doing to test a length mechanism.

### ITS2: the screen-dominated case, and a failure mode 16S never showed

NovaSeq ITS2, 30 sample pairs, cutadapt-trimmed, per-sample `dada`, 3,808 ASVs.
Control channel bit-identical (0 of 8,395 cells); timing control 2.6%.

**On ITS2 the k-mer screen is the single largest cost in `b_compare`.** Measured
over 30 samples (`docs/findings/data/its2-phase-split.txt`):

| arm | screen share | screen ns/comp | align share | pass rate |
|---|---|---|---|---|
| **k-mer, k=5** | **43.9%** | **308 ns** | 40.6% | 1.93% |
| minimizer k=8 @0.62 | **7.5%** | **30 ns** | 63.8% | 1.92% |

A **36-point swing in screen share at a matched pass rate**, and the source of the
15% wall-clock win below.

Two effects compound to get there, and only the first was predicted:

1. **The pass rate collapses.** ITS2 is diverse enough that the screen rejects
   nearly everything, and the screen is paid on all of it while the aligner is
   paid only on survivors:

   | dataset | k-mer pass rate | measured screen share |
   |---|---|---|
   | MiSeq SOP 16S | 26.8% | 0.9% |
   | PacBio HiFi | 9.9% | 4.0% |
   | **ITS2** | **1.74%** | **43.9%** |

2. **The screen's own per-comparison cost rises 7x** — 308 ns against the 44 ns
   measured on a 1,979-unique 16S sample. Same k, same 1 KB vector. This is the
   working-set effect: ITS2 samples are individually diverse enough that their
   k-mer vectors leave cache *in per-sample mode*, without any pooling.

`share = 1 / (1 + pass x align_ns / screen_ns)` predicts 43.9% from those inputs.
An earlier revision of this page predicted ~12% from the same formula using the
16S screen cost — the formula was right and the input was wrong, and the error
was in failing to apply a cache effect this page had already identified for
pooled runs to a dataset that reaches it per-sample.

| cutoff | churn | count L1 | pass% | aligned | **wall time** | reads vs base |
|---|---|---|---|---|---|---|
| 0.40 | 58 | 5.569% | 0.96% | 57.9% | 72.0% | **-122,492 (-5.3%)** |
| 0.50 | 30 | 2.592% | 1.20% | 70.2% | 75.2% | -56,670 |
| 0.60 | 7 | 0.642% | 1.61% | 92.9% | 81.6% | -13,287 |
| **0.62** | 7 | 0.365% | **1.74%** | **100.2%** | **84.9%** | -7,145 |
| 0.65 | 6 | 0.349% | 2.02% | 115.7% | 92.6% | +3,958 |
| 0.80 | 29 | 1.525% | 8.26% | 470.4% | 231.2% | +31,925 |

**At 0.62 the two arms perform the same number of alignments (100.2%) and the
minimizer runs in 84.9% of the time.** Alignment work held constant, wall clock
down **15%** — a difference that can only be the screen, and the first clean
isolation of the screen's contribution on this branch, against a 2.6% control
channel. The phase split above shows where it comes from: 43.9% -> 7.5% of
`b_compare` busy time.

#### Read retention, not reassignment

`net/gross` on shared ASVs is **0.868** at 0.62 — systematic, where 16S and
PacBio both sat near 0.15 (cancelling). The cause is not bias in *where* reads go
but in *whether they survive*:

| cutoff range 0.40 - 0.80 | worst change in total reads |
|---|---|
| 16S, 26.8% pass | **0.073%** |
| ITS2, 1.74% pass | **5.3%** |

Read retention is ~70x more cutoff-sensitive on ITS2. Same shroud ->
`lambda = 0` chain as the [reassignment mechanism](#the-mechanism-a-mis-set-screen-moves-reads-it-does-not-lose-asvs),
but with a different endpoint: a raw shrouded from *every* cluster fails the
abundance test and is **dropped from the table entirely**, taking its reads with
it. At 26.8% pass almost every raw finds a cluster and this is negligible; at
1.74% pass a large population sits near the shroud boundary and the cutoff
directly governs how many reads survive. Retention crosses the k-mer screen's at
**~0.63**.

This is the metric to watch on any high-diversity pool, and neither ASV churn nor
count L1 exposes it — the ASV count moves only 3,800-3,832 across a range that
loses 5.3% of reads.

#### Ecological magnitude

Bray-Curtis on shared ASVs, 30 samples: median 0.0016, max **0.0043** at 0.62
(no sample above 0.01). At 0.65 the median is lower but the max is **0.0186**,
with one sample above 0.01 — the same tail-versus-centre trade seen on 16S, and
another case where the aggregate ranks the two settings the other way round.

**Recommended for ITS2-like pools: k=8, cutoff 0.62** — 15% faster at matched
alignment count, retention within 0.31% of the k-mer screen, worst-sample
Bray-Curtis 0.0043.

### The k-mer screen's cost varies 1.8x with thread placement; the minimizer's does not

Two pooled ITS2 sweeps, **same node, same binary, same 48 threads**, identical
SBATCH settings (`--cpus-per-task=48 --threads-per-core=1 --mem=256g`) except that
the second added `--exclusive`:

| | job 2367726 | job 2367778 |
|---|---|---|
| host | `compute-5-6` | `compute-5-6` (same node) |
| `--exclusive` | no | **yes** |
| CPUs allowed the process | **256 logical / 128 physical** | **96 logical / 48 physical** |
| threads | 48 | 48 |
| **k-mer screen ns/comp** | **1300** | **732** |
| k-mer screen share | 76.5% | 71.3% |
| **minimizer screen ns/comp** | **33** | **33** |
| k-mer wall | 176.6s | 119.3s |
| minimizer @0.62 | 125.5s (71.1%) | 95.3s (**79.9%**) |
| timing control | 0.5% | 5.5% |

`cpu_allocation()` reads `Cpus_allowed_list` from `/proc/self/status`, so those
are the process's **affinity mask**, not the hardware. Without `--exclusive` the
mask spanned the whole 128-core node; with it, only the 48 allocated cores.

**What is established:** at identical thread count on identical hardware, the
k-mer screen cost **1.8x more** under the wider mask, and the minimizer screen
cost the same to the nanosecond.

**What is inference:** the likely cause is thread *placement*. A 256-CPU mask lets
the scheduler spread 48 threads across many CCDs and NUMA domains of an EPYC node,
so the memory-bound `4^k` sweep crosses fabric; the narrower mask packs them. That
is [measuring on a NUMA node](measuring-on-numa.md)'s effect arriving
unintentionally, through mask width rather than deliberate pinning. **Co-tenancy
is ruled out** — no other jobs were on the node. But two runs cannot separate
placement from any other consequence of `--exclusive`, so this is the leading
explanation, not a demonstrated one. `dev/numa_pin.sh` is the tool for settling it.

Either way, the reportable consequences do not depend on which:

- **Quote ~20% for pooled ITS2**, from the `--exclusive` run. The 29% figure an
  earlier revision gave was measured against a baseline degraded by whatever this
  is, which flatters the minimizer.
- **The k-mer screen's cost is not a property of the workload alone**; it varies
  1.8x with how the scheduler places threads. The minimizer's does not. A speedup
  ratio cannot express that, and for sizing jobs on a shared cluster the
  predictability may matter as much as the mean.
- **Report the `cpu allocation` line and the exclusivity setting with any timing.**

### The limitation: the minimizer arm's compare phase is 56% serial

Pooled ITS2, 48 threads, timing control **0.5%**. Full split in
`docs/findings/data/its2-pooled-phase-split-rerun.txt`:

```text
k-mer      compare=110.76s (map= 94.01 parallel, store=15.44 serial)               -> 14% serial
minimizer  compare= 63.98s (setup=15.19 + map=28.22 parallel + store=20.57 serial) -> 56% serial
```

| arm | setup (serial) | map | map efficiency | screen share |
|---|---|---|---|---|
| k-mer | **0.0%** | 84.9% | **90%** | 76.5% |
| minimizer @0.62 | **23.7%** (15.19s) | 44.1% | **76%** | 7.8% |

Two separate problems, and the first is a design fault in this branch.

**1. `shared_counts` is serial, `O(nraw)` per cluster — 23.7% of compare.** The
inverted index moves screen work *out* of the parallel map and into
single-threaded setup. That is what makes the per-pair cost `O(1)`, and it is why
the screen drops from 76.5% to 7.8% — but the work does not vanish, it becomes
serial. An earlier revision of this page predicted ~11% at 48 threads from an
8-thread measurement; the actual is **23.7%**, so the prediction was optimistic by
2x.

Two candidate fixes have now been tested, and **both failed**.

*The buffer clear was not the cost.* `shared_counts` zeroes the whole `nraw`
array per cluster, and tracking touched indices instead looked obvious. Built, it
was a **regression**: 5.97 -> 7.70s of setup on PacBio, +29%, because it adds a
branch and a push per increment. The arithmetic never supported the clear anyway
-- a 3.3 MB memset is ~165 us, ~0.6s across 3,414 clusters, not the 15s being
explained. The **scatter** is the cost, not the clear.

*Removing the index does not help either — it costs more.* Without the index each
pair does its own merge-join over two sketches, which is fully parallel and has
no serial setup at all. Measured as a first-class arm on per-sample PacBio
(`DADA2RS_MINIMIZER_INDEX=0`, 3 reps, control 0.2%):

| arm | screen | ns/comp | screen share | `setup` (serial) | wall |
|---|---|---|---|---|---|
| k-mer | 322.11s | 2408 | 12.3% | 0.00s | 72.94s |
| minimizer @0.45, index **on** | 4.29s | **32** | 0.2% | 33.44s | **59.29s** |
| minimizer @0.45, index **off** | 321.92s | **2407** | 15.6% | 0.00s | 62.00s |
| minimizer @0.42, index **on** | 4.28s | **32** | 0.3% | 31.13s | **48.17s** |
| minimizer @0.42, index **off** | 303.98s | **2272** | 19.1% | 0.00s | 51.22s |

The merge-join costs **2272-2407 ns/comp against the k-mer sweep's 2408** — the
same, within noise. Both screens are **memory-latency-bound, not op-count-bound**: two
~475-entry sketches fetched from scattered memory cost what sweeping a `4^k`
vector costs, however few operations you perform once the lines are resident.
This is the third negative result on this branch from costing a scan by
operations rather than by access pattern.

So ~31-33s of serial setup **buys back ~300s of parallel screen work**, and
index-off is slower end to end at every cutoff tested — by 4.4% at 0.45 and
**6.0% at 0.42**, both with disjoint replicate ranges. The 2.0% at 0.48 is
**not a result**: it is below the minimizer arm's own noise floor (see below). The
index earns its serial cost, and `DADA2RS_MINIMIZER_INDEX` is a diagnostic, not a
tuning knob.

**The margin growing as the cutoff tightens is the mechanism confirming
itself.** Tightening removes alignment work, so the screen's share of
`b_compare` rises — 12.3% for the k-mer screen, 15.6% for index-off at 0.45,
**19.1% at 0.42** — and a cheaper screen is worth proportionally more. An index
whose serial cost were the dominant term would trend the other way.

**Pooled ITS2 settles it the other way, by a wider margin.** This was predicted
to be the case that *favoured* index-off — sketches are ~4x smaller on 300 bp
reads, and a single pool at 48 threads has no cross-sample concurrency to hide
serial time behind. Measured (3 reps, control 0.5%):

| arm | screen ns/comp | screen share | `setup` | compare | wall |
|---|---|---|---|---|---|
| k-mer | 1383 | 77.2% | 0.00s | 116.20s | 176.57s |
| minimizer @0.63, index **on** | **33** | 9.1% | 16.14s | **60.04s** | **126.59s** |
| minimizer @0.63, index **off** | 1144 | 79.5% | 0.00s | 96.11s | 146.99s |

**The index wins by 13.9%** (126.59s vs 146.99s, replicate ranges 125.99-127.10
and 145.56-148.23, disjoint).

The half of the prediction about sketch size was right and irrelevant. The
merge-join *is* much cheaper here than on HiFi — **1144 ns/comp against 2272-2407**
— and it even beats the k-mer screen it replaces (1383). But that was never the
comparison that decides anything. **Index-off is 35x more expensive per
comparison than index-on**, and on a workload where the screen is 77% of
`b_compare`, 16s of serial setup against ~2650s of parallel screen work is not a
close call. Reasoning about the merge-join versus the *k-mer screen* instead of
versus the *index* is what made this look open.

So the serial scatter is not a design fault to be parallelised away on the
strength of its share alone. An earlier revision of this page said it was, and
projected the arm from 71% to ~60% of k-mer wall clock on that basis. That
projection is withdrawn.

> **REOPENED by pooled PacBio, which reverses the sign.** 3 reps, control
> **0.03%**:
>
> | arm | median | range | vs k-mer |
> |---|---|---|---|
> | k-mer | 335.69s | 335.53-341.69 | — |
> | minimizer @0.45, index **on** | 423.23s | 420.20-435.89 | **+26.1%** |
> | minimizer @0.45, index **off** | **289.30s** | 285.98-294.92 | **−13.8%** |
>
> Index-off beats index-on by **31.6%**, and the indexed arm is *slower than the
> k-mer screen it replaces* — the only configuration where that happens.

**The serial scatter is confirmed as the cause.** Phase split, index-on vs the
two arms that have no setup at all:

| arm | screen ns/comp | align | map | **`setup`** | compare |
|---|---|---|---|---|---|
| k-mer | 3000 | 9926.95s | 263.28s | 0.00s | 272.88s |
| minimizer, index **on** | **32** | 7316.50s | **169.74s** | **162.69s (46.9%)** | 347.09s |
| minimizer, index **off** | 3440 | 7251.39s | 213.46s | 0.00s | 223.18s |

The index does exactly what it is supposed to inside the parallel region — map
drops 213.46s -> 169.74s — and then pays **162.69s of serial time**, 46.9% of the
compare phase, for a 43.72s parallel saving. That is also why core usage visibly
collapses on this arm.

**Two things the data falsified about my own prediction.** I predicted mean
sharing would reach 10^4-10^6 because conserved 16S saturates a small minimizer
alphabet. Measured: **6,509x**, and the alphabet did *not* saturate — distinct
minimizers grew from ~10.5k per sample to **40,177** in the pool. The
conservation story was wrong. What survives is the cost ordering:

| config | entries x sharing | `setup` | map saved | **setup / saved** | index verdict |
|---|---|---|---|---|---|
| pooled ITS2 | 93,092 | 16.14s | 56.14s | **0.29** | wins 13.9% |
| per-sample PacBio | 207,240 | 33.44s | 76.26s | **0.44** | wins 4.4% |
| **pooled PacBio** | **3,111,302** | **162.69s** | 43.72s | **3.72** | **loses 31.6%** |

`entries x sharing` ranks all three correctly and predicts the total scatter work
within ~3x. The decision quantity is `setup / map-saved` crossing 1, and a proxy
computable **before any scatter runs** — `sharing x threads / nraw`, from
`n_postings/n_keys` at index-build time — is monotone in the same order (0.073,
0.185, 0.571). Three points bracket the crossover between **0.19 and 0.57**;
that is a weakly determined threshold, not a calibrated one.

**Concordance is unaffected by any of this.** Under pooling, cutoffs 0.45, 0.48
and 0.50 are identical to the k-mer screen in **all 95 samples** on ASV calls and
reads, index-on and index-off alike; 0.42 differs in 9 samples. So the exactly-
concordant plateau at >=0.45 holds in pooled mode as well as per-sample, and the
index question is purely about speed.

### The index is now chosen per workload, not always built

Four configurations, and the index wins three and loses one badly — so "always
index" is wrong and so is "never". `minimizers::decide_index` scores each pool
and picks:

```text
score = sharing × threads / nraw          sharing = entries / distinct
```

The derivation: serial scatter costs `nclusters × entries × sharing`, the
parallel saving it buys is `ncomps × entries × c / threads`, and
`ncomps ≈ nclusters × nraw`, so `nclusters` and `entries` cancel. `threads` is
`rayon::current_num_threads()` **inside** `B::new`, which is the sub-pool (4
under per-sample concurrency, 48 pooled) — the global count would inflate the
score ~12x and decline the index exactly where it wins.

| configuration | score | measured setup/saved | index |
|---|---|---|---|
| pooled ITS2 | 0.073 | 0.29 | wins 13.9% |
| per-sample PacBio | 0.185 | 0.44 | wins 4.4% |
| pooled PacBio | **0.564** | 3.72 | **loses 31.6%** |

Default threshold **0.30**, overridable with `DADA2RS_MINIMIZER_INDEX_MAX`.
**Three points determine a bracket, not a threshold** — the crossover lies in
0.19..0.57 — and 0.30 is deliberately biased low because the errors are not
symmetric: indexing wrongly cost **+26%** against the k-mer screen, declining
wrongly costs 2-4%.

**Validated on pooled PacBio**, the case it exists to catch. 3 reps, control
**0.03%**, with both forced modes run alongside:

| arm | median | range | vs k-mer | `setup` |
|---|---|---|---|---|
| k-mer | 335.69s | 335.53-341.69 | — | 0.00s |
| forced **on** | 423.23s | 420.20-435.89 | +26.1% | 160.09s |
| forced **off** | 289.30s | 285.98-294.92 | −13.8% | 0.00s |
| **auto (rule)** | **290.35s** | 288.38-297.82 | **−13.5%** | **0.00s** |

```text
[dada] minimizer index: declined (per-pair merge-join) — score 0.564
  (= sharing 6430.3× × 48 threads / 547273 raws) vs threshold 0.300;
  40177 distinct minimizers, 258350362 postings
```

The rule lands **+0.36% from the right choice** (overlapping ranges) and −31.4%
from the wrong one, and `setup` is 0.00s on the auto arm — so the decision is
real, not a fast run that happened to coincide.

**The auto arm is a second control channel, and it reads higher than the first.**
`decide_index` runs whenever the sketches are present, *including* when the mode
is forced — so `_auto` and a forced-on arm execute identical work: same
distinct-count pass, same index build, same scatter. Any gap between them is
noise **by construction**. On pooled soil 16S that gap is **1.65%** (auto 612.06s
vs forced-on 622.34s, overlapping ranges), against a k-mer control channel of
0.75% on the same run.

That is worth more than it looks. Every control channel on this page so far has
been k-mer-vs-k-mer, which measures the *baseline* arm's variance and not the
arm the claims are about. The minimizer arm's own floor is roughly **double**
it — which is why the 2.0% index margin at PacBio cutoff 0.48 is withdrawn
above, and why the 15-28% Illumina figures are unaffected.

**The decision costs at most ~1s of 290s.** Counting 40,177 distinct minimizers
over 258 M postings without building them was the part of this most likely to be
too expensive to be worth it. The cached `_noidx` replicates predate the rule, so
`auto − noidx = +1.05s` is an upper bound — inside the noise, and it also
absorbs any other difference between the two binaries.

`DADA2RS_MINIMIZER_INDEX` remains an **explicit override** (0/false/off/no,
1/true/on/yes, or auto), and warns rather than silently defaulting on an
unrecognised value — it used to be that any value but `0` meant on, so
`=true` would otherwise have changed meaning when the rule landed.

**Still unvalidated outside the fitted points.** Pooled soil 16S is the
discriminating case: a large pool of *conserved* sequence like pooled PacBio, but
with short reads and small sketches like ITS2, so the two terms pull opposite
ways. If the auto arm there matches the slower forced arm with a score near 0.30,
the threshold is miscalibrated and tunable; if it picks wrong with a score far
from 0.30, the score is the wrong quantity.

**The minimizer's pooled-PacBio win is not a cheaper screen.** At 3440 ns/comp
the merge-join is *more expensive* than the k-mer sweep's 3000 here. Index-off
still finishes 13.8% ahead purely because 0.45 passes 9.00% of pairs against the
k-mer screen's 12.10% — a 27% cut in alignment work. On this workload the
minimizer is a better *filter*, not a cheaper *screen*.

### Calibrating on read retention: the cutoff is ~0.64

Read retention is the best-behaved signal for choosing a cutoff on
screen-dominated data: unlike ASV churn (discrete, a handful of events) or count
L1 (a flat-bottomed U), it is smooth, monotone in cutoff, and has a **true zero
crossing** — the point where the minimizer arm recovers the same number of reads
as the k-mer screen. Interpolating it:

| dataset | retention-neutral cutoff | matched-pass cutoff |
|---|---|---|
| pooled 16S | **0.637** | 0.65 |
| pooled ITS2 | **0.636** | 0.62 |
| per-sample ITS2 | **0.639** | 0.62 |
| per-sample soil 16S | 0.663 | 0.65 |
| MiSeq SOP 16S | none — retention flat (<0.02%) | ~0.80 |

**Three of four land at 0.636-0.639**, across two amplicons and both pooling
modes — a much tighter target than matched-pass-rate, which spans 0.50 (PacBio)
to 0.80 (MiSeq SOP). So **~0.64 is a reasonable single default for
screen-dominated Illumina data**, with a per-dataset sweep to confirm.

The two targets are systematically offset, retention-neutral sitting **~0.02
above** matched-pass in every case. That follows from the screens being different
*sets* rather than different *rates*: matching the pass rate matches alignment
count, but the minimizer's misses are not the k-mer screen's misses, so it drops
slightly more reads and needs a little more permissiveness to break even on
recovery. Which of the two to target is a real choice — matched-pass maximises
speed at equal work, retention-neutral preserves read recovery — and they differ
by about one grid step.

MiSeq SOP has no crossing at all: its retention moves <0.02% across the whole
0.40-0.80 range, because at 26.8% pass rate almost every raw finds a cluster. The
signal exists only where the screen dominates, which is where the cutoff matters.

### Soil 16S: the same effect, from diversity alone

NovaSeq soil 16S, 30 samples, per-sample `dada`, **17,398 ASVs** — the most
diverse table measured here. Control bit-identical (0 of 45,934 cells); timing
control 3.0%. Split in `docs/findings/data/soil16s-phase-split.txt`.

| arm | screen share | screen ns/comp | pass rate |
|---|---|---|---|
| **k-mer, k=5** | **29.8%** | **454 ns** | 5.52% |
| minimizer k=8 @0.65 | **2.7%** | **30 ns** | 5.64% |

This is the control for the diversity hypothesis: same instrument and prep as
ITS2, a *16S* amplicon, and the k-mer screen still costs **29.8%** of `b_compare`
against 0.9% on the mouse-gut MiSeq SOP. So the driver is not the ITS2 amplicon —
it is **how diverse the pool is**, which sets both the pass rate and the working
set. Note soil 16S has a *higher* per-comparison cost than per-sample ITS2 (454 vs
308 ns, a larger working set) but also a higher pass rate, so a lower share: both
terms in `share = 1/(1 + pass x align_ns/screen_ns)` matter.

| cutoff | churn | count L1 | aligned | reads vs base | wall time |
|---|---|---|---|---|---|
| 0.40 | 3816 | 14.677% | 21.5% | -67,726 | 34.4% |
| 0.55 | 949 | 3.818% | 52.3% | -13,110 | 52.9% |
| 0.62 | 273 | 0.942% | 82.6% | -3,953 | 71.2% |
| **0.65** | 277 | **0.842%** | **102.4%** | **-911** | **83.4%** |
| 0.70 | 483 | 1.703% | 155.5% | +2,476 | 114.5% |
| 0.80 | 730 | 2.762% | 483.8% | +2,743 | 312.2% |

**16.6% faster at matched alignment work**, with retention within 0.06% and count
L1 at its minimum. Churn is 277 of 17,398 ASVs (1.6%) — larger in absolute terms
than any other dataset here, and proportionally larger than ITS2's 0.6%, which is
what a table with 17k mostly-rare ASVs should do.

### ITS2 under full pooling: the screen becomes 76.5% of the phase

Same 30 ITS2 samples, `dada-pooled` (R `pool=TRUE`), 3,028 ASVs, 2.39 **billion**
comparisons. Control bit-identical (0 of 13,651 cells); timing control 2.8%.
Raw split in `docs/findings/data/its2-pooled-phase-split.txt`.

| arm | screen share | screen ns/comp | align share | pass rate |
|---|---|---|---|---|
| **k-mer, k=5** | **76.5%** | **1305 ns** | 6.7% | 0.70% |
| minimizer k=8 @0.62 | **9.3%** | **33 ns** | 34.0% | 0.70% |

Pooling raises the k-mer screen from 43.9% to **76.5%** of busy time and its
per-comparison cost from 308 to **1305 ns**, while the minimizer's stays at 33.
Nothing about the algorithm changed — only how much memory the screen structures
occupy, which is the point.

| cutoff | churn | count L1 | aligned | reads vs base | wall time |
|---|---|---|---|---|---|
| 0.40 | 91 | 11.829% | 30.2% | -260,246 | 94.4% |
| 0.50 | 56 | 6.090% | 49.4% | -133,444 | 81.2% |
| 0.60 | 23 | 1.351% | 88.5% | -28,930 | 73.7% |
| **0.62** | 18 | 0.769% | **101.3%** | -15,830 | **69.3%** |
| 0.65 | 14 | 0.810% | 128.2% | +14,024 | 71.6% |
| 0.80 | 53 | 3.669% | 818.4% | +78,855 | 82.6% |

**At 0.62 alignment work is matched (101.3%) and wall clock is 69.3% — a 30.7%
speedup attributable to the screen alone**, against a 2.8% control channel, and
double the 15% seen per-sample.

**Wall time is non-monotonic**, which per-sample runs did not show: 0.40 is
*slower* (94.4%) than 0.62 despite performing only 30.2% of the alignments,
because it yields 3,069 ASVs against 3,026 and more clusters means more
`b_compare` invocations and more shuffle/bud work. Under pooling, over-tight
screening costs more in cluster proliferation than it saves in alignment — so the
cost curve has a genuine interior optimum rather than trading monotonically
against accuracy.

Read retention behaves as per-sample, scaled up: -260k reads (-11.6%) at 0.40,
crossing the k-mer screen's retention near 0.63.

### Is the count matrix *distorted*, or just perturbed?

Count L1 is a global sum, and three different situations produce the same number:
a small error spread over every cell, a few samples wrecked while the rest are
exact, or a systematic bias that shifts every composition. `dev/compare_counts_correlation.py`
separates them. Against the k-mer baseline on 362 samples (control arm reads
1.000000 / 100% exact, as it must):

| | k=8 @ 0.62 | k=8 @ 0.70 |
|---|---|---|
| Pearson r, **raw** counts | 0.999976 | 0.999969 |
| Pearson r, log10 counts | 0.995613 | 0.995445 |
| Spearman rho | 0.998261 | 0.997911 |
| cells exactly equal | 98.578% | 98.836% |
| cells within 10% | 99.615% | 99.672% |
| **net / gross change** | **0.150** | **0.107** |
| Bray-Curtis median | 0.000070 | 0.000053 |
| **Bray-Curtis max** | **0.0092** | **0.0267** |
| **samples with BC > 0.01** | **0** | **3** |

**Raw-count Pearson is not diagnostic** — 0.99998 for an arm with 19 churned
ASVs and 0.5% of reads misplaced, because abundant ASVs pin it. Any correlation
check here has to be on log counts or ranks.

**The disagreement is reassignment, not bias.** `net/gross` of 0.10-0.15 means the
changes very nearly cancel: reads move *between* ASVs rather than appearing or
vanishing. That is the
[shuffle mechanism](#the-mechanism-a-mis-set-screen-moves-reads-it-does-not-lose-asvs)
confirmed from the count matrix, independently of the `map` diff that first found
it.

**And a result L1 concealed: 0.70 has a worse tail than 0.62.** Their L1s are
nearly equal (0.0653% vs 0.0612%), but 0.70 — the *zero-churn* setting —
concentrates its error into three samples, one at Bray-Curtis 0.027, while 0.62
never exceeds 0.0092 and leaves no sample above 0.01. **0.62 is better on
worst-case per-sample dissimilarity and 5% faster; its cost is 5 tail ASVs.** An
aggregate metric ranked these the other way round.

Error by ASV abundance decile at 0.62 confirms the shape is benign: the most
abundant decile holds 84% of all reads and moves **0.0421%**, while relative
deviation rises to ~2% only in decile 9 (644 reads total), which is counting
noise at that magnitude.

### Illumina MiSeq — the full cutoff curve, 362 samples

The complete 0.40-0.80 sweep at k=8, against a bit-identical control (0 of 30,105
count cells) and a **2.4% timing control channel**. Timing arms all share the
baseline error model, so the wall-clock column isolates the screen.

| cutoff | ASVs | churn | count L1 | max churned abundance | wall time vs k-mer |
|---|---|---|---|---|---|
| 0.40 | 415 | 19 | 0.4986% | 15 | **42.8%** |
| 0.45 | 413 | 11 | 0.3140% | 13 | 52.3% |
| 0.50 | 413 | 9 | 0.1864% | 13 | 65.8% |
| 0.55 | 414 | 6 | 0.1132% | 13 | 78.7% |
| 0.60 | 416 | 6 | 0.0706% | 77 | 89.6% |
| **0.62** | 413 | 5 | **0.0612%** | 13 | **94.8%** |
| 0.65 | 415 | 3 | 0.0660% | 9 | 103.8% |
| **0.70** | **416** | **0** | 0.0653% | — | 116.9% |
| 0.75 | 416 | 2 | 0.0700% | 78 | 133.6% |
| 0.80 | 417 | 1 | 0.0843% | 78 | 149.6% |

Four things fall out of this table.

**1. `n_asv` is worse than useless here.** ASV *count* spans 413-417 across the
whole range while count L1 varies **8x**. At cutoff 0.40 the table holds 415 ASVs
against the baseline's 416 — by ASV count, indistinguishable — while 19 ASVs have
churned and 0.5% of reads sit on the wrong variant. Any sweep judged on ASV counts
would have called 0.40 fine.

**2. A mis-set screen never loses abundant ASVs.** The maximum churned abundance
is **15 reads at cutoff 0.40**, where the sketch retains only ~53% of
near-neighbours. Damage is confined to the low-abundance tail at every setting;
what changes is *where the reads go*, per
[the mechanism section](#the-mechanism-a-mis-set-screen-moves-reads-it-does-not-lose-asvs).

**3. Accuracy is U-shaped, not monotone.** L1 bottoms out at 0.62-0.70 and rises
again to 0.0843% at 0.80. Above 0.70 the minimizer is *more* permissive than
k-mer@0.42 and starts aligning pairs k-mer rejects, diverging in the opposite
direction. **The optimum is where the two screens' effective behaviour matches,
not where recall is maximal** — which is the same lesson as
[calibration bracketing](#calibration-brackets-only-a-sweep-decides), arrived at
from the output side.

**4. Illumina does have a real operating point, and it is not free.** Wall time is
monotone in cutoff and far outside the 2.4% control:

- **0.62: 5% faster, at the L1 minimum, for 5 churned ASVs of <=13 reads.**
- 0.60: 10% faster, 6 churned ASVs.
- 0.70: zero churn, but **17% slower**.

An earlier revision of this page called Illumina "accuracy-neutral, cost-negative".
That was wrong — it was inferred from a 20-sample subset without timing. The real
trade is explicit: **a handful of tail ASVs for 5-10% wall clock**, or exact
reproduction for +17%.

### Illumina MiSeq — full run, 362 samples, 416 ASVs, 2.97 M reads

The `(k, cutoff)` grid, against a control channel (k-mer run twice) that is
bit-identical: **0 of 30,105 count cells**.

| k=8, cutoff | ASV churn | count L1 | aligned vs k-mer |
|---|---|---|---|
| 0.55 | 6 | 0.1132% | 75.4% |
| 0.58 | 4 | 0.0938% | 84.0% |
| 0.60 | 6 | 0.0706% | 90.3% |
| 0.62 | 5 | 0.0612% | **95.7%** |
| 0.65 | 3 | 0.0660% | 104.6% |
| **0.70** | **0** | 0.0653% | **121.3%** |
| 0.72 | 1 | 0.0730% | 128.4% |

k=9 is worse than k=8 on set agreement at every cutoff and is not recommended.

**There is no cutoff on Illumina that buys both.** Zero ASV churn requires 0.70,
which costs **+21% alignment work**; cost parity sits near 0.63, which churns 3-5
ASVs of 416. Compare PacBio, where zero churn arrives *with* 7-11% **fewer**
alignments. That is the platform split in one line: **long reads get accuracy and
speed together, short reads must trade one against the other.**

It also corrects a recommendation made from the 20-sample SOP subset, which
suggested 0.58-0.65. At 362 samples the answer is 0.70: more samples means more
marginal ASVs with a chance to churn, so more recall is needed. The subset
under-estimated the cutoff, which is the
[fixture lesson](#how-the-fixtures-misled) one scale up.

> **The alignment column carries a confound.** 14 of the 16 arms ran their own
> `learn-errors`, so they do not share an error model, and the screen shapes the
> model. The monotone trend with cutoff is plainly cutoff-driven, but the exact
> percentages would move under a shared model — which is what
> `dev/run_screen_sweep_pacbio.sh` enforces structurally and what an `ERR_DIR=`
> re-run would give for Illumina. The ASV columns are unaffected: those are
> legitimate end-to-end comparisons of complete configurations.

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

## Cost: the screen is 33x cheaper, and that is not where the win is

`--verbose` splits `b_compare`'s busy time (#127 instrumentation). PacBio HiFi,
one 224k-read sample, k=7, **shared error model** so only the screen differs:

| | k-mer | minimizer k=8 @0.50 |
|---|---|---|
| screen | 38.41 s (4.0%), **660 ns/comp** | 1.15 s (0.1%), **20 ns/comp** |
| align total | 897.1 s (93.2%) | 989.7 s (97.0%) |
| aligned pairs | 5,764,121 (9.9% passed) | 6,315,520 (10.9% passed) |

The sketch plus inverted index make the screen **33x cheaper per comparison** —
the index turns the per-pair merge-join into an array lookup. But the k-mer screen
is only **4.0%** of busy time, so making it free caps out at 4%. **The lever is
alignment count**, at 93-97%, and that is set purely by the cutoff.

At cutoff 0.50 the minimizer aligns ~10% *more* pairs and loses. At 0.42-0.45 it
aligns 7-11% fewer:

| arm (shared error model) | aligned | vs k-mer | ASVs | ASV set |
|---|---|---|---|---|
| k-mer k=7 @0.42 | 5,764,121 | 100.0% | 899 | — |
| minimizer @0.45 | 5,337,487 | **92.6%** | 899 | identical, abundance L1 = 2 reads |
| minimizer @0.42 | 5,136,392 | **89.1%** | 899 | identical, abundance L1 = 2 reads |

**The last 1.5% of recall bought nothing and cost ~20% of the alignment work.**
899 ASVs at 0.42, 0.45 and 0.50 alike. This is why the calibration must *bracket*
the region and cost must choose within it — optimising recall alone lands at the
expensive end every time, and did here.

### Wall clock: not resolvable on the development machine

| arm | n | median s | spread | vs k-mer |
|---|---|---|---|---|
| k-mer | 4 | 125.22 | 9.0% | 100.0% |
| **control:** k-mer again | 4 | 129.79 | 5.3% | **103.6%** |
| minimizer @0.45 | 4 | 123.59 | 8.7% | 98.7% |

**Control channel 3.6%; measured effect 1.3%. Inside the noise floor, therefore
not a result.** From the split above the mechanism predicts ~11% (4.0% screen +
7.4% x 93.2% alignment), and the deterministic counts support it — but this box
cannot distinguish it, and a 2-replicate run earlier *did* produce a "~10% faster"
figure that four replicates erased. Same failure as
[measuring on a NUMA node](measuring-on-numa.md), where a genuine 3.5-6% gain was
twice reported as flat: **a speed claim is a claim about the rig first.**

What is safe to state: the screen is 33x cheaper per comparison, the calibrated
cutoff performs 7-11% fewer alignments, the ASV set is bit-identical, and the
resident screen is 74% smaller. Whether that adds up to wall-clock time needs a
quiet, pinned node.

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

Holding the error model fixed and varying one factor at a time decomposes it.
**Measured after the gapless fix** — an earlier revision reported this pre-fix,
where ~80% of the "denoising" term was actually the gapless method-selection bug
rather than screening (20-sample MiSeq SOP, k=8, cutoff 0.65):

| arm | ASV churn | count L1 | share | `err_out` vs baseline |
|---|---|---|---|---|
| k-mer model + **minimizer denoising** | 2 | 0.0555% | **85%** | 1.00x (uses baseline model) |
| **minimizer model** + k-mer denoising | **0** | 0.0105% | 16% | 1.05x |
| both changed | 2 | 0.0652% | — | 1.05x |

**The denoising-stage screen dominates ~5:1.** A 5% error-model difference yields
*zero* ASV churn and 0.0105% L1; swapping the denoising screen yields churn 2 and
0.0555%. So **whichever screen you denoise with determines the answer, almost
regardless of which screen trained the model.**

The gapless correction did move the balance: pre-fix the split read 94%/8%,
post-fix it is 85%/16%, so the model's relative contribution doubled even as the
absolute numbers shrank 5x. Worth noting because it is the direction that would,
if it continued, eventually make the model the dominant term — and this is
low-diversity data, where the screen is 0.9% of runtime and shapes little of what
reaches `build_trans_mat`. On a pooled ITS2 pool the screen is 76.5% and the model
difference should be much larger, so the decomposition is worth repeating there
before this is treated as general.

Two further points worth keeping:

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

## Calibration brackets; only a sweep decides

Full 362-sample MiSeq, **169.2 M sampled pairs per curve**, `--min-core 200`
(16,699 pairs dropped, 0.01%):

| | k-mer | minimizer k=8 |
|---|---|---|
| p99 / p99.9 near-neighbour distance | 0.4785 / 0.5255 | 0.5730 / 0.6180 |
| recall at 0.42 | **88.45%** | 61.06% |
| recall at 0.50 | 99.68% | 89.88% |
| recall at 0.60 | 100.00% | 99.73% |
| cutoff matching k-mer@0.42 | — | **0.50 by recall AND by pass-rate** |

Two results here, and the second is the one that matters.

### Why cutoff decoupling helps here for a different reason

[KDIST cutoff decoupling](kdist-cutoff-decoupling.md) found that keeping
`learn-errors` lenient while tightening `dada` was safe and fast. That works
because it is **one metric at two thresholds**: the strict set is a *nested
subset* of the lenient one, so the extra pairs denoising drops are by
construction pairs the lenient learn phase already saw and judged irrelevant.
That nesting is what makes "headroom" a meaningful idea.

Across two *mechanisms* the nesting breaks — but only because the calibrated
cutoffs differ. The [pair-level audit](#the-pair-level-audit) found the minimizer
screen is a strict subset of the k-mer screen **at equal cutoff**
(minimizer-only = 0 in every configuration tested). Running them at their own
operating points — 0.64 against 0.42 — makes them two subsets of *different*
parents: overlapping, neither containing the other. Hence
["a different 11%"](#matching-recall-does-not-reproduce-the-table), and hence no
headroom argument. A lenient learn cutoff cannot recover the k-mer answer, because
the difference lives in the denoising *mechanism*, not in its threshold.

**But there is a distinct reason to decouple here.** As the screen opens, both
mechanisms converge on the same thing — with both screens fully open they are
**bit-identical** ([control](#the-fix-and-what-it-says-about-the-k-mer-screen)),
because an unscreened comparison has no mechanism. So a *lenient* learn cutoff
makes the fitted error model approach mechanism-neutral, removing the model term
from the difference and leaving only the denoising-stage effect. That is worth
having for two reasons: the resulting model is directly comparable to a
k-mer-trained one, and it isolates the screen for measurement without the
`ERR_DIR` gymnastics used above.

The payoff is bounded — the model term is only ~16% of the difference on the SOP —
so this buys interpretability more than accuracy. `LEARN_CUTOFF` in
`dev/run_screen_sweep.sh` exposes it.

### DADA2's own 0.42 is not lossless at scale

**The k-mer screen at `KDIST_CUTOFF = 0.42` retains only 88.45%** of pairs within
10% true divergence, and does not reach 100% until ~0.60. On a 200k-pair sample it
measured as exactly 100% — 200k pairs simply do not reach the tail. This is a
finding about the *production* screen, independent of minimizers, and it sits
beside the other one on this page: 0.42 is also not lossless against
[unscreened denoising](#which-arm-is-actually-more-faithful).

### Matching recall does not reproduce the table

Both calibration criteria agree the minimizer should use **0.50** to match
k-mer@0.42. The sweep says **0.70** for zero ASV churn; 0.55-0.60 churn 6.

The calibration is not wrong — it answers a different question. **Aggregate recall
matching is not set matching.** The pair-level audit established the minimizer
screen is a strict *subset* of the k-mer screen at equal cutoff, so at 0.50 it
rejects ~11% of near-neighbours, but a **different** 11% than k-mer rejects at
0.42. Reproducing k-mer's ASV output requires never dropping a pair k-mer kept —
which means being strictly **more permissive than k-mer**, not equally permissive.
Hence 0.70, where minimizer recall reaches 100% against k-mer's 88.45%.

**So the rule is:** calibrate to bracket the region, then sweep against the actual
table, and expect the sweep to land *above* what recall-matching predicts. Two
screens with the same aggregate recall are not interchangeable; only a superset
is safe.

## The mechanism: a mis-set screen moves reads, it does not lose ASVs

A too-tight screen degrades **abundance**, not richness, and the stage
responsible is **`b_shuffle`'s raw-to-cluster assignment**, not budding.

20-sample MiSeq SOP, shared error model, k=8:

| cutoff | uniques reassigned to a different ASV | Abundance births (base -> arm) | churn | count L1 |
|---|---|---|---|---|
| 0.42 (too tight, ~60% recall) | **201 of 35,722 (0.563%)** | 2187 -> 2198 (**+11**) | 13 | **0.4257%** |
| 0.80 (too loose) | 8 (0.022%) | 2187 -> 2185 (-2) | 1 | 0.0201% |

Reassignment outnumbers budding **18:1**. The causal chain:

> shrouded pair -> `sub = None` -> `lambda = 0` -> the raw's best cluster is
> invisible to it -> it joins its **second-best** -> its reads land on a
> different ASV, and **both ASVs survive**.

The set survives because reassignment moves a raw between clusters that both
already exist and both keep other members. A cluster only disappears if it loses
*every* member, which is why the churn that does occur is confined to 2-7-read
ASVs — few enough members that losing a couple empties them. At cutoff 0.42, 12
of 13 churned ASVs hold 2-7 reads and only one exceeds 13.

Budding does respond in the expected direction, just weakly. It takes over only
at far worse recall: the original k=11 default shrouded raws from *every* cluster,
so `lambda` collapsed everywhere, abundance p-values went to zero and clusters
fragmented — which is when the ASV *set* genuinely broke (241 vs 232).

**Consequence for how these sweeps are read:** a screen sweep judged on ASV counts
would rate cutoff 0.42 as nearly fine — 13 churned ASVs of 232, all but one in the
low-abundance tail — while reads are being reassigned at 7x the rate of the good
range. Only the count matrix sees it. This is the same abundance-not-richness
shape as the [NovaSeq binned-quality result](binned-quality-illumina-novaseq.md),
arrived at through a completely different mechanism, and it is why
`dev/compare_seqtab_matrix.py` exists alongside `compare_asvs.py`.

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

## Three claims this falsified

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


### 3. "The serial scatter is the index's price and should be parallelised away." It is the index's bargain.

The inverted index makes the screen `O(1)` per pair by scattering posting lists
into an `nraw` array once per cluster — serially. At 23.7% of `b_compare` on
pooled ITS2 that looked like an obvious design fault, and this page said so.

Both ways out were measured and both cost more:

* **Tracking touched indices instead of clearing the array** was a 29%
  regression (5.97 -> 7.70s), and the arithmetic never supported the clear as the
  cost in the first place — a 3.3 MB memset is ~165 us.
* **Removing the index entirely**, so each pair does its own fully-parallel
  merge-join with no serial setup, measured **2272-2407 ns/comp against the
  k-mer sweep's 2408** — no better than the screen it replaces, and slower end to
  end at every cutoff tested (control 0.2%, disjoint ranges at 0.42 and 0.45):

| cutoff | index on | index off | index better by |
|---|---|---|---|
| 0.42 | **48.17s** | 51.22s | **6.0%** |
| 0.45 | **59.29s** | 62.00s | **4.4%** |
| 0.48 | 70.15s | 71.60s | 2.0% — *below the noise floor, not a result* |

Both screens are **memory-latency-bound, not op-count-bound**. Two ~475-entry
sketches fetched from scattered memory cost what a `4^k` sweep costs, however few
operations run once the lines are resident. So ~31-33s of serial setup buys back
~300s of parallel screen work, and the index's margin *grows* as the screen
tightens — because tightening removes alignment work and raises the screen's
share (12.3% k-mer, 19.1% index-off at 0.42). An index whose serial cost
dominated would trend the other way.

This is the third negative result on this branch from costing a scan by
operations rather than by access pattern, which is a general lesson about this
codebase and not about minimizers.

**Not closed — pooled PacBio reverses it.** Pooled ITS2 was expected to favour
index-off and the index won there by **13.9%**; pooled PacBio then went the other
way by **32.8%** (index-off 287.64s vs index-on 428.04s, 2 reps, control 1.0%).
The answer is workload-dependent, and the working model is that the scatter costs
`nraw x sketch entries` — a product only pooled PacBio maximises. The merge-join is indeed
cheaper on short reads (1144 ns/comp vs 2272-2407) and even beats the k-mer
screen, but it is still **35x** the index's 33 ns/comp, and the screen is 77% of
`b_compare` on that workload.
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

## Status: experimental, and what promotion would require

**This stays experimental and opt-in.** `--screen-backend minimizer` is off by
default, on the same footing as the [WFA alignment backend](https://github.com/HPCBio/dada2-rs/issues/51),
and nothing here argues for flipping the default. The evidence base is five
workloads, each against a bit-identical control, but every one of them is a single
dataset per configuration.

What has been settled:

- ✅ Per-backend cutoff calibration exists (`kdist-calibrate --screen-backend`)
  and its recommendation is derived, not guessed.
- ✅ The screen is recorded in the error model's `params` and checked on mismatch.
- ✅ The gapless method-selection coupling is fixed, and the k-mer backend is
  bit-identical through that change.
- ✅ PacBio is measured at scale — **95 samples, 2.43 M reads, one shared error
  model, 3 timed replicates, control 0.2%** — and is **exactly concordant** at
  cutoff ≥0.45, with 22% fewer alignments than the k-mer screen.
- ✅ The inverted index is verified **exact** at scale: two independent index-off
  arms reproduced their index-on twins cell for cell (`0 / 0 / 0`).
- ✅ Wall clock is measured with control channels on four workloads. The earlier
  claim that "no timing claim is made anywhere on this page" is superseded.
- ✅ Memory is measured at calibrated parameters (−74% resident screen on HiFi).

What promotion would require, in order:

1. **A judgement about the accuracy cost, which is real and not zero *on
   Illumina*.** At the recommended cutoffs: 473-580 churned ASVs of 22,359 on
   pooled 16S, 1.2% count L1; 18 of 3,028 and 0.77% on pooled ITS2. **On PacBio
   HiFi the cost is exactly zero** at ≥0.45, so this requirement is
   platform-specific and is already met on one platform. Churn is confined to the rare tail
   (≤15 reads in every dataset checked) and Bray-Curtis stays under 0.02 per
   sample, but whether that is acceptable depends on what the tables are used
   for. **This is not a question the data answers.**
2. **Replication.** One dataset per configuration. A second diverse pool per
   platform would say whether ~0.64 is a default or a coincidence.
3. **The serial `setup` phase — addressed in code, pending one more workload.**
   Touched-index tracking regressed 29%, and removing the index lost on
   per-sample PacBio (up to 6.0%) and pooled ITS2 (13.9%) — but **won on pooled
   PacBio by 32.8%**, where the indexed arm is 26.4% SLOWER than the k-mer screen
   it replaces
   ([claim 3](#3-the-serial-scatter-is-the-indexs-price-and-should-be-parallelised-away-it-is-the-indexs-bargain)).
   So the index is a bargain on three configurations and a liability on the
   fourth. `minimizers::decide_index` now picks between them from
   `sharing × threads / nraw`, and it is
   [validated on pooled PacBio](#the-index-is-now-chosen-per-workload-not-always-built)
   — the case it exists to catch — landing 0.36% from the right choice and 31.4%
   from the wrong one, for a decision cost of at most ~1s of 290s. What is left
   is **one more workload**: pooled soil 16S, which is outside the bracket the
   threshold was fitted to. Separately, **map efficiency** (90% -> 75%) sits inside the parallel region
   with its own knob in `DADA2RS_PAR_GRAIN`, and index-off holds 89%, so that
   cost is specific to the regime where per-item work collapses to a single array
   read.
4. **A default-selection story.** The right cutoff varies with pass rate
   (0.45-0.80 across workloads), so a fixed default cannot be right everywhere —
   the shipped 0.63 is an Illumina value and is 0.18 too loose on HiFi.
   Either ship the calibration as a required step, or auto-derive the cutoff from
   a cheap pass-rate probe at run start.
5. **Dropping `kord` on this path**, if possible — 2,982 B/raw on HiFi, *larger
   than the sketch itself*, and retained only for the no-indel predicate, which
   now derives both of its counts from it.

Worth splitting out separately: **the gapless decoupling is a fix to existing
code** and stands whether or not this backend ever ships. `kdist` was serving an
undocumented second role as an input to alignment method selection, and any future
screen replacement inherits that trap.

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
