# Template damage makes C>T/G>A variants

**Verdict:** on deep PacBio HiFi amplicon data, almost every low-frequency
single-substitution variant that survives denoising is a **C>T or G>A change at
0.03–0.2% of its parent**. That is the signature of cytosine deamination in the
template. A damaged molecule is copied faithfully, so every read of it reports
the change at high quality, and DADA2 correctly calls it a distinct sequence.
These variants are real molecules but not real alleles. Flagging them after
denoising ([`flag_artifacts.py`](../scripts/flagging.md), class `deamination`)
**doubles precision on a truth-set mock with no loss of recall**:

| ATCC MSA-1003, 20 samples | precision before → after | recall |
|---|---|---|
| Sequel IIe, full qualities | 0.27 → **0.61** | 0.79 (unchanged) |
| Revio, binned qualities | 0.38 → **0.61** | 0.83 (unchanged) |

The same signature is present in an independent library prepared three years
earlier, and accounts for about 17% of the ASVs in a 95-sample real-community
PacBio run. It is what binning quality scores was absorbing in
[Binned quality on PacBio](binned-quality-error-model.md). On Illumina the
deamination rate is similar but hidden under other errors, and the damage that
survives denoising is guanine oxidation (G>T/C>A) instead.

Found while validating the alignment and indel-model work in issues
[#291](https://github.com/HPCBio/dada2-rs/issues/291) and
[#309](https://github.com/HPCBio/dada2-rs/issues/309).

## Data

- **ATCC MSA-1003** (20 strains, staggered abundances), PacBio Kinnex
  full-length 16S (27F–1492R), 20 samples. One set of libraries sequenced
  twice: on a Sequel IIe (full qualities, and a copy binned to PacBio's seven
  levels) and on a Revio (binned). The full-quality Sequel IIe set was
  restricted to the binned copy's read IDs, so the two Sequel IIe arms see
  identical reads.
- **A 2020 Sequel II MSA-1003 run** from a different library prep (9 samples,
  full qualities).
- **PacBio 95**, a 95-sample real-community full-length 16S run.
- **PRJNA674596 V7–V9**, Illumina MiSeq paired reads, 41 samples including
  three Zymo mock samples, with a V7–V9 truth set cut from the Zymo reference.

Denoised with `dada-pooled --band 32 --kmer-size 7` (`--errfun pacbio` for
full qualities, `binned-qual` for binned), consensus chimera removal, scored
with `reference-eval` against the MSA-1003 full-length truth set (52 alleles).

## Binning hid it

On identical Sequel IIe reads:

| | binned | full quality |
|---|---|---|
| chimera-free ASVs | 60 | 154 |
| true alleles recovered | 41 | 41 |
| 1–3 edits from a true allele | 9 | 102 |
| further from any allele | 10 | 11 |

The 93 extra near-miss ASVs in the full-quality run are **all single
transitions, all C>T or G>A**, at a median 0.11% of their parent and spread
evenly along the gene, with no hotspots. The deeper Revio run of the same
libraries shows the same class in binned data too: 43 of its 45 such variants
(96%) are C>T or G>A.

## Why DADA2 keeps them

A deaminated cytosine reads as T, and the opposite strand of the same molecule
gives G>A. The change is in the molecule, not in one read, so the consensus
quality at the site is high: around Q90 on full-quality data. That is where the
error model expects the fewest errors:

| transition rate | binned model, Q40 | full model, Q60 | Q80 | Q90 |
|---|---|---|---|---|
| A>G | 1.8 × 10⁻⁴ | 8.4 × 10⁻⁴ | 3.1 × 10⁻⁴ | **1.3 × 10⁻⁵** |
| C>T | 4.0 × 10⁻⁴ | 7.7 × 10⁻⁴ | 3.6 × 10⁻⁴ | **1.6 × 10⁻⁵** |

A variant at 0.1% of a 100k-read parent is far above 10⁻⁵, so it splits off.
Binning pools every quality at Q40 and above into one rate about 10× higher, and
most of these variants are absorbed instead. The denoiser is behaving
correctly; the variant is a real but artifactual molecule, which no per-read
error model can distinguish from a low-frequency allele.

Which taxa show it is a matter of depth. Detection needs a deep single allele:
*Porphyromonas gingivalis* (78k reads, four identical rRNA copies) has 57 such
variants, the 33–40k-read alleles 8–13, and alleles under ~20k reads almost
none. *E. coli* carries similar depth split across six alleles and shows almost
none. Composition plays no part: *P. gingivalis*, *S. mutans* and
*S. epidermidis* are 51–53% GC.

## It is not one library

The 2020 library is too shallow to resolve these variants as ASVs, so the test
is per site. Within *P. gingivalis* reads at Q≥80, with the Kinnex reads
subsampled to the 2020 depth (12.5k reads), count the positions where several
reads carry the same substitution:

| | positions with the same change in ≥3 reads: C>T/G>A | A>G/T>C | transversions |
|---|---|---|---|
| Kinnex, Sequel IIe (subsample 1) | 266 | 51 | 16 |
| Kinnex, Sequel IIe (subsample 2) | 282 | 49 | 12 |
| 2020 Sequel II, other library | **243** | **31** | **11** |

Both libraries carry a 5–8× excess of C>T/G>A sites. The same test pooled over
all taxa shows only 1.1×, because one-base differences from alleles missing
from the truth set swamp the signal; restricting it to one well-defined allele
is what makes it work.

## Specific to C>T/G>A on PacBio

Low-fraction (≤0.5% of the parent) single-substitution children that track
their parent across samples, counted by type under identical criteria:

| | C>T / G>A | A>G / T>C | transversions (8 types) |
|---|---|---|---|
| PacBio 95, real samples | **365** | 2 | 11 |
| MSA-1003, Sequel IIe full | 94 | 0 | 0 |
| MSA-1003, Revio | 43 | 1 | 1 |
| Illumina V7–V9 (PRJNA674596, 41 samples) | 15 | 9 | 59 |

Polymerase error would not single out C>T/G>A. On the Illumina run the
low-fraction variants are mostly transversions, a different artifact class:
all of them G>T or C>A (next section).

## On Illumina: oxidation

Illumina reads are not consensus reads, so the test is the overlap of a read
pair: a substitution in both R1 and R2 is in the sequenced molecule, one in a
single read is a sequencing error. On the PRJNA674596 V7–V9 Zymo mock (three
samples, 112k pairs; sites carried by more than 1% of pairs excluded, which
removes unlisted alleles and position-specific sequencing errors), against the
*P. gingivalis* allele on full-quality Sequel IIe reads (Q≥80):

| molecular rate, × 10⁻⁵ per base | C>T + G>A | A>G + T>C | G>T + C>A |
|---|---|---|---|
| Illumina V7–V9 (both reads) | 18.0 | 15.7 | **36.7** |
| PacBio Sequel IIe | **23.2** | 5.3 | 2.3 (all transversions) |

Deamination occurs at a similar rate on both platforms, consistent with it
arising in the amplicon PCR rather than in PacBio sequencing; on Illumina it is
hidden under a 3× higher background of other molecular transitions. The excess there
is G>T/C>A, the signature of 8-oxoguanine, a known artifact of Illumina
library preparation. It is what survives denoising: the run's low-fraction
single-substitution ASVs are G>T 41, C>A 21, C>T 22, G>A 9, A>G/T>C 11.

Flagged as `oxidation`, the run gets 46 calls on three abundant parents in the
real samples, at 0.19–0.5% of the parent. The Zymo samples are too shallow per
allele to confirm the class: their five G>T near-misses, at 0.1–0.2%, have fewer
than three siblings and are not called.

## Why the platforms differ

**The denoiser treats every substitution type alike.** What differs is how far
each damage type sits above the learned error model, and how deep its parent
is. The class that becomes ASVs is the one the model underestimates most, at
the qualities where its variants sit.

| molecular rate, × 10⁻⁵ per base | Illumina | PacBio | ASVs on |
|---|---|---|---|
| deamination (C>T/G>A) | 18 | 23 | PacBio only |
| oxidation (G>T/C>A) | 37 | ≈2 | Illumina only |
| other transitions | 16 | 5 | neither |

The backtrace followed every unique one substitution from a truth allele
through a re-run of each mock's pooled `dada` with `--pooled-record`; every
re-run reproduces its original ASV set exactly. A Poisson test of each child's
reads against its expected count, at `OMEGA_A` (10⁻⁴⁰), predicts 167 of the
172 children that split. The expected count is the model's rate for the type
at the site's mean quality times the parent's reads. The split rule is the same
for every type; only the inputs differ.

| | children ≥2 reads | ≥2× expected | split |
|---|---|---|---|
| Sequel IIe full, deamination | 5,688 | 1,438 | 97 |
| Sequel IIe full, other transitions | 1,934 | 349 | 2 |
| Sequel IIe binned (same reads), deamination | 5,688 | 1,385 | 4 |
| Revio binned, deamination | 10,490 | 2,537 | 47 |
| Illumina R1, deamination | 763 | 172 | 5 |
| Illumina R1, oxidation | 971 | **508** | 2 |

- **PacBio, full qualities: the model's shape.** The full-quality model
  expects C>T at about 10⁻⁵ at Q89–92 but 6 × 10⁻⁴ at Q93, the CCS cap where
  most bases sit. 93 of the 97 deamination splits are children whose site
  averages Q89–92; at Q93 2 of 4,718 split. Deamination children are also more
  abundant (1,042 with ≥10 reads, against 70 for other transitions), so at the
  same quality they cross the line and other types do not.
- **PacBio, binned: depth.** Binning judges every base at Q40, about
  3 × 10⁻⁴, so the same children sit only about 3× above expectation. At that
  excess a split needs a parent of roughly 150k reads; Sequel IIe parents
  (median 9k) rarely reach it, the deeper Revio run's do. Every Revio split is
  a child of 50 or more reads. Binned data is now the PacBio norm, so on it
  damage ASVs will grow with depth.
- **Illumina: the model matches deamination, not oxidation.** At Q37–38 the
  learned C>T rate (about 2 × 10⁻⁴) matches the molecular rate, so deamination
  children sit at expectation. G>T is learned at about 1.5 × 10⁻⁴, below its
  molecular rate: half the oxidation children are at least 2× above
  expectation, at a median 5.7×, against a fifth of deamination children at
  3.5×. Splitting at that excess needs parents of roughly 125–200k reads; the
  mock alleles (about 9k reads each) are too shallow, so its oxidation
  children are absorbed. The 46 mouse-sample calls have parents of 14–34k
  merged reads, at the same per-site fractions as the mock's children.
- **Paired-end denoising: no detectable effect.** Children in the read overlap
  split as often as those outside it (4 of 1,510 against 7 of 2,113 in R1), but
  too few split to measure a difference.

### The Q89–92 dip on full-quality PacBio

The dip is in the counts, not the fit. learn-errors tabulates each difference
at the unique's mean quality at that site, as DADA2 does. On full-quality
Sequel IIe reads:

| | Q90 | Q91 | Q92 | Q93 (cap) |
|---|---|---|---|---|
| share of bases | 2% | 4% | 18% | 70% |
| C>T, counted rate | 2.9 × 10⁻⁵ | 3.8 × 10⁻⁵ | 4.5 × 10⁻⁶ | 5.8 × 10⁻⁴ |

A molecular difference (damage, a PCR error, a member lumped into the wrong
cluster) is present in every pass of the molecule, so CCS gives it full
quality: Q93 collects them, and its counted rate is the molecular rate, not
the sequencing rate. DADA2's `PacBioErrfun` estimates Q93 directly from these
counts and fits LOESS to Q0–92 only, weighted by bases observed, so the 22M
Q92 bases also pull Q85–91 2–4× below their own counts.

Averaging then moves abundant variants into the dip. Among deamination
children, the share whose mean site quality falls below Q93 rises from 4% for
single-read children to 42% at 10–50 reads and 72% above 50, and every split
is a child of at least 10 reads. An abundant molecular variant is therefore
judged against a sequencing-error rate about 100× below the molecular rate.
Quality on HiFi reads measures agreement between passes, not whether the
molecule was copied faithfully.

Binned models have no dip: rates fall monotonically to Q40, and qualities
between bins (from averaging) get higher rates, not lower (C>T 6 × 10⁻⁴ at
Q39 against 4 × 10⁻⁴ at Q40), so averaging makes a split harder.

Why learn-errors underestimates G>T on Illumina is not established.

**What these data cannot answer.** The PacBio and Illumina mocks differ in more
than platform: mock (ATCC vs Zymo), lab, amplicon length (1450 vs 400 nt),
polymerase, cycle number, and Illumina's index PCR. The backtrace explains why
a variant was absorbed or split; whether the different oxidation rates come
from the platform or the preparation needs one library sequenced on both.

## Flagging

`flag_artifacts.py` labels a candidate `deamination` when it is one C>T or G>A
substitution from a more abundant ASV, at or below `--damage-max-frac` (0.005)
of it, its parent has at least `--damage-min-siblings` (3) such candidates, and
it does not vary independently of the parent across samples; `oxidation` is the
same for G>T or C>A. Damage hits many sites of the same template, so the
siblings are the evidence, as recurrence across parents is for homopolymer
slippage. `--action remove` or `merge` with `--include-damage` removes them, or
merges them into their parent.

On the tables above, every `deamination` call on the mocks was a near-miss of a
true allele, never a true allele, and the Illumina run gets none. PacBio 95
gets 365 (about 23k reads), in every alignment configuration tested.

Both classes are **highly experimental**. The evidence on the mocks is strong,
but damage rates are expected to vary between preparations, facilities and
labs, and applying the classes to a data set without a mock needs more
assessment.

## The error model also counts them

`learn-errors` counts every cluster member's differences from its centre as
errors, so damage adds C>T/G>A events to the model. Against the truth set the
learned models are calibrated overall at Q40 but wrong by type: transversions
are 4–30× too high, from members far from their centre
([#319](https://github.com/HPCBio/dada2-rs/issues/319)), and two transitions
are about 2× too low. Damage is one more reason a per-read error model cannot
describe what is in an amplicon library.

## Consequences

- **Flag damage after denoising, like chimeras.** Lowering the denoiser's
  sensitivity would also lose real low-frequency alleles; the type, fraction
  and sibling pattern separate damage at no cost to recall on the mocks.
- **Expect it more on full-quality and deep data.** Binning hides most of it;
  added depth reveals more.
- **Expect it to vary by lab.** Damage depends on how samples and libraries
  are handled, so a mock's rates are evidence for that preparation, not a
  constant.
- **Not tested:** whether it varies with sample storage, extraction or PCR
  conditions; other polymerases and library kits; whether the 2020 library
  shows it at the ASV level when sequenced deeper; and oxidation on a deep
  Illumina mock.

## Reproducing

- Truth scoring: `dada2-rs reference-eval` against
  `msa1003_FL_27F-1492R.fasta`.
- Flagging: `scripts/flag_artifacts.py nochim.json --out flags.tsv
  [--action remove --include-damage --out-seqtab filtered.json]`.
- Data: `mocks/PacBio-SQII-MSA-1003-{binned,full}`,
  `mocks/PacBio-Kinnex-MSA-1003` (Revio) and `mocks/PacBio-MSA-1003` (2020) in
  the project's data root.
