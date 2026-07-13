# Diagnostics

This page documents **experimental diagnostic tooling** built into `dada2-rs`, along with related helper scripts, for interrogating the internals of the DADA2 algorithm. 

---

## `--failed-uniques`: which sequences failed to denoise

A recurring question (e.g. [benjjneb/dada2#1899](https://github.com/benjjneb/dada2/issues/1899)) is *which* unique sequences fail to denoise — particularly for high-diversity samples (soil, seawater) that shed substantial reads through denoising. DADA2's answer lives in the per-unique `map`: a unique whose final abundance p-value falls below `OMEGA_C`, and that is not corrected to any cluster center, gets a `null` map entry (R's `$map == NA`). Those nulls are the uniques that failed to denoise, and they can be traced back to the reads they cost.

`dada2-rs` exposes this directly. The `--failed-uniques <FILE>` flag — available on `dada`, `dada-pseudo`, and `dada-pooled` — writes a TSV of the dropped uniques and their abundances, so you don't have to hand-join the `map` to the derep uniques.

### Output

A tidy long-format TSV with a header, one row per failed unique per sample it appears in:

```
sequence	sample	reads
TACGAAGG…AACA	soil_A	7
TACGGAGG…AATC	soil_A	3
TACGGAGG…AAAC	soil_B	2
```

| Column | Meaning |
|---|---|
| `sequence` | The unique that failed to denoise (`map == null`). |
| `sample` | Sample the unique (and its reads) belongs to. |
| `reads` | Reads this unique contributed in that sample (i.e. reads lost to the failure). |

Rows are sorted by `sample`, then descending `reads`, then `sequence`, so output is deterministic regardless of how many samples were denoised concurrently. Sum the `reads` column (optionally per sample) to get the total reads lost to failed denoising.

### Semantics follow the subcommand

The notion of "failed" matches how each subcommand denoises:

| Subcommand | "Failed" decided | `reads` is |
|---|---|---|
| `dada`, `dada-pseudo` | **per sample** (`map == null`; for pseudo, in the round-2 pass) | the unique's in-sample abundance |
| `dada-pooled` | **globally** — pooled denoising runs once on the merged unique table, so failure is a property of the merged index (`result.map == null`) | the failed merged unique's read count *in that sample* (one row per sample it appears in) |

!!! note "Per-sample `null` is unambiguous"

    In the pooled per-sample output JSON, a `map` entry could in principle be
    `null` either because the unique failed denoising or because its cluster has
    zero reads in that sample. The latter never actually happens for a unique
    that is present in the sample (a present unique always carries ≥1 read into
    its cluster), so a per-sample `map == null` always means a genuine denoising
    failure. The `--failed-uniques` TSV reports only real failures.

### Invocation

```bash
# single sample
dada2-rs dada sampleA.fastq.gz --error-model err.json \
    -o sampleA.dada.json --failed-uniques sampleA.failed.tsv

# multi-sample (per-sample failures, one combined TSV with a sample column)
dada2-rs dada *.fastq.gz --error-model err.json \
    --output-dir dada_out/ --failed-uniques failed.tsv

# pooled (global failures, expanded per sample)
dada2-rs dada-pooled *.fastq.gz --error-model err.json \
    -o pooled_out/ --failed-uniques failed.tsv
```

To trace a failed unique's divergence from the nearest surviving ASV (is it a distant artifact, or a real low-abundance variant the abundance test shed?), feed the same `dada` run to [`kdist-calibrate --from-dada`](#3-post-inference-mode-from-dada), whose `failed` class is exactly this population.

---

## `kdist-calibrate`: k-mer screen, `KDIST_CUTOFF`, `BAND_SIZE`

Currently the k-mer **screen** (`KDIST_CUTOFF`) and the alignment **band size** 
(`BAND_SIZE`) are two constants that ship with set defaults or recommendations, but without much documentation on how these were derived (which data were used, the tools for calibration, etc). Some of the tooling, e.g., ESPRIT, also are no longer available online. 

!!! note "Diagnostic subcommand"

    `kdist-calibrate` is a **diagnostic** subcommand: it exists to *measure* and
    *characterize* algorithm behaviour, not to denoise, and is not part of the
    denoising pipeline itself. Treat its results as exploratory — see the
    [caveats](#caveats) before drawing conclusions, especially the small sample
    sizes behind the preliminary numbers below.


### What it measures

DADA2 avoids most pairwise alignments with an inexpensive **k-mer distance screen**:
pairs whose k-mer distance exceeds `KDIST_CUTOFF` (default **0.42**, nominally
"~10% nucleotide divergence", calibrated on Illumina 16S) are assumed too
different to be linked by amplicon error and skipped. Surviving pairs are
aligned within a diagonal **band** (`BAND_SIZE`, default **16**, with PacBio HiFi data recommendations being **32**).

Both constants raise questions this tool answers empirically:

- **What divergence does `kdist = 0.42` actually correspond to** on *your* data,
  platform, `k`, and pooling regime? (as noted, the original ESPRIT reference
  implementation that defined the k-mer distance is no longer available.)
- **How much headroom does 0.42 have** above the real error-copy distances —
  i.e. how far could it be tightened without dropping true error copies?
- **Is `BAND_SIZE = 16` the right size** — does it cover real error-copy
  alignments, and is it over-provisioned for short reads?

It answers these by re-deriving, on real sequences, the relationship between the
k-mer distance and the **true unbanded NW alignment divergence**.

---

### Invocation

The input is one or more **derep JSON** files (`.json` / `.json.gz`, as produced
by the `derep` subcommand):

```bash
dada2-rs kdist-calibrate sampleA.derep.json [sampleB.derep.json ...] \
    --k 5 --threads 8 --verbose -o kdist.csv
```

| Flag | Default | Meaning |
|---|---|---|
| `--k` | `5` | k-mer size (DADA2 fixed k-mer size is 5; full-length PacBio wants 7). |
| `--cutoff` | `0.42` | Screen cutoff used for the `screened_in` flag and summaries. |
| `--band` | `-1` | Alignment band radius; **negative = unbanded** (the correct default — a band would truncate the divergence of distant pairs). |
| `--max-pairs` | `200000` | Max pairs computed **per population** (random-subsample above this to bound the O(n²) cost). |
| `--max-uniques` | `0` (all) | Randomly subsample each sample to at most this many uniques before pairing. |
| `--per-sample` | off | Compute pairs **within** each sample (independent regime) instead of pooling all uniques into one set (full-pool regime). |
| `--nearest-parent` | off | Abundance-aware mode (see below). |
| `--from-dada` | off | Post-inference mode (see below); positional inputs are `dada` output JSONs. |
| `--from-dada-pooled` | off | Pooled post-inference mode (see below); positional inputs are the pooled record(s) written by `dada-pooled --pooled-record <path>`. Self-contained — no `--derep-dir`. |
| `--derep-dir` | — | With `--from-dada`: directory of the derep JSONs that fed `dada` (matched by sample name). |
| `--threads` | `1` | Threads for the parallel Needleman–Wunsch alignments. |
| `--seed` | fixed | RNG seed for subsampling (reproducible). |
| `-o`, `--output` | stdout | Write the CSV here. |
| `--verbose` | off | Print per-population progress and summaries to stderr. |

!!! note "Pooling mode = which sequences you feed it"

    The pair *population* the screen sees depends on the denoising mode.
    `--per-sample` mirrors per-sample/independent denoising; the default pools
    all inputs (full-pool). Pseudo-pooling's screen population is per-sample
    (priors change the partition, not which pairs are screened), so model it
    with `--per-sample`.

---

### Modes and outputs

#### 1. All-pairs mode (default)

For sampled unique-sequence pairs, emits the k-mer distance alongside the true
alignment divergence. CSV columns:

| Column | Meaning |
|---|---|
| `sample` | Population label (`pool`, or the sample name under `--per-sample`). |
| `kdist` | k-mer screen distance (see formula below). |
| `edits` | Substitution + indel columns in the aligned core. |
| `core_len` | Aligned-core length (terminal overhang trimmed). |
| `pct_div` | `100 · edits / core_len` — true percent divergence. |
| `band_req` | Minimum diagonal band that reproduces this alignment (see [band size](#band-size)). |
| `screened_in` | `1` if `kdist < cutoff` (DADA2 would align this pair). |
| `ab_i`, `ab_j` | The two sequence abundances. |

The k-mer distance itself (the ESPRIT metric DADA2 ports) is:

$$
\text{kdist} = 1 - \frac{\sum_i \min(c^{1}_{i},\, c^{2}_{i})}{\min(L_1, L_2) - k + 1}
$$

- $c^{1}_{i}, c^{2}_{i}$ — the **count** (within-sequence multiplicity) of k-mer
  $i$ in sequence 1 and sequence 2. The metric is **multiplicity-aware**, not
  presence/absence: the numerator is a multiset intersection, so a k-mer
  occurring 3× in one sequence and 2× in the other contributes
  $\min(3, 2) = 2$ shared occurrences.
- **numerator** — total shared k-mer occurrences between the two sequences.
- **denominator** — total k-mer positions in the shorter sequence,
  $\min(L_1, L_2) - k + 1$.

!!! note "Two kinds of "abundance""

    The counts $c^{1}_{i}/c^{2}_{i}$ are **within-sequence k-mer
    multiplicities** — a property of the two sequences' composition only. They
    are unrelated to the sequences' read abundances (the `ab_i`/`ab_j` columns),
    which do not enter the screen distance. (The `kmer8` storage saturates a
    per-k-mer count at 255, which no realistic amplicon reaches.)

#### 2. Abundance-aware mode (`--nearest-parent`)

For each unique, links it to its nearest **more-abundant** neighbour — its
candidate error-copy "parent", mirroring DADA2's greedy center-based comparison
— and aligns that one pair. The distribution of parent-link distances is the
empirical **error-copy distance ceiling**; `cutoff − ceiling` is the screen's
**headroom**. CSV columns:

| Column | Meaning |
|---|---|
| `sample` | Population label. |
| `ab` | Abundance of the (child) unique. |
| `parent_ab` | Abundance of its nearest more-abundant neighbour. |
| `ab_ratio` | `parent_ab / ab` (larger ⇒ more plausibly an error copy). |
| `kdist`, `edits`, `core_len`, `pct_div`, `band_req`, `screened_in` | As above, for the child↔parent link. |

This mode replaces the O(n²) *alignments* of all-pairs mode with only O(n)
alignments — one per unique — but it still does O(n²) *cheap* k-mer comparisons
to find each parent, and (unlike all-pairs mode) it **ignores `--max-pairs`**:
every unique is scanned against its entire more-abundant prefix, uncapped.

!!! warning "Run `--nearest-parent` with `--per-sample` (or `--from-dada-pooled` for pooled)"

    Because the scan is uncapped, feeding it a **pooled raw-derep union** is a
    trap: `n` is the sum over all samples (e.g. ~750k uniques across 384 MiSeq
    samples), so the O(n²) k-mer scan is enormous, and the "nearest more-abundant
    parent" of a unique may live in a **different sample** — a cross-sample link
    that is not a real error copy. Always pair `--nearest-parent` with
    `--per-sample`, where `n` is per-sample small and every parent link is a
    genuine within-sample error-copy candidate. For a *pooled* regime, don't use
    the raw-derep union at all — use [`--from-dada-pooled`](#4-pooled-post-inference-mode-from-dada-pooled),
    which scores the merged pool once with pooled abundances. `--max-pairs` does
    **nothing** here (passing it alongside `--nearest-parent` is rejected); to
    bound the scan on a large population, cap it with `--max-uniques` instead.

```bash
# Abundance-aware, per sample (the correct regime for --nearest-parent)
dada2-rs kdist-calibrate derep/*.derep.json.gz \
    --k 5 --nearest-parent --per-sample --threads 24 -o kdist.abundance.csv
```

#### 3. Post-inference mode (`--from-dada`)

The two modes above operate on **derep** inputs — the screen's view, *before*
denoising. This mode operates on `dada` **output**, so every comparison is
labelled by *what inference actually decided*. Each input unique gets one of
three fates, and is aligned against the relevant cluster center:

| `class` | Compared against | What it is |
|---|---|---|
| `member` | its own cluster center | a real error copy denoising **corrected** (the within-cluster cloud). |
| `failed` | its **nearest** center | a unique the abundance test **shed** (`map == null`) but did not assign. |
| `center_pair` | another surviving center | two ASVs that both **survived** — the inter-ASV resolution floor. |

Because labels come from dada's actual abundance-p-value partition (not the
nearest-more-abundant *proxy* of mode 2), this is the ground-truth version of
the headroom question — and the `failed` class is the population mode 2 can't
see at all.

Invocation pairs each `dada` output JSON with its derep input (located in
`--derep-dir` by sample name, so indices line up with dada's `map`):

```bash
dada2-rs kdist-calibrate --from-dada \
    dada_out/sampleA.json [dada_out/sampleB.json ...] \
    --derep-dir derep/ \
    --k 5 --threads 8 --verbose -o post.csv
```

CSV columns:

| Column | Meaning |
|---|---|
| `sample` | Sample name (from the dada output). |
| `class` | `member` / `failed` / `center_pair` (above). |
| `cluster` | Cluster id of the partner center. |
| `ab` | Abundance of the query unique (or center *a* for `center_pair`). |
| `center_ab` | Abundance of the partner center. |
| `ab_ratio` | `center_ab / ab`. |
| `birth_type` | How the partner ASV was born: `Abundance`, `Prior` (pseudo-pool prior), `Initial`, or `Singleton`. |
| `birth_pval` | The abundance p-value at that birth — small ⇒ a confident split; near `OMEGA_A` ⇒ a **borderline** ASV. |
| `kdist`, `edits`, `core_len`, `pct_div`, `band_req`, `screened_in` | As above, for the query↔partner alignment. |

**Tracing priors / pseudo-pooling.** Run `dada-pseudo`, then filter the table on
`birth_type == Prior`: those ASVs exist *only* because a prior from another
sample rescued them past `OMEGA_P`. Their `center_pair` rows show how close each
sits to the nearest abundance-born survivor — small divergence there means the
prior recovered a real low-abundance variant; large divergence is worth a closer
look. `birth_pval` lets you sort every ASV by how borderline its split was,
independent of the prior question.

#### 4. Pooled post-inference mode (`--from-dada-pooled`)

`--from-dada` runs on the **per-sample** `dada` outputs. For a `dada-pooled` run
that is subtly wrong for a *pool-level* question: pooled denoising decides every
unique's fate **once**, on the merged unique table, and then writes each
per-sample JSON as a *projection* of that single global partition. Feeding those
per-sample files back to `--from-dada` re-fragments the pool — a sequence shared
by *N* samples is loaded *N* times (once per sample's derep, at that sample's
local read count), so it contributes *N* CSV rows and has its `failed` singleton
/ multi-read split re-derived from each local count instead of its pooled
abundance. In a two-sample run, summing the per-sample populations gives
1805 unique-rows; the merged table has only 1730 — the 75 difference is exactly
the sequences present in both samples, double-counted.

To assess the pool as it was actually denoised, ask `dada-pooled` for a
self-contained record with `--pooled-record <path>`. It carries the merged
unique table — pooled abundances, the global map, and the global ASV list — in
one file. It is **off by default**, and you give it an explicit path (gzip
follows a `.gz` extension); keep that path **outside** `--output-dir` so it does
not join the per-sample `*.json.gz` glob that downstream sequence-table /
aggregation steps rely on:

```bash
dada2-rs dada-pooled sampleA.derep.json sampleB.derep.json ... \
    --error-model err.json --output-dir dada_out/ --gzip \
    --pooled-record kdist/pooled.json.gz
```

```
dada_out/                 kdist/
├── sampleA.json.gz       └── pooled.json.gz   # the merged unique table,
└── sampleB.json.gz                            # kept clear of dada output
```

`--from-dada-pooled` reads that record directly — no `--derep-dir`, because the
merged uniques are carried inline — and screens the whole run as **one**
population (`sample = __pooled__`) with pooled abundances and global fate:

```bash
dada2-rs kdist-calibrate --from-dada-pooled \
    kdist/pooled.json.gz \
    --k 7 --threads 8 --verbose -o post.pooled.csv
```

The CSV columns are identical to `--from-dada` (§3); only the population differs
(one deduplicated merged table vs. one block per sample). `member` / `failed` /
`center_pair` mean the same thing, and the `failed` singleton split now reflects
each unique's **pooled** read count. Use `--from-dada` when you want the
per-sample view (e.g. which samples a failure appears in); use
`--from-dada-pooled` when you want the pool's true screen behaviour without
double-counting shared sequences.

> **Note.** This is a `dada-pooled` concept. `dada-pseudo`'s screen population is
> genuinely per-sample (priors reshape each sample's partition), so it has no
> `--pooled-record` — use `--from-dada` there and filter on `birth_type == Prior`
> (§3) to trace prior-born ASVs.

#### `--verbose` summaries (stderr)

`--verbose` adds per-population summary lines that are usually what you want
before touching the CSV:

**All-pairs mode:**

```
[kdist] pool : 449 uniques, 100576 pairs -> 100576 computed (k=5, band=-1, 4 threads)
[kdist] 100576 pairs: screened-in (kdist<0.42) 81470 (81.0%); of those 44555 are >5% divergent (leakage)
[kdist] all pairs band-fit (100576, max_req 3): ≤2:99.9% ≤4:100.0% ≤8:100.0% ≤16:100.0% ...
[kdist] screened-in band-fit (81470, max_req 2): ≤2:100.0% ...
```

- **screened-in** — fraction of pairs the screen would align.
- **leakage** — of those, the fraction too divergent (`> --leak-pct`, default
  5%) to be an error copy = wasted alignments.
- **band-fit** — for candidate bands `[2,4,8,16,32,64,128]`, the fraction of
  alignments whose true path fits (i.e. a banded aligner at that size would
  compute correctly), plus the maximum band required.

**Abundance-aware mode:**

```
[kdist] pool : 448 children | nearest-parent kdist median 0.021 p90 0.064 | 445 (99.3%) within cutoff 0.42 | clear-error-copy ceiling 0.127 -> headroom 0.293
[kdist] pool : clear-error-copy band-fit (435, max_req 1): ≤2:100.0% ...
```

- **ceiling / headroom** — the max k-mer distance among clear error-copy links
  (≤3% divergent), and how far that sits below the cutoff.
- **clear-error-copy band-fit** — the band-fit curve restricted to real error
  copies (the safety-relevant question for `BAND_SIZE`).

**Post-inference mode:**

```
[kdist] sam1F : 896 uniques (9 centers, 828 members, 59 failed), 9 ASVs, 923 jobs (k=5, band=-1, 4 threads)
[kdist] sam1F : 59 failed | singletons 59 (14 within cutoff) | multi-read 0 (0 within cutoff) — failed singletons are the --detect-singletons tradeoff, not distance
[kdist] sam1F : 2/9 ASVs born from priors (pseudo); filter the table on class=center_pair,birth_type=Prior to see their nearest survivor
```

- the fate breakdown (centers / members / failed) per sample, plus the job
  count actually aligned.
- the failed class split by abundance. A unique that fails the abundance test
  for being a **singleton** (the default `≥2 reads` rule, toggled by
  [`--detect-singletons`](parameters.md)) is a different thing from one that is
  genuinely distant from every center — so failed singletons are reported
  separately, with how many sit *within* the screen cutoff (near a center, i.e.
  plausible error copies / real low-abundance variants that just lacked a second
  read) vs beyond it (the distant tail).
- a prior line appears only when some ASV was born from a pseudo-pool prior.

**Pooled post-inference mode (`--from-dada-pooled`):** the same lines, but one
`__pooled__` block over the merged unique table — the fate counts and the
singleton split are pool-wide (pooled abundances), not per-sample:

```
[kdist] __pooled__ : 1730 uniques (11 centers, 1601 members, 118 failed), 11 ASVs, 1774 jobs (k=7, band=-1, 8 threads)
[kdist] __pooled__ : 118 failed | singletons 118 (30 within cutoff) | multi-read 0 (0 within cutoff) — failed singletons are the --detect-singletons tradeoff, not distance
```

!!! note "`failed` ≠ distant noise"

    Most `failed` uniques are typically **singletons**: under the default a
    singleton cannot seed a new ASV regardless of distance, so it lands in
    `failed` for the read-count tradeoff, not because the screen judged it far
    from everything. Use the singleton split (and re-running `dada` with
    `--detect-singletons`) to tell the two apart before reading anything into the
    failed-class divergence.

---

### Processing the CSV

The CSV is deliberately raw, one row per pair, so any tabular tool works. A
common first cut — the k-mer-distance ↔ divergence calibration curve — bins by
`kdist` and reports the median divergence:

```python
import pandas as pd
df = pd.read_csv("kdist.csv")
bins = [0, .1, .2, .3, .4, .42, .44, .5, .6, .8, 1.01]
df["bin"] = pd.cut(df.kdist, bins)
print(df.groupby("bin", observed=True).pct_div.median())

# What divergence does the cutoff correspond to here?
near = df[(df.kdist >= 0.41) & (df.kdist < 0.43)]
print("kdist≈0.42 ->", near.pct_div.median(), "% divergence")

# Leakage: screened-in but clearly not an error copy
si = df[df.screened_in == 1]
print("leakage:", (si.pct_div > 5).mean())
```

For `--nearest-parent` output, the headroom and band questions fall out
directly:

```python
np_df = pd.read_csv("kdist_np.csv")
ec = np_df[np_df.pct_div <= 3]            # clear error copies
print("error-copy kdist ceiling:", ec.kdist.max())
print("band needed for error copies:", ec.band_req.max())
```

---

### Preliminary outcomes

!!! danger "Preliminary — single-sample illustrations"

    The single-sample numbers in this section come from **one Illumina sample**
    (449 uniques) and **one PacBio sample** (259 uniques). They illustrate the
    *kind* of result the tool produces; they are **not** a basis for retuning any
    default. The [384-sample MiSeq run](#at-scale-a-384-sample-miseq-run) below is
    the first at-scale check — it *confirms* the single-sample trends rather than
    overturning them, but is still one platform / amplicon. See [caveats](#caveats).

#### The screen cutoff (`KDIST_CUTOFF = 0.42`)

On an Illumina V4 sample (sam1F, 240 bp, k=5):

| k-mer distance bin | median divergence |
|---|---|
| 0.00–0.10 | 1.2% |
| 0.30–0.40 | 10.0% |
| **0.40–0.42** | **14.2%** |

`kdist = 0.42` corresponds to **~14.6%** divergence in this case — *not* the 
nominal 10% (which sits at kdist ≈ 0.29). The cutoff is safe (every ≤3%-divergent pair
is screened in) but **measurably looser** than its documented calibration, and
~55% of screened-in pairs are too divergent to be error copies (pure leakage).

#### Screen saturation on long reads

The k-mer distance is `1 − shared/(L−k+1)`. When sequence length `L` ≫ `4ᵏ`,
every sequence contains nearly the whole k-mer vocabulary, so even very
divergent pairs share most k-mers and the distance **saturates** below the
cutoff:

| PacBio samPB (1464 bp) | max kdist reached | screened-in | kdist=0.42 means |
|---|---|---|---|
| **k = 5** | 0.339 (never hits 0.42) | **100%** (blind) | unreachable |
| **k = 7** | 0.688 | 43% (prunes 57%) | **~11% divergence** |

At the stock k=5, the screen on full-length 16S prunes *nothing* in this case — every pair is
aligned (correctness is unaffected; only compute is wasted). At k=7
(`4⁷ = 16384 ≫ 1460`) the screen de-saturates **and** the 0.42 cutoff lands back
near its intended ~10%.

!!! note "k is a `dada2-rs` setting"

    Tunable `k` is a `dada2-rs` feature; upstream R DADA2 hard-codes `k = 5`.
    These measurements are the mechanistic basis for our PacBio `k = 7`
    recommendation.

#### Screen headroom (abundance-aware)

| sample (k) | error-copy kdist ceiling | headroom below 0.42 |
|---|---|---|
| Illumina sam1F (k=5) | 0.127 | **0.293** |
| PacBio samPB (k=7) | 0.111 | 0.309 |

Real error copies sit far below the cutoff (median parent-link kdist ≈ 0.02 on
Illumina). On these samples the cutoff could be roughly **3× tighter** and still
capture every clear error copy — though that figure is a per-sample *lower
bound*, not a safe global value.

#### Band size

`band_req` is the minimum band that reproduces a pair's true alignment:

| population | max band needed | covered by band 8 | covered by band 16 |
|---|---|---|---|
| Illumina, all pairs | **3** | 100% | 100% |
| Illumina, error copies | **1** | 100% | 100% |
| PacBio (k=7), error copies | **10** | 98.8% | 100% |

The default `BAND_SIZE = 16` is **over-provisioned for Illumina** (error copies
need ≤1) and **appropriately sized for PacBio** (CCS homopolymer indels push
band_req to 10; band 8 would miss ~1% of real error copies). Since DP cost is
O(L·band), a smaller short-read band is a direct speed-up. Like the cutoff, 16
looks like a single worst-case constant applied uniformly.

#### At scale: a 384-sample MiSeq run

The first multi-sample check — a 384-sample MiSeq 16S run (2×250, `k = 5`),
scored independently for R1 and R2. The all-vs-all geometry used `--max-pairs
400000` over the ~753k-unique pool; abundance mode used `--per-sample`; the
pooled partition came from `dada-pooled --from-dada-pooled`. Every single-sample
trend above reappears — and holds across **both** reads.

**Cutoff (`kdist = 0.42`) is looser than its "~10%" label.**

| read | uniques | kdist 0.42 → divergence | 10% divergence → kdist |
|---|---|---|---|
| **R1** | 753,645 | **14.9%** | 0.309 |
| **R2** | 752,747 | **13.1%** | 0.327 |

Same direction and magnitude as the single Illumina sample (~14.6%). The nominal
10% actually sits at kdist ≈ 0.31–0.33. R2's distance saturates slightly faster
per unit divergence (noisier tail), but the conclusion is identical for both reads.

**The cutoff is very safe — real error copies sit far below it** (abundance mode,
`--per-sample`):

| read | nearest-parent kdist (median / p90) | within cutoff 0.42 | clear error copies missed by screen |
|---|---|---|---|
| **R1** | 0.021 / 0.123 | 99.2% | 1 / 538,762 (**0.000%**) |
| **R2** | 0.032 / 0.147 | 99.3% | 4 / 408,977 (**0.001%**) |

The handful of "missed" clear error copies are `core_len ≤ 2` degenerate
alignments (the `core == 0 → pct_div = 0` fallback), not real full-length error
copies. In practice the screen misses **nothing**, confirming ~0.28 of headroom —
0.42 could be tightened substantially before touching a genuine error copy.

**Band 16 is over-provisioned for MiSeq.** Across 400k random pairs, `band_req`
median = 1 and p99 = 3–4; **band ≤ 4 covers 100%** of screened-in pairs and all
clear error copies (the max_req 158–240 outliers are distant junk pairs that
never survive inference). Consistent with the single Illumina sample (error
copies needed ≤ 1).

**The screen causes no denoising failures.** In the `--from-dada-pooled`
partition, every *multi-read* failed unique is *within* cutoff (1000/1000 R1,
835/835 R2) — those failures are the pooled abundance p-value, not distance.
Failed singletons are the `--detect-singletons` tradeoff. The screen is not
responsible for any shed sequence.

This is a strong at-scale confirmation on **short-read MiSeq 16S**. Long reads
behave differently — the PacBio run below is the complementary check.

#### At scale: a 74-sample PacBio HiFi run (`k = 7`)

A 74-sample PacBio HiFi full-length 16S run (~1.5 kb, ~857k uniques, `k = 7`).
This is the long-read counterpart to the MiSeq set, and it validates the two
things that are *specific* to long reads: that `k = 7` de-saturates the screen,
and that the band must be larger.

**`k = 7` de-saturates the screen and lands the cutoff near ~10%.** Where the
single-sample PacBio at `k = 5` never exceeded kdist 0.339 (screen blind, prunes
nothing), at `k = 7` the pool reaches **kdist 0.955** and the screen prunes
**85.6%** of pairs (14.4% screened-in). The cutoff `0.42` corresponds to **12.5%**
divergence (10% sits at kdist 0.374) — the closest to its nominal label of any
dataset here, and only 60% of screened-in pairs are >5%-divergent (vs 89% on
MiSeq): at `k = 7` the screen is both *functional* and *less leaky* on long reads.

**The cutoff is safe with even more headroom than MiSeq.** Real error copies sit
extremely close to their parent — nearest-parent kdist **median 0.005**, 99.6%
within cutoff, and **0.000%** of clear error copies (with a non-degenerate core)
are screened out. Per-sample clear-error-copy ceilings run ~0.10–0.16, so ~0.28
of headroom, matching MiSeq.

**Band 16 is appropriately sized here — not over-provisioned.** This is the sharp
contrast with MiSeq — but reading it correctly matters. Looking at the *real*
inference alignments (member→center pairs from `--from-dada-pooled`, 528k of
them), genuine error copies (≤3% divergence) need `band_req` median 1, p99 **7**,
and **band ≤ 8 covers 99.5%** of them. The long tail (`band_req` into the
hundreds) is **not** homopolymer indels — HiFi consensus reads are not
homopolymer-error-prone (unlike CLR / 454 / Nanopore). It is dominated by:

- **absorbed noise** (>3% divergence, ~0.3% of members): off-target / chimeric
  **singletons** parked in the nearest cluster because they are not abundant
  enough to split off — median ~23% divergence, needing bands of 500–1600;
- **truncation / offset artifacts**: partial-length or off-target reads that are
  *identical over their overlap* but shifted by hundreds of bp (`edits ≈ 0`,
  short `core_len`), so the ends-free alignment carries a large terminal gap.

Crucially, these tail cases need `band_req` of 179–1676 — which the default `16`
(or `32` for HiFi) **already does not reach**, so they are aligned sub-optimally
under any practical band regardless. They therefore do **not** constrain
`BAND_SIZE`, and they are not error copies of a real center. (The `--nearest-parent`
proxy inflates this tail further because it pairs singleton↔singleton, a
comparison that never occurs in real center-based correction.) So for the
*alignment-reach* question, genuine HiFi error copies need ≈ the same small band
as Illumina (≤ 8).

!!! warning "`band_req` is necessary but NOT sufficient — the end-to-end A/B overrules it"

    An earlier version of this page concluded from the above that `BAND_SIZE`
    could be lowered on HiFi too. **An end-to-end band A/B refutes that.**
    Sweeping `--band` through both error learning and `dada-pooled` (dev
    `run_band_sweep.sh`) and diffing the ASV set:

    - **MiSeq (384-sample, k=5):** band 16→8 preserves the ASV set on both R1/R2
      (≤ 0.0013% read drift); band 4 starts to break (loses a borderline ASV on
      R2). Band 8 is a safe floor.
    - **PacBio HiFi (74-sample, k=7):** band 32→16 **already changes the ASV
      catalogue** — −14/+7 ASVs, 0.02% read drift, and it worsens at band 8
      (−9/+11, 0.06%). The lost ASVs are `Abundance`-born, often strongly
      significant, and merge into a **Hamming-1 neighbour** already in the set.

    Why the divergence, when HiFi error copies need only band ≤ 8? Because the
    band perturbs two things: alignment reach (mechanism 1, weak here) **and the
    learned error model** (mechanism 2). On 1.5 kb reads the error-model
    perturbation is ~10–20× larger than on MiSeq (`err_out` max Δ ~2e-3 vs
    ~1e-4), and that shift flips the abundance test for close (Hamming-1) ASV
    *pairs*, merging variants that band 32 resolves. So band 32 earns its keep on
    HiFi not through alignment reach but by keeping the error model stable enough
    to separate near-identical variants. **The platform-aware default (16
    Illumina / 32 HiFi) is correct**; `band_req` alone understates the impact.

**The screen causes no denoising failures.** In the `--from-dada-pooled`
partition, all 16,235 failed uniques are *within* cutoff (15,955 singletons +
280 multi-read) — abundance-driven, not the screen.

!!! note "Two platforms — cutoff agrees, band diverges"

    MiSeq and PacBio agree on the screen: `KDIST_CUTOFF = 0.42` is safe with large
    headroom on both, and `k = 7` is what makes it meaningful on full-length 16S.
    They **disagree on the band**, which is the point of the platform-aware
    default. The end-to-end A/B (see the warning above) shows `BAND_SIZE` can drop
    to 8 on MiSeq without touching the ASV set, but on HiFi even 32→16 changes the
    catalogue — so `16` (Illumina) / `32` (HiFi) is correct, and the long-read
    band is doing real work via the error model, not alignment reach. This is
    still 16S bacterial amplicon on two chemistries; newer Illumina chemistries
    (i100, NovaSeq, binned qualities) remain to be tested and may move the MiSeq
    boundary. No default change is warranted beyond confirming the existing ones
    are well-placed.

---

#### Caveats!

!!! danger "Read before using these numbers"

    - **Sample sizes and marker coverage.** The single-sample outcomes are from
      one Illumina and one PacBio sample (hundreds of uniques each); the at-scale
      confirmation spans a 384-sample **MiSeq 16S** run and a 74-sample **PacBio
      HiFi 16S** run. All four point the same way, but every dataset here is 16S
      bacterial amplicon — retuning `KDIST_CUTOFF` or `BAND_SIZE` still needs
      validation across other markers (ITS, 18S), organisms, and chemistries
      before any global default change.
    - **Per-sample lower bound.** The headroom and band-fit ceilings are the
      maximum over a *single* sample's error copies. Deeper or more diverse
      samples can push real error copies further out, so the measured slack is a
      lower bound on safe tightening, not a target.
    - **Abundance dependence.** The distance at which a real error copy can
      appear scales with the abundance of its parent and sequencing depth; the
      `--nearest-parent` proxy does not run the full abundance p-value, so it
      approximates rather than reproduces DADA2's actual linkage decision.
      `--from-dada` removes this caveat — its labels *are* dada's decision — but
      then sees only the comparisons inference reached (it cannot show pairs the
      screen wrongly merged, since those never became two ASVs).
    - **Long-read saturation confounds k=5.** PacBio k=5 distances are
      compressed by saturation; use k=7 figures for long-read interpretation.
    - **Unbanded by design.** The tool aligns unbanded (`--band -1`) so distant
      pairs report true divergence; this is why it is slow on long reads, and
      why `band_req` is meaningful (it is derived from the *true* path).
    - **Subsampling.** `--max-pairs` / `--max-uniques` random-subsample to bound
      the O(n²) cost; results are statistical, and pooled inputs in particular
      should be run with a cap.
