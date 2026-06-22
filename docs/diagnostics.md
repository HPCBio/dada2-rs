# Diagnostics

This page documents **experimental diagnostic tooling** built into `dada2-rs`, along with related helper scripts, for interrogating the internals of the DADA2 algorithm. 

---

## `kdist-calibrate`: k-mer screen, `KDIST_CUTOFF`, `BAND_SIZE`

Currently the k-mer **screen** (`KDIST_CUTOFF`) and the alignment **band size** 
(`BAND_SIZE`) are two constants that ship with set defaults or recommendations, but without much documentation on how these were derived (which data were used, the tools for calibration, etc). Some of the tooling, e.g., ESPRIT, also are no longer available online. 

!!! warning "Hidden, experimental subcommand"

    `kdist-calibrate` is a **hidden** subcommand: it does not appear in
    `dada2-rs --help` and is **not part of the supported pipeline**. It exists
    to *measure* and *characterize* algorithm behaviour, not to denoise. Its
    flags, output columns, and defaults may change without notice. Treat its
    results as exploratory — see the [caveats](#caveats) before drawing
    conclusions, especially the small sample sizes behind the preliminary
    numbers below.


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

This mode is also far cheaper — O(n²) *cheap* k-mer comparisons to find each
parent plus only O(n) alignments (vs O(n²) alignments in all-pairs mode) — so it
scales to pooled/multisample inputs much better.

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

!!! danger "Preliminary — tiny datasets"

    All numbers below come from **single, small samples** (one Illumina sample,
    449 uniques; one PacBio sample, 259 uniques). They illustrate the *kind* of
    result the tool produces; they are **not** a basis for retuning any default.
    See [caveats](#caveats).

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

---

#### Caveats!

!!! danger "Read before using these numbers"

    - **Tiny sample sizes.** The preliminary outcomes are from one Illumina and
      one PacBio sample (hundreds of uniques each). They demonstrate the
      *method*, not population-level truth. Any retuning of `KDIST_CUTOFF` or
      `BAND_SIZE` must be validated across **many** samples spanning depth,
      diversity, and chemistry.
    - **Per-sample lower bound.** The headroom and band-fit ceilings are the
      maximum over a *single* sample's error copies. Deeper or more diverse
      samples can push real error copies further out, so the measured slack is a
      lower bound on safe tightening, not a target.
    - **Abundance dependence.** The distance at which a real error copy can
      appear scales with the abundance of its parent and sequencing depth; the
      `--nearest-parent` proxy does not run the full abundance p-value, so it
      approximates rather than reproduces DADA2's actual linkage decision.
    - **Long-read saturation confounds k=5.** PacBio k=5 distances are
      compressed by saturation; use k=7 figures for long-read interpretation.
    - **Unbanded by design.** The tool aligns unbanded (`--band -1`) so distant
      pairs report true divergence; this is why it is slow on long reads, and
      why `band_req` is meaningful (it is derived from the *true* path).
    - **Subsampling.** `--max-pairs` / `--max-uniques` random-subsample to bound
      the O(n²) cost; results are statistical, and pooled inputs in particular
      should be run with a cap.
