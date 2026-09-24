# The DADA2 SOP tutorial, in dada2-rs

A step-for-step mirror of the [DADA2 MiSeq SOP
tutorial](https://benjjneb.github.io/dada2/tutorial.html), on the same data, with
the R function named at every step.

Every command below was run to produce the numbers shown, and the results were
checked against an end-to-end R run on the same data — **identical ASV sets,
abundances within 0.003%**. The dataset is the
tutorial's own: 20 samples of mouse gut V4 16S (19 timepoints plus a mock),
available from the [DADA2 tutorial
page](https://benjjneb.github.io/dada2/tutorial.html).

!!! info "Which walkthrough do you want?"
    This one follows the tutorial, so it is the right page if you already know
    DADA2 in R and want the translation. For the general shape of a paired-end
    run with your own data, see [Illumina MiSeq](walkthrough-illumina.md).

## The mapping, at a glance

| Tutorial step | R | dada2-rs |
|---|---|---|
| Inspect read quality | `plotQualityProfile()` | `summary` + [`plot_quality_profile.R`](../scripts/plotting.md#plot_quality_profiler) |
| Filter and trim | `filterAndTrim()` | `filter-and-trim` |
| Learn error rates | `learnErrors()`, `plotErrors()` | `learn-errors` + [`plot_errors.R`](../scripts/plotting.md#plot_errorsr) |
| Sample inference | `dada()` | `dada` |
| Merge paired reads | `mergePairs()` | `merge-pairs` |
| Construct sequence table | `makeSequenceTable()` | `make-sequence-table` |
| Remove chimeras | `removeBimeraDenovo()` | `remove-bimera-denovo` |
| Track reads | the tutorial's `track` table | [`track_reads.py`](../scripts/pipeline-helpers.md#track_readspy) |
| Assign taxonomy | `assignTaxonomy()` | `assign-taxonomy` |
| Species | `addSpecies()` | `assign-species` |

Set up once:

```bash
RAW=/path/to/MiSeq_SOP
mkdir -p filtered/fwd filtered/rev dada_fwd dada_rev
```

## 1. Inspect read quality

The tutorial opens here, and so does this: look at the quality profile *before*
choosing truncation lengths.

```bash
dada2-rs summary "$RAW/F3D0_S188_L001_R1_001.fastq.gz" \
    --sample-name F3D0_R1 --report -o summary/F3D0_R1.json
```

```text
=== FASTQ quality metrics: F3D0_R1 ===
total sequences:   7793
read length:       251 cycles
quality range:     Q2–Q39 (27 distinct values)
binned:            no
```

`--report` is the quick look; the JSON is what gets plotted:

```bash
Rscript scripts/plot_quality_profile.R --out=quality_fwd.pdf summary/*_R1.json
```

`binned: no` matters more than it looks — on a binned-quality instrument it
would say so, and that changes which error model you should use. See
[binned quality scores](../findings/binned-quality.md).

## 2. Filter and trim

R's `filterAndTrim()` takes vectors of files; `filter-and-trim` handles one
sample, so this is a loop. **Name the outputs by sample, not by read direction**
— put forward and reverse in separate directories. Step 8 depends on it, and
getting it wrong is not obvious until then.

```bash
for r1 in "$RAW"/*_R1_001.fastq.gz; do
    s=$(basename "$r1" _L001_R1_001.fastq.gz)
    dada2-rs filter-and-trim \
        --fwd "$r1" --rev "${r1/_R1_/_R2_}" \
        --filt "filtered/fwd/${s}.fastq.gz" --filt-rev "filtered/rev/${s}.fastq.gz" \
        --sample-name "$s" \
        --trunc-len 240 160 --max-n 0 --max-ee 2 2 --trunc-q 2 \
        --phix-genome data/dada2/phix_genome.fa \
        --compress -o "filtered/${s}.filter.json"
done
```

Same parameters as the tutorial: `truncLen=c(240,160)`, `maxN=0`,
`maxEE=c(2,2)`, `truncQ=2`, `rm.phix=TRUE`.

**152,360 reads in → 139,642 out (91.7%).**

## 3. Learn the error rates

```bash
dada2-rs learn-errors filtered/fwd/*.fastq.gz --nbases 100000000 --threads 10 -o errors_fwd.json
dada2-rs learn-errors filtered/rev/*.fastq.gz --nbases 100000000 --threads 10 -o errors_rev.json
```

Converged in **5 self-consistency rounds forward (11.4 s)** and **6 reverse
(7.7 s)** on 10 threads. `--nbases 1e8` is R's default — and worth knowing that
[which samples get drawn at that budget matters about as much as the budget
itself](../findings/learn-errors-nbases-convergence.md).

The equivalent of `plotErrors(errF, nominalQ=TRUE)`:

```bash
Rscript scripts/plot_errors.R errors_fwd.json errors_fwd.pdf
```

## 4. Sample inference

R loops `dada()` over samples; here one invocation takes all of them.

```bash
dada2-rs dada filtered/fwd/*.fastq.gz --error-model errors_fwd.json \
    --output-dir dada_fwd/ --threads 10 --sample-jobs 3
dada2-rs dada filtered/rev/*.fastq.gz --error-model errors_rev.json \
    --output-dir dada_rev/ --threads 10 --sample-jobs 3
```

**2.8 s forward, 1.6 s reverse.** `--sample-jobs` sets how many samples run at
once, which is the knob that matters here rather than `--threads` alone — see
[threading](../findings/threading-serial-steps.md).

## 5. Merge paired reads

```bash
dada2-rs merge-pairs \
    --fwd-dada dada_fwd/*.json --rev-dada dada_rev/*.json \
    --fwd-fastq filtered/fwd/*.fastq.gz --rev-fastq filtered/rev/*.fastq.gz \
    --threads 10 -o merged.json
```

## 6 and 7. Sequence table, then chimeras

```bash
dada2-rs make-sequence-table merged.json -o seqtab.json
dada2-rs remove-bimera-denovo seqtab.json --method consensus --threads 10 \
    -o seqtab.nochim.json
```

| | ASVs | reads |
|---|---:|---:|
| `seqtab.json` | **293** | 128,895 |
| `seqtab.nochim.json` | **232** | 124,245 |

**96.4% of reads survive chimera removal** while 61 of 293 ASVs are removed —
the tutorial's point that chimeras are a large fraction of unique sequences and
a small fraction of reads.

### Checked against R, on the same data

R DADA2 1.40.0 was run end to end on these 20 samples with the same parameters:

| | R | dada2-rs |
|---|---:|---:|
| reads in → filtered | 152,360 → 139,642 | **identical** |
| sequence table | 293 ASVs | **293** |
| after chimera removal | 232 ASVs | **232** |
| ASV sequences | — | **232 shared, 0 unique to either side** |
| total reads retained | 124,249 | 124,245 (−0.0032%) |

The ASV *sets* are identical. Abundance differs on 6 of the 232, by at most 6
reads.

### Where that last 0.003% comes from — and how to remove it

All of it is the **LOESS fitting surface**, and nothing downstream of the error
model. Three arms, same filtered reads, same everything else:

| arm | ASVs vs R | abundances differing | read delta |
|---|---|---:|---:|
| default (`--loess-preset default`, direct surface) | 232, none unique either way | 6 / 232 | −4 |
| **`--loess-preset r-dada2`** (R's interpolate surface) | 232 | **0** | **0** |
| R's own error model, via [`learnerrors_to_dada2rs.R`](../scripts/pipeline-helpers.md#learnerrors_to_dada2rsr) | 232 | **0** | **0** |

So **`--loess-preset r-dada2` reproduces R's chimera-filtered table exactly on
this dataset** — every ASV, every count — and handing our pipeline R's own
fitted model does no better.

!!! note "Exact here, not exact everywhere"
    Reaching zero on this dataset is a best case. Earlier work on other data
    found small count-level differences persisting *even when supplied with R's
    exact error models* — much closer than without, but not identical. The SOP is
    20 samples of clean, deeply-characterised data.

    There is also a limit to what "matching R" can mean at all. R DADA2's output
    shifts slightly across package versions, R releases and platforms (Linux
    against macOS), so exact agreement is a statement about *one* R build on
    *one* machine, not a property two implementations can hold in general.
    Reproducing a specific published table means pinning the version and the
    environment, in R as much as here.

Note what is *not* a factor here. `--nbases 1e8` against ~33 Mbases of forward
reads means neither tool subsamples, so both fit on identical input; the
difference is purely how the surface is fitted. On a run large enough to
trigger subsampling, [which reads get drawn becomes its own
term](../findings/learn-errors-nbases-convergence.md).

The background is [the LOESS
page](../findings/loess-error-model-correctness.md): our native LOESS matches R's
`loess(surface = "direct")` to machine precision, while R DADA2 takes R's
default `surface = "interpolate"`. The default here stays `direct`; `r-dada2` is
for when bit-parity with R is what you want.

## 8. Track reads through the pipeline

The tutorial's sanity table — where a sample loses reads, if it does.

```bash
dada2-rs make-sequence-table dada_fwd/*.json -o seqtab_R1.json
dada2-rs make-sequence-table dada_rev/*.json -o seqtab_R2.json

python3 scripts/track_reads.py \
    -f filtered/*.filter.json \
    -d seqtab_R1.json seqtab_R2.json \
    -m merged.json -s seqtab.nochim.json -o track.tsv
```

```text
sample       input  filtered  denoisedF  denoisedR  merged  nochim
F3D0_S188    7793   7113      7113       7113       6540    6528
F3D141_S207  5958   5463      5463       5463       4987    4864
F3D142_S208  3183   2914      2914       2914       2595    2521
...
TOTAL        152360 139642    139642     139642     128895  124245
```

!!! bug "The denoised columns do not match R yet"
    `denoisedF` and `denoisedR` above equal `filtered`, because our `dada`
    assigns every read to a cluster. R attributes ~1.6–1.9% fewer — for F3D0 it
    reports 6,976 / 6,979 against our 7,113, and R's are the figures the
    tutorial publishes. R produces the *same clusters* (128 for F3D0 forward, as
    we do) and simply declines to place some reads in them.

    It does not propagate: `merged` differs by 11 reads across all 20 samples
    and `nochim` by 4, and the ASV set is identical. But it means these two
    columns cannot currently do the job they exist for, which is showing a
    sample that loses reads at denoising. Tracked in
    [#204](https://github.com/HPCBio/dada2-rs/issues/204).

!!! warning "One row per sample depends on step 2's naming"
    `track_reads.py` joins the steps on sample name. If the filtered files carry
    an `_R1` / `_R2` suffix, the denoised outputs inherit it and each sample
    splits into three rows that never join — totals still correct, per-sample
    rows useless. Naming the filtered files `<sample>.fastq.gz` in separate
    `fwd/` and `rev/` directories is what keeps the join intact.

## 9. Assign taxonomy

```bash
dada2-rs assign-taxonomy seqtab.nochim.json \
    --ref-fasta silva_nr99_v138.2_toGenus_trainset.fa.gz \
    --threads 10 -o taxonomy.json

dada2-rs assign-species taxonomy.json \
    --ref-fasta silva_v138.2_assignSpecies.fa.gz -o taxonomy_species.json
```

| rank | assigned |
|---|---|
| Kingdom – Order | 232/232 (100%) |
| Family | 210/232 (91%) |
| Genus | 127/232 (55%) |
| Species | 16/232 (7%) |

Taxonomy took **41 s**, species **4.8 s** (352,047 reference sequences).

Two things worth knowing:

- **`assign-species` implements R's `addSpecies()`, not `assignSpecies()`.** It
  applies the genus-consistency check, so an ASV whose classifier genus
  disagrees with the binomial genus gets no species. Raw `assignSpecies` returns
  17 here; the check drops one, giving 16 — identical to R's `addSpecies`, ASV
  for ASV.
- **The species step is an Illumina-workflow habit.** PacBio runs skip it.

`--seed` has a fixed default, so these assignments are reproducible without
setting it; vary it to see which calls sit near `--min-boot`.

## Exporting

```bash
dada2-rs seq-table-to-tsv seqtab.nochim.json -o seqtab.tsv
dada2-rs tax-to-tsv taxonomy_species.json -o taxonomy.tsv
```

From there, [`examples/import_mia.R`](../scripts/pipeline-helpers.md#examplesimport_miar)
loads both into a `TreeSummarizedExperiment`.

## What is different from R, and why

- **`filter-and-trim` is one sample per call**, where `filterAndTrim()` is
  vectorised. Hence the loop.
- **`remove-primers` has no separate step here.** The SOP does not use primers;
  a run that needs them fuses primer removal and filtering into one pass.
- **Intermediates are JSON files, not R objects.** Every step is inspectable and
  restartable, at the cost of naming things carefully — see step 8.
- **The map from reads to clusters is always emitted.** R's `dada()` carries it
  too; here it is unconditional because downstream tools rely on it.
