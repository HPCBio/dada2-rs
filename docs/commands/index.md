# Command reference

`dada2-rs <subcommand> --help` prints a one-line summary for every flag,
grouped into categories. These pages carry the detail that does not fit on one
line: what a parameter actually does, when to change it, and what the evidence
behind its default is.

## Reading `--help`

Flags are grouped under a fixed set of headings, used consistently across every
subcommand:

| Heading | What it holds |
|---|---|
| **Input** | Input files, sample naming, Phred offset |
| **Error model** | The error model file, fitting function, and its knobs |
| **Denoising** | The DADA2 algorithm parameters (R's `setDadaOpt()` surface) |
| **Pseudo-pooling** | Prior selection between the two `dada-pseudo` rounds |
| **Alignment** | Needleman-Wunsch scoring, band radius, aligner backend |
| **Screening** | The pre-alignment k-mer screen that decides which pairs are aligned |
| **Trimming** | Truncation and end-trimming of reads |
| **Filtering** | Read- or sequence-level discard rules |
| **Primers** | Primer sequences and matching tolerance |
| **Merging** | Paired-end overlap and concatenation rules |
| **Chimera** | Bimera and trimera detection parameters |
| **Classification** | Taxonomic assignment thresholds |
| **Metrics** | Optional per-read metrics (`summary`, `summary-merge`) |
| **Pooling regime** | Which population is scored, and how it is sampled (`kdist-calibrate`) |
| **Evaluation** | Truth-set classification thresholds (`reference-eval`) |
| **Performance** | Threads and concurrency; never affects results |
| **Diagnostics** | Extra output for QC and debugging; never affects ASVs |
| **Experimental** | Unstable flags. Off by default, and not drop-in replacements |

Two headings are worth reading as a pair. **Denoising** parameters define the
clusters; **Screening** parameters only decide which pairs are worth aligning
in the first place. A screen flag can change runtime by orders of magnitude and
can, if set too aggressively, cost you real ASVs — but it is not part of the
algorithm's definition of an ASV.

Anything under **Performance** or **Diagnostics** is safe to change without
affecting the ASV table.

## Pages

**Quality assessment**

- [`summary`](summary.md) — per-position quality metrics from a FASTQ file
- [`summary-merge`](summary-merge.md) — union per-sample summaries into a run report

**Read preparation**

- [`filter-and-trim`](filter-and-trim.md) — filter and trim one sample
- [`remove-primers`](remove-primers.md) — primer detection and trimming
- [`derep`](derep.md) — dereplicate a FASTQ file
- [`sample`](sample.md) — dereplicate and subsample, one JSON per sample

**Error models**

- [`learn-errors`](learn-errors.md) — fit an error model from FASTQ or JSON
- [`errors-from-sample`](errors-from-sample.md) — fit from pre-sampled JSON

**Denoising**

- [`dada`](dada.md) — per-sample (R `pool=FALSE`)
- [`dada-pooled`](dada-pooled.md) — full pooling (R `pool=TRUE`)
- [`dada-pseudo`](dada-pseudo.md) — pseudo-pooling (R `pool="pseudo"`)

**Tables and chimeras**

- [`merge-pairs`](merge-pairs.md) — merge denoised read pairs
- [`make-sequence-table`](make-sequence-table.md) — sample x sequence table
- [`remove-bimera-denovo`](remove-bimera-denovo.md) — bimera removal
- [`chimera-diagnostics`](chimera-diagnostics.md) — trimera screen

**Taxonomy**

- [`assign-taxonomy`](assign-taxonomy.md) — Naive Bayes k-mer classifier
- [`assign-species`](assign-species.md) — exact-match species assignment

**Converters**

- [`seq-table-to-tsv`](seq-table-to-tsv.md)
- [`seq-table-to-fasta`](seq-table-to-fasta.md)
- [`tax-to-tsv`](tax-to-tsv.md)

**Calibration and evaluation**

- [`kdist-calibrate`](kdist-calibrate.md) — calibrate the k-mer screen
- [`reference-eval`](reference-eval.md) — evaluate ASVs against a truth set


## See also

- [Parameters: `setDadaOpt()` parity](../parameters.md) — the R-equivalency table
- [Tuning for your data](../tuning-for-your-data.md)
- [Findings](../findings/index.md) — the measurements behind the defaults
