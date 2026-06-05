# Illumina MiSeq (paired-end) workflow

An end-to-end paired-end example. The steps below follow the
[DADA2 MiSeq SOP](http://benjjneb.github.io/dada2/tutorial.html). All
intermediate outputs are JSON and can be inspected or plotted independently. Run
`dada2-rs <subcommand> --help` for the full parameter set of any step.

## 1. Filter and trim

```bash
dada2-rs filter-and-trim \
  --fwd  raw/sample_R1.fastq.gz --filt filtered/sample_R1.fastq.gz \
  --rev  raw/sample_R2.fastq.gz --filt-rev filtered/sample_R2.fastq.gz \
  --trunc-len 240 160 \
  --max-n 0 --max-ee 2 2 --trunc-q 2 \
  --compress --verbose
```

## 2. Learn the error model

**Option A — from pre-computed derep JSON files** (via the `sample` subcommand):

```bash
# Dereplicate and subsample forward reads across all samples
dada2-rs sample filtered/*_R1.fastq.gz \
  --output-dir sample_json/ --nbases 100000000 --verbose

# Learn error model
dada2-rs errors-from-sample sample_json/*.json \
  --errfun loess -o errors_fwd.json --verbose
```

**Option B — directly from FASTQ** (single command):

```bash
dada2-rs learn-errors filtered/*_R1.fastq.gz \
  --nbases 100000000 --errfun loess \
  -o errors_fwd.json --verbose
```

### Visualise cluster diagnostics during error learning

Pass `--diag-dir` to emit an `iter_NNN.json` file for each self-consistency
iteration, then plot with the bundled R script:

```bash
dada2-rs errors-from-sample sample_json/*.json \
  --errfun loess --diag-dir diag_fwd/ -o errors_fwd.json --verbose

Rscript comparison/plot_cluster_diag.R diag_fwd/ diag_fwd/cluster_diag.pdf
```

The plot shows cluster counts, birth-type breakdown, convergence trace, and
alignment work across iterations.

## 3. Denoise each sample

```bash
dada2-rs dada filtered/sample_R1.fastq.gz \
  --error-model errors_fwd.json --show-map \
  -o dada/sample_R1.json --verbose
```

Repeat for reverse reads using the reverse error model.

!!! tip "Multiple samples at once"
    `dada` accepts more than one input. Pass several filtered FASTQs and an
    `--output-dir` to denoise them in one invocation; use `--sample-jobs N` to
    control how many run concurrently. See
    [Performance & Benchmarking](benchmarking.md) for the concurrency model.

## 4. Merge paired reads

```bash
dada2-rs merge-pairs \
  --fwd-dada dada/fwd/*.json \
  --rev-dada dada/rev/*.json \
  --fwd-fastq filtered/fwd/*.fastq.gz \
  --rev-fastq filtered/rev/*.fastq.gz \
  -o merged.json --verbose
```

## 5. Build sequence table and remove chimeras

```bash
dada2-rs make-sequence-table merged.json -o seqtab.json
dada2-rs remove-bimera-denovo seqtab.json --method consensus -o seqtab_nochim.json
```

## Visualise the error model

```bash
Rscript scripts/plot_errors.R errors_fwd.json errors_fwd.pdf
```
