# PacBio HiFi (single-end) workflow

PacBio HiFi 16S reads are long (~1.5 kb), single-end, and dominated by
homopolymer indels rather than substitution errors. The workflow therefore
differs from [Illumina](walkthrough-illumina.md) in a few key ways:

- **Primers are removed, reads oriented, and quality/length filtered in one
  step** — see step 1 below.
- **No reverse reads and no merge step** — denoised ASVs go straight into the
  table.
- **`--errfun pacbio`** is used instead of `loess`.
- The `dada` step uses PacBio-tuned alignment parameters: a wider band
  (`--band 32`), a relaxed homopolymer gap penalty (`--homo-gap-p -1`), and a
  larger k-mer screen (`--kmer-size 7`).

!!! warning "Do not use k=5 on HiFi reads"
    On ~1.5 kb reads the Illumina default `--kmer-size 5` makes the
    pre-alignment screen a no-op (nearly every pair is fully aligned), so
    denoising runs several times slower. Use `--kmer-size 7`.

!!! important "Keep learn and denoise params consistent"
    Pass the same `--band`, `--homo-gap-p`, and `--kmer-size` to **both**
    `learn-errors` and `dada`. If the error model is learned with different
    alignment params than denoising uses, `dada` will warn, and results change
    subtly. (See the param-mismatch note in
    [Performance & Benchmarking](benchmarking.md#4-built-in-instrumentation-the-binarys-own-logs).)

## 1. Remove primers, orient, and filter

HiFi reads come in both orientations and still carry the amplification primers.
`remove-primers` (mirrors R's `removePrimers()`) trims them and, with `--orient`
(on by default), flips reads that match only in the reverse-complement direction
so every read points the same way — essential for PacBio, since the k-mer screen
in later steps relies on consistently-oriented reads. It also accepts the same
length/quality filters as `filter-and-trim`, so primer removal and filtering are
done in a single pass (no separate `filter-and-trim` step needed):

```bash
dada2-rs remove-primers raw/sample.fastq.gz \
  --fout filtered/sample.fastq.gz \
  --primer-fwd AGRGTTYGATYMTGGCTCAG \
  --primer-rev RGYTACCTTGTTACGACTT \
  --trim-fwd --trim-rev --orient \
  --min-len 1000 --max-len 1600 \
  --max-n 0 --max-ee 2 --trunc-q 0 \
  --compress -o filtered/sample.json
```

The `--primer-rev` sequence is given 5′→3′ (catalog direction) and is
reverse-complemented automatically before matching. Reads with no primer match
are discarded.

!!! note "Alternative: cutadapt"
    This step consolidates the `removePrimers` / `filterAndTrim` steps into one.
    You can also use a tool like `cutadapt`, which trims, filters, and reorients
    in one pass. We see slightly more reads passing filter with `cutadapt` at the
    same general parameters, likely due to more flexible primer matching.

## 2. Learn the error model

```bash
dada2-rs learn-errors filtered/*.fastq.gz \
  --nbases 100000000 --errfun pacbio --band 32 --homo-gap-p -1 --kmer-size 7 \
  -o errors_pacbio.json --verbose
```

## 3. Denoise each sample

```bash
dada2-rs dada filtered/sample.fastq.gz \
  --error-model errors_pacbio.json \
  --band 32 --homo-gap-p -1 --kmer-size 7 \
  -o dada/sample.json --verbose
```

Alternatively, store the PacBio parameters in the error-model JSON and pass
`--inherit-err-params` so every `dada` call picks them up automatically.

## 4. Build sequence table and remove chimeras

There is no merge step — feed the per-sample `dada` JSON files directly into the
table:

```bash
dada2-rs make-sequence-table dada/*.json -o seqtab.json
dada2-rs remove-bimera-denovo seqtab.json --method consensus -o seqtab_nochim.json
```
