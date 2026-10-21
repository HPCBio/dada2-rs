# `remove-primers`

Remove primer sequences from a FASTQ file. Mirrors R's `removePrimers()`.

```bash
dada2-rs remove-primers reads.fastq.gz -f trimmed.fastq.gz \
  --primer-fwd AGRGTTYGATYMTGGCTCAG \
  --primer-rev RGYTACCTTGTTACGACTT
```

Detects and trims forward (and optionally reverse) primers from each read using
mismatch-tolerant, IUPAC-aware matching. **Reads lacking a primer match are
discarded.** With `--orient` (on by default), reads that match primers only in
the reverse-complement direction are flipped before trimming.

Outputs a trimmed FASTQ; JSON stats (`reads_in` / `reads_out`) go to stdout or
`-o`.

As a port of R's `removePrimers()`, `remove-primers` is currently designed for
single-end reads, mainly PacBio HiFi. Run on each file of a pair, it would
discard different reads from R1 and R2 and break the pairing. For paired-end
data use `cutadapt`; see
[Removing primers from paired-end reads](../primer-removal.md).

!!! tip "Verify the primers are actually there"
    Before trusting any primer-removal result on data you did not prepare,
    grep the raw reads for the primer — and for the reverse complement of the
    reverse primer at the 3′ end. Metadata is not evidence, and a run that
    removes 0% or 100% of reads usually means the primer string, its
    orientation, or `--rc-primer-rev` is wrong.

## Input

**`<INPUT>`** — input FASTQ, uncompressed or gzipped.

**`--sample-name`** — sample identifier for the output JSON. Defaults to the
input filename stem.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.
Only relevant when the quality-based filters below are used.

## Primers

**`--primer-fwd`** (required) — forward primer in its 5′→3′ (catalog /
synthesis) direction. IUPAC ambiguity codes are accepted, e.g.
`AGRGTTYGATYMTGGCTCAG`.

**`--primer-rev`** — reverse primer, also 5′→3′. Omit to skip reverse primer
detection.

**`--rc-primer-rev`** (default true) — reverse-complement `--primer-rev` before
matching. Primers are conventionally written 5′→3′, and the reverse primer must
be reverse-complemented to match the orientation it actually appears in reads.
Pass `--rc-primer-rev=false` only when you are supplying `--primer-rev` already
reverse-complemented.

!!! warning "The opposite of R's convention"
    R's `removePrimers()` expects `primer.rev` already reverse-complemented
    (its PacBio tutorial passes `rc(R1492)`). Here `--primer-rev` takes the
    catalog orientation. A primer copied from an R script is therefore
    reverse-complemented twice, matches nothing, and every read is discarded
    as primerless; add `--rc-primer-rev=false`, or give the catalog sequence
    ([#296](https://github.com/HPCBio/dada2-rs/issues/296)).

**`--max-mismatch`** (default 2) — maximum mismatches allowed when matching each
primer.

**`--allow-indels`** — also allow insertions and deletions when matching, using
Levenshtein edit distance where each mismatch or indel counts 1 toward
`--max-mismatch`. Significantly slower than the default mismatch-only mode.

**`--trim-fwd`** / **`--trim-rev`** (both default true) — trim the matched
primer from the 5′ / 3′ end.

**`--orient`** (default true) — detect and correct read orientation: reads that
match primers only in the reverse complement are flipped before trimming.

These four flags, and `--compress` below, are on by default as in R. The bare
flag keeps them on; turn one off with `=false`, e.g. `--orient=false`. The `=`
is required, so `--orient false` is an error.

## Trimming and Filtering

These are applied *after* primer trimming and have the same meaning as in
[`filter-and-trim`](filter-and-trim.md), except that here each takes a single
value rather than a forward/reverse pair:

`--trunc-q`, `--trunc-len`, `--trim-left`, `--trim-right`, `--max-len`,
`--min-len`, `--max-n`, `--min-q`, `--max-ee`, `--phix-genome`,
`--rm-lowcomplex`.

All are off by default here, unlike `filter-and-trim`.

## Performance

**`--threads`** (default 1) — threads for primer matching and bgzf output
compression. Values above 1 switch the output to bgzf.

## Output

**`--fout` / `-f`** (required) — output FASTQ file.

**`--compress`** (default true) — gzip-compress the output FASTQ;
`--compress=false` writes plain FASTQ.

**`--output` / `-o`** — write JSON stats here instead of stdout.

**`--compact`** — minified JSON.

## Diagnostics

**`--verbose`** — progress to stderr.

## See also

- [PacBio HiFi walkthrough](../walkthroughs/walkthrough-pacbio.md) — primer
  removal in context
- [Reading the prep first](../findings/reading-the-prep.md)
