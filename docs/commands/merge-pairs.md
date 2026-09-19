# `merge-pairs`

Merge denoised forward and reverse reads into full-length amplicons.

```bash
dada2-rs merge-pairs \
  --fwd-dada  fwd_dada/*.json \
  --rev-dada  rev_dada/*.json \
  --fwd-fastq fwd_fastq/*.fastq.gz \
  --rev-fastq rev_fastq/*.fastq.gz \
  -o merged.json
```

For each sample the forward and reverse FASTQ files are re-dereplicated to
reconstruct the read → unique mapping, which is composed with the unique → ASV
mapping from the dada JSONs to count every (forward ASV, reverse ASV) pair. Each
distinct pair is then aligned — ends-free Needleman-Wunsch of the forward ASV
against the reverse complement of the reverse ASV — and accepted or rejected on
overlap length, mismatches and indels.

!!! warning "Files are matched by position"
    The first `--fwd-dada` corresponds to the first `--rev-dada`,
    `--fwd-fastq` and `--rev-fastq`. Shell globbing sorts consistently, so
    matched directories work, but mismatched naming will silently pair the
    wrong files. Pass `--check-sample-ids` to have that verified.

## Input

**`--fwd-dada`** / **`--rev-dada`** (both required) — dada JSON files.

**`--fwd-fastq`** / **`--rev-fastq`** (both required) — the FASTQ files those
dada runs came from, re-dereplicated here to recover the read → unique map.

**`--sample-names`** — override sample names; defaults to the `--fwd-dada`
filename stems.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Merging

**`--min-overlap`** (default 12) — minimum overlap between the forward and
RC(reverse) ASVs.

**`--max-mismatch`** (default 0) — maximum mismatches allowed in the overlap
region.

**`--just-concatenate`** — concatenate forward and RC(reverse) with an N spacer
instead of merging.

**`--rescue-unmerged`** — concatenate pairs that *fail* to merge rather than
dropping them. Useful for variable-length amplicons such as ITS, where the reads
may genuinely not overlap. Rescued reads are marked `concatenated: true`, and
this takes precedence over `--return-rejects`.

**`--concat-nnn-len`** (default 10) — number of `N` characters in the
concatenation spacer.

**`--trim-overhang`** — trim overhanging portions of the reads past the overlap.

## Performance

**`--threads`** (default 1) — used within each sample for dereplication.

## Output

**`--output` / `-o`** — write JSON here instead of stdout.

**`--return-rejects`** — include rejected merges, with `accept: false`, in the
output.

**`--compact`** — minified JSON.

## Diagnostics

**`--check-sample-ids`** — verify that the forward and reverse dada JSONs carry
the same `sample` field, that it matches the resolved sample name, and that both
FASTQ filenames contain the sample name as a substring. Cheap insurance against
positional mismatching; worth using on every multi-sample run.

**`--verbose`** — per-sample progress to stderr.

## See also

- [Illumina MiSeq walkthrough](../walkthroughs/walkthrough-illumina.md)
