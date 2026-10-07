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
    matched directories work, but mismatched naming will pair the wrong files.
    `merge-pairs` always warns when a FASTQ differs from the one the dada JSON
    records as its `source_fastq`; `--check-sample-ids` makes a mismatch fatal.

## Recommended: name each pair upstream

Give both reads of a pair the same sample name before denoising, rather than
relying on file names:

```bash
for s in sampleA sampleB; do
  dada2-rs derep filtered/fwd/${s}_R1.fastq.gz --sample-name $s -o derep/fwd/$s.json
  dada2-rs derep filtered/rev/${s}_R2.fastq.gz --sample-name $s -o derep/rev/$s.json
done
dada2-rs dada derep/fwd/*.json --error-model errors_fwd.json --output-dir dada/fwd/
dada2-rs dada derep/rev/*.json --error-model errors_rev.json --output-dir dada/rev/

dada2-rs merge-pairs --check-sample-ids \
  --fwd-dada  dada/fwd/*.json \
  --rev-dada  dada/rev/*.json \
  --fwd-fastq filtered/fwd/*_R1.fastq.gz \
  --rev-fastq filtered/rev/*_R2.fastq.gz \
  -o merged.json
```

A derep JSON carries the sample name and its source FASTQ into the dada JSON,
and `dada --output-dir` names each output after the sample, so the forward and
reverse lists line up. If you denoise FASTQ directly, set the same names with
`--sample-name` (`dada`) or `--sample-names` (`dada`, `dada-pooled`,
`dada-pseudo`).

Without explicit names, the sample name is the FASTQ file name minus its
extension, which keeps the read-direction suffix (`sampleA_R1` vs `sampleA_R2`),
so `--check-sample-ids` rejects correct input.

## Input

**`--fwd-dada`** / **`--rev-dada`** (both required) — dada JSON files.

**`--fwd-fastq`** / **`--rev-fastq`** (both required) — the FASTQ files those
dada runs came from, re-dereplicated here to recover the read → unique map.

**`--sample-names`** — comma-separated sample names, one per input set, in
`--fwd-dada` order; defaults to the `--fwd-dada` filename stems.

**`--phred-offset`** — 33 for Sanger / Illumina 1.8+, 64 for Illumina 1.3–1.7.

## Merging

**`--min-overlap`** (default 12) — minimum number of **matching** bases in the
overlap between the forward and RC(reverse) ASVs, as R's `minOverlap`.

**`--max-mismatch`** (default 0) — maximum mismatches **plus indels** allowed in
the overlap, as R's `maxMismatch`.

The overlap is found as R's `mergePairs` finds it: an unbanded, ends-free
alignment scored 1 / −64 / −64 (match / mismatch / gap) at `--max-mismatch 0`,
and 1 / −8 / −8 otherwise. These are not the denoising scores: the heavy
penalties make a short perfect overlap win over a longer imperfect one at
another offset, which matters on length-variable amplicons such as ITS
([#272](https://github.com/HPCBio/dada2-rs/issues/272)). Where the two reads
disagree at an accepted mismatch, the forward base is used; R instead takes the
base from whichever read's cluster has more error-free reads
([#274](https://github.com/HPCBio/dada2-rs/issues/274)).

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
FASTQ filenames contain the sample name as a substring; any failure is fatal.
Off by default, because it needs the samples to have been named upstream (see
[above](#recommended-name-each-pair-upstream)). The substring test is loose: `sam1`
also matches `sam10_R1.fastq.gz`.

**`--verbose`** — per-sample progress to stderr.

## See also

- [Illumina MiSeq walkthrough](../walkthroughs/walkthrough-illumina.md)
