# Dereplicate per sample, keep JSON

**Feeding `dada-pooled` pre-dereplicated JSON instead of FASTQ gives identical
output and cuts the serial load front by 44–90%.** On a full MiSeq run that is
25–31% of pooled wall time; on NovaSeq and PacBio, where denoising dominates, it
is 1.5–4%. Peak RSS also falls 14–30%.

What is left for a cheaper intermediate format is the JSON parse: **at most 3.9%
of pooled wall** on any dataset measured. That closes
[#1](https://github.com/HPCBio/dada2-rs/issues/1) (move off JSON) as not worth
the change. The real cost on the write side is gzip at level 6, which is most of
a per-sample `derep` run, tracked in
[#331](https://github.com/HPCBio/dada2-rs/issues/331).

## The question

A workflow that dereplicates each sample in its own job, then feeds the JSON to
`learn-errors`, `dada`, `dada-pseudo` and `dada-pooled`, moves FASTQ
dereplication out of the pooled run's serial load front. The concern was whether
JSON on the read side would give that back, and whether this was the point to
change format.

## Setup

[`dev/derep_ab.sh`](https://github.com/HPCBio/dada2-rs/blob/191f8ff/dev/derep_ab.sh)
runs `dada-pooled` on the FASTQs, then `derep` on each FASTQ to compact
`.json.gz` (8 processes at a time), then `dada-pooled` on the JSONs, with the
same input order and sample names. Per-sample outputs must match field for
field, ignoring only `input_file` and the version string.

Cluster node compute-5-6, 48 threads, `numactl --interleave=all`
([why](measuring-on-numa.md)), release-native build at 795e3cc, one replicate
per arm. PacBio used `--band 32 --kmer-size 7`.

## Results

All five runs: **outputs identical.**

| dataset | wall, FASTQ → JSON | load front, FASTQ → JSON | peak RSS |
|---|---|---|---:|
| MiSeq F, 362 samples | 99.2 → 74.6 s (−24.7%) | 27.9 → 3.5 s | −20% |
| MiSeq R | 79.8 → 55.2 s (−30.9%) | 26.9 → 2.6 s | −21% |
| NovaSeq 16S F, 30 samples | 702.6 → 675.1 s | 16.2 → 5.6 s | −30% |
| NovaSeq 16S R | 626.2 → 627.6 s | 16.1 → 6.1 s | −27% |
| PacBio, 95 samples | 403.7 → 386.9 s (−4.1%) | 37.9 → 21.1 s | −14% |

**The NovaSeq wall deltas are noise.** The denoising phase, unchanged between
arms, differed by ±2% (683 vs 666 s forward, 606 vs 618 s reverse), more than the
10 s saved in loading. The load-front numbers are well clear of it.

**The peak RSS drop is not explained.** Both arms end with the same pooled
uniques, so the difference is in how the FASTQ path builds them. It is one
replicate, but RSS is usually stable between runs; treat it as real until
measured otherwise.

## What a format change could still buy

With JSON input, the load front splits into read + gunzip and parse:

| dataset | read + gunzip | parse | parse, share of wall |
|---|---:|---:|---:|
| MiSeq F | 1.2 s | 2.3 s | 3.0% |
| NovaSeq 16S F | 1.4 s | 4.0 s | 0.6% |
| PacBio | 5.8 s | 15.0 s | 3.9% |

A binary format removes most of the parse and some of the read. The PacBio parse
matches the post-fix figure in [the derep load page](derep-load-format-bound.md)
(the same 4.95 GB of JSON), so the ceiling held across the cluster migration.
**Under 4% of wall is not worth giving up files that can be opened and read**,
which several findings here depended on.

## The cost moved to the write side

Running `derep` as separate processes cost 3–5× the CPU of dereplicating inside
`dada-pooled`:

| dataset | in-process derep | sum of per-sample `derep` | wall, 8 at a time |
|---|---:|---:|---:|
| MiSeq F | 27.8 s | 86.9 s | 19.5 s |
| NovaSeq 16S F | 15.6 s | 66.6 s | 11.5 s |
| PacBio | 37.6 s | 179.9 s | 30.6 s |

On one NovaSeq sample (laptop, 234 MB of uncompressed JSON, 3 replicates),
`derep` takes 2.4 s writing plain `.json` and 6.3 s writing `.json.gz`: **gzip at
level 6 is ~60% of the run.** Uncompressed is not the fix; it is 15× larger, and
a workflow stages these files between jobs. Compressed, the derep JSON is smaller
than the FASTQ it came from (PacBio 339 vs 826 MB; NovaSeq F 101 vs 212 MB).
Faster levels and backends are assessed in
[#331](https://github.com/HPCBio/dada2-rs/issues/331).

## What this dictates

- **Dereplicate once per sample in workflows.** The derep JSON is produced once
  and read by every downstream command, and the pooled run no longer
  dereplicates serially.
- **Keep JSON for intermediates.** Do not reopen #1 unless parse becomes a large
  share of wall at some scale not measured here.
- **Fix gzip write speed before anything else on I/O.** It is the largest I/O
  cost found, and a compression level is a smaller change than a format.
