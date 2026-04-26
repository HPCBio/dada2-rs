# Cluster trace dumps

`dada2-rs` can optionally write a full description of cluster structure
to JSON files for off-line plotting and inspection. This is finer-grained
than `--diag-dir` (which only emits cluster counts and birth-type
breakdowns) — every cluster's center, members, hamming distances, λ
values, and final p-values are dumped.

## Producing traces

For learning runs, write one file per (sample × iteration):

```bash
dada2-rs learn-errors *.fastq.gz \
    --cluster-trace-dir traces/ \
    --trace-min-abund 2 \
    -o err.json
```

For one-off `dada` ASV calls, write a single trace for the final result:

```bash
dada2-rs dada sample.fastq.gz \
    --error-model err.json \
    --cluster-trace traces/sample1_clusters.json
```

Knobs (apply to both subcommands):

| Flag | Effect |
| --- | --- |
| `--trace-no-members` | Skip the `members` array; emit only cluster centers and birth metadata. ~10× smaller files. |
| `--trace-min-abund N` | When `members` is included, drop members with abundance < N. Cluster's full member count is still recorded as `n_members`; the dropped count is in `members_excluded`. |

## File layout

Each file is tagged JSON:

```json
{
  "dada2_rs_command": "cluster-trace",
  "dada2_rs_version": "0.1.0",
  "sample": "sample_001",
  "iteration": 3,
  "nclust": 142,
  "total_reads": 7228,
  "trace_no_members": false,
  "trace_min_abund": 1,
  "nq": 40,
  "err_in": [[...]],
  "sequences": ["ACGT...", ...],
  "clusters": [
    {
      "id": 0,
      "center_seq_id": 12,
      "abundance": 4231,
      "n_members": 318,
      "birth_type": "Initial",
      "birth_from": 0,
      "birth_hamming": 0,
      "members_excluded": 0,
      "members": [
        { "raw_seq_id": 12, "abundance": 1234, "hamming": 0,
          "lambda": 1.0, "e_reads": 4231.0, "pval": 1.0 }
      ]
    }
  ]
}
```

Centers and members reference sequences by integer index into the
top-level `sequences` array (deduplicated). `iteration` is omitted for
the single-shot `dada` trace.

## Plotting examples

| File | Language | Notes |
| --- | --- | --- |
| `cluster_size_dist.R` | R + ggplot2 + jsonlite | Cluster abundance vs rank, faceted per file, colored by birth type. Useful for spotting convergence across iterations. |
| `lambda_distribution.R` | R + ggplot2 + jsonlite | Member λ vs hamming, jittered, sized by abundance. Reveals which members sit near the acceptance boundary. |
| `birth_genealogy.py` | Python (stdlib only) | Renders the cluster birth tree as a Graphviz DOT file: which cluster spawned which, with edge color = birth p-value and edge label = birth hamming distance. Pipe through `dot -Tpdf`. |

Invoke them as:

```bash
Rscript examples/cluster_trace/cluster_size_dist.R \
    out.pdf traces/cluster_iter_*.json

python3 examples/cluster_trace/birth_genealogy.py traces/cluster_iter_007_sample_001.json \
    | dot -Tpdf -o tree.pdf
```

## Cost

Per-file size scales with `n_members`. For a typical Illumina sample
with a few thousand uniques across ~150 final clusters, expect ~100 KB
per file with default flags, dropping to ~10 KB with
`--trace-no-members`. With `--trace-min-abund 2` you'll cut singleton
noise without losing structural detail.

The dump itself is cheap (one disk write per sample-iter); the algorithm
is unchanged. Disabling the trace (omitting the flag) bypasses the
writer entirely.
