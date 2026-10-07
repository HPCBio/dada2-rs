# Per-read parity with R

A sequence-table comparison says a read moved between two ASVs, or a
transition count is off by one. It cannot say which read, or whose
alignment. These scripts run per-sample `dada` on both sides under one fixed
error model and compare what the table discards (#277):

- each unique's cluster, joined by sequence;
- the transition counts (R's `$trans`, our `aux.transitions`).

With the model held fixed, any difference comes from screening, alignment or
the divisive loop, not the fit.

| Script | What it does |
|---|---|
| `r_dada_maps.R` | R `dada()` per sample; writes `<s>.r.uniques.tsv`, `<s>.r.trans.tsv`, `<s>.r.dada.rds` |
| `rs_dada_maps.sh` | `derep`, then `dada --aux-outputs` on the derep JSON, so `map` indexes a file on disk |
| `compare_maps.py` | per sample: moved uniques and differing transition cells; exit 1 on any difference |

```bash
Rscript scripts/learnerrors_to_dada2rs.R err.rds err.json
Rscript dev/r_parity_maps/r_dada_maps.R err.rds out/r filt/*.fastq.gz --threads=4
THREADS=4 dev/r_parity_maps/rs_dada_maps.sh target/release/dada2-rs err.json out/rs filt/*.fastq.gz
python3 dev/r_parity_maps/compare_maps.py out/r out/rs -o diffs.tsv
```

Defaults are the PacBio settings: `BAND_SIZE=32`, k=5. Set `--band=` and
`BAND=` together for other platforms.

## First result (#277)

On 4 PacBio samples under R's pinned model, every unique went to the same
cluster, but 37 transition cells differed in 3 of 4 samples. The cause was the
gapless shortcut firing in the final-subs pass, where R never takes it
(`use_kmers = false` leaves `kodist = -1`). That guard was lost in 447c1f3.
After the fix: 0 cells.

Per-sample runs cannot reproduce a pooled difference. For that, run
`write_reference.R --pool=true --save-dada` and re-derep each FASTQ to recover
the sequences each `$map` indexes.
