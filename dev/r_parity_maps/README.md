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
| `r_pooled_maps.R` | the same uniques table from a saved pooled run (`write_reference.R --save-dada`) |
| `rs_pooled_maps.sh` | lays out a `dada-pooled` output directory for `compare_maps.py` |
| `export_err_in.R` | writes each learnErrors pass's input model (`$err_in[[k]]`) as `err_in_<k>.json` |
| `run_replay.sh` | replays every learnErrors pass on both sides from R's per-pass inputs |

Both sides check that map + derep rebuild every ASV's reads, so a derep in a
different order than the run's own fails rather than scrambling the join.

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

This path feeds `dada --aux-outputs` and cluster quality only. learn-errors
builds its transitions separately (`build_trans_mat`, Raws without k-mers), so
it never took the shortcut, and the learned model is unchanged.

## Replaying learnErrors

A learned model that differs from R's says nothing about which pass first
diverged. `learnErrors` keeps the input model of every self-consistency pass, so
each pass can be rerun on both sides from identical input:

```bash
THREADS=8 dev/r_parity_maps/run_replay.sh target/release/dada2-rs ref.err.rds replay train/*.fastq.gz
```

The FASTQs must be exactly the training set, all of it used. Each pass runs
per-sample `dada` with `OMEGA_C = 0`, as `learnErrors` does. At the last pass,
R's per-sample transitions must sum to the model's own `trans`. If they do not,
the run is not replaying `learnErrors`.

- **Every pass identical:** denoising is the same given the same model, so the
  learned models drift through the error-model fit.
- **A pass differs:** `pass_<k>/diffs.tsv` names the samples, uniques and cells.

Both ends must be exact for this to mean anything: the converter writes 17
significant digits and `serde_json` parses with `float_roundtrip`.

## Tracing one read

When one read lands in a different cluster, `trace_read.sh` follows it through
both implementations. A patched R DADA2 (`r_trace.patch`) and dada2-rs
(`DADA2RS_TRACE_RAW`) print the same `TRACE` lines at the same points: every
cluster birth, the read's comparison with each new cluster (λ, stored or not,
lock, `E_minmax`), and its cluster after each shuffle and p-update. Doubles are
printed as raw IEEE-754 bits, so equal means bit-equal. The first differing
line names the quantity.

```bash
dev/r_parity_maps/install_r_trace.sh ~/R/dada2-trace          # once, with the reference's R
R_TRACE_LIB=~/R/dada2-trace ITER=6 dev/r_parity_maps/trace_read.sh \
    target/release/dada2-rs ref.err.rds sample.fastq.gz <read sequence> trace_out
```

Build the traced DADA2 with the same R, on the same machine, as the reference
being compared: compiler flags are part of what is under test. Birth p-values
differ in their last bits unless dada2-rs is built with `--features
rmath-ppois`.

## Pooled runs

Per-sample runs cannot reproduce a pooled difference. For a pooled run:

```bash
Rscript dev/r_parity_maps/r_pooled_maps.R ref.dd.rds filt/ out/r
dev/r_parity_maps/rs_pooled_maps.sh target/release/dada2-rs run/dada out/rs filt/*.fastq.gz
python3 dev/r_parity_maps/compare_maps.py out/r out/rs -o moved.tsv
```

On 4 local PacBio samples pooled under R's model: 0 moved uniques. As a
negative control, R pooled against our per-sample maps reports 305-733 moved
uniques per sample.
