# External indel-rate fits

`dada2-rs learn-errors --indel-rates-out` learns indel rates for `dada`'s
experimental indel-aware λ ([#291](https://github.com/HPCBio/dada2-rs/issues/291)).
By default each iteration fits its rates natively, pooling each quality with
its neighbours until it holds 20 events. `--indel-fit-cmd "<command>"` hands
the fit to a command instead, the way `--errfun-cmd` does for substitutions
(see [`../external_errfun`](../external_errfun/README.md)).

## Wire format

`dada2-rs` runs `<command> <counts-tsv-path> <rates-tsv-path>` once per
indel-learning iteration. The command string is whitespace-split into argv (no
shell interpolation).

The counts TSV is the `dada --indel-opportunities` layout: a header
`sample  table  key  opportunities  events`, with one sample, `learned`. The
command must write the rates TSV with a header `table  key  rate` and these
tables:

| table | key | rate |
| --- | --- | --- |
| `ins_q_novel` | quality of the inserted base | insertion rate per base, for bases matching neither neighbour |
| `ins_q_run` | quality of the inserted base | insertion rate per base, for bases lengthening a run |
| `del_fq_L1` … `del_fq_L5` | lower quality of the read bases either side | deletion rate per centre base in a run of 1, 2, 3, 4, 5 or more |
| `del_L` (optional) | run length 1-5 | deletion rate when qualities are unknown; default the median of `del_fq_L<n>` |
| `ins_r`, `del_r` | 0 | r of the geometric length model: a k-base indel costs r^(k-1) beyond its first base |

Keys need not be contiguous: a missing key takes the nearest key's rate. Rates
must be finite and non-negative. A non-zero exit fails `learn-errors` with the
command's stderr.

## Examples

| File | Language | Notes |
| --- | --- | --- |
| `indel_gam.R` | R + mgcv | Poisson GAMs with REML smoothing: `s(Q)` per insertion context, with the capped top quality (PacBio Q93) fitted apart; deletions in runs of 1 and 2 on a shared flanking-quality curve plus run-length deviations penalised towards it. Runs of 3 or more keep dada2-rs's native pooled rates: their events are too sparse on Illumina for the shared curve to fit them. |

```bash
dada2-rs learn-errors *.fastq.gz -o err.json \
    --indel-rates-out indel_rates.tsv \
    --indel-fit-cmd "Rscript examples/external_indel_fit/indel_gam.R"
dada2-rs dada sample.fastq.gz --error-model err.json --indel-rates indel_rates.tsv -o sample.json
```

`--indel-rates-out` records the command as `# indel_fit`.
