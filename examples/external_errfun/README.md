# External error functions

`dada2-rs` lets you plug in a user-supplied error-fitting function via
`--errfun external --errfun-cmd "<command>"`. Use this when:

- a new sequencing technology needs a custom model and you don't want to
  rebuild dada2-rs,
- you want to compare dada2-rs's native fits against R DADA2's
  `loessErrfun` (or any other R function from the literature), or
- you're prototyping an alternative model in R / Python before
  considering whether to port it to Rust.

## Wire format

`dada2-rs` runs `<command> <trans-tsv-path> <err-tsv-path>` once per
self-consistency iteration. The command string is whitespace-split into
argv (no shell interpolation — wrap in a shell script if you need
quoting). Both files use a layout compatible with R's
`read.table(..., row.names = 1, header = TRUE, check.names = FALSE)`:

- The header is a leading tab followed by `nq` integer column labels
  (`0`, `1`, …, `nq-1`).
- Each of the 16 data rows starts with a row label
  (`A2A`, `A2C`, …, `T2T`) followed by `nq` whitespace-separated values.

Your script must:
1. Read the trans matrix (16 × `nq` non-negative integers) from `argv[1]`.
2. Compute a 16 × `nq` matrix of error rates (each value in `[0, 1]`).
3. Write it to `argv[2]` in the same shape. R's
   `write.table(err, args[2], sep = "\t", quote = FALSE, col.names = NA)`
   produces the expected format directly.

`dada2-rs` validates dimensions and value range, then enforces the
diagonal self-transition probabilities post-hoc, so you can either
reconstruct them yourself (as the R examples do) or ignore the diagonal
rows entirely.

## Provenance

The exact `--errfun-cmd` string is recorded in the output JSON's `params`
block as `errfun_cmd`, alongside `errfun: "external"`. This is purely
informational — the script contents are not hashed.

## Examples

| File | Language | Notes |
| --- | --- | --- |
| `loess_reference.R` | base R + stats | Verbatim port of `dada2:::loessErrfun`. Useful as a parity check against the native Rust `loess_errfun`. |
| `loess_modified.R` | base R + dplyr + magrittr | Adds the Salazar/Ruscheweyh `span = 2` + `log10(tot)` weights and the Holland-Moritz Q40-floor monotonicity fix from issues #791, #938, #1307. |
| `noqual.py` | Python 3, stdlib only | Pure-Python `noqualErrfun` equivalent. Demonstrates the contract for users without R. |

Invoke them via:

```bash
dada2-rs learn-errors *.fastq.gz \
    --errfun external \
    --errfun-cmd "Rscript examples/external_errfun/loess_reference.R" \
    -o err_R_reference.json
```
