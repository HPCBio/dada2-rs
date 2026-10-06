# Rename-pair tracing

For a member-order floor run (#157), how far is each **renamed** ASV's losing
sibling from budding on its own? That distribution is what issue #246 asks
for: whether `OMEGA_A = 1e-40` is close to separating the pairs it currently
merges.

A rename pair is one organism reported under two names in two arms: in one
sample, an ASV only the baseline reports and an ASV only the other arm
reports, within a few edits of each other.

| Script | What it does |
|---|---|
| `find_pairs.py` | lists the rename pairs, per sample, between a baseline arm and the other arms |
| `trace_pairs.sh` | re-runs each sample under the arms its pairs need, with `--cluster-trace` |
| `summarize.py` | per pair: is it reproduced, are both in one cluster, and how many orders of magnitude the sibling's pA × nraw sits above `OMEGA_A` |
| `trace_asv_pair.sh` | the same trace for one pair given by sequence, printed for reading |

## Running it

The floor run must be per-sample `dada`. Pooled and pseudo clusters span
samples, so a one-sample trace would not reproduce them.

```bash
# 1. the floor (dev/run_member_order.sh), if not already run
MODE=per-sample dev/run_member_order.sh "$BIN" errF.json floor filt/*.fastq.gz
# 2. the pairs
python3 dev/rename_pairs/find_pairs.py floor/insertion floor/shuffle_{1..5} -o pairs.tsv
# 3. the traces: same binary, model and EXTRA flags as the floor run
EXTRA="..." dev/rename_pairs/trace_pairs.sh "$BIN" errF.json filt pairs.tsv trace
# 4. the table (-o) and the distribution (stderr)
python3 dev/rename_pairs/summarize.py pairs.tsv trace -o summary.tsv
```

Step 3 is the only expensive one: one `dada` run per sample for the baseline,
plus one per other arm that a pair in that sample needs.

## Reading the numbers

**pA is the budding p-value, recomputed.** It is `calc_pA` conditioned on the
sequence being present, P(X ≥ reads | E) / (1 − e^−E), which is what the
budding test uses. Every `Abundance` birth's `birth_pval` in a trace equals
this × `nraw`.

The trace's per-member `pval` is something else: the final `OMEGA_C` pass's
**unconditioned** value (`prior = true`, as in R). For a doubleton that is
smaller by about a factor of E, so reading it as pA makes a sibling look much
closer to budding than it is. An earlier version of `trace_asv_pair.sh` did
that. The `pval` field is now used only as a check: `summarize.py` recomputes
it from reads and `e_reads` and warns on any disagreement.

**A pair that is not reproduced is excluded.** It was traced under a different
binary, model, FASTQ or flags than the floor run, so its numbers describe some
other run.

**Assumes no `--prior-file` and no singleton detection.** A raw with a prior
buds on the unconditioned value.
