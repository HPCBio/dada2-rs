# Inside `b_compare`: the screen, the DP kernel, and where the time actually goes

*Issue [#127](https://github.com/HPCBio/dada2-rs/issues/127). Measurement only —
no algorithm or output changed. 362-sample MiSeq (R1 + R2) and 95-sample PacBio
HiFi, pooled `dada`, 24 threads, both platforms pinned to the same node.*

## Why this measurement

Closing [#124](shuffle-build-scan.md) left `b_compare` as the only phase big
enough to matter: **77% / 63% / 88% of `run_dada`** on MiSeq R1, MiSeq R2 and
PacBio respectively. Its parallel map runs at **95–99% efficiency on 24
threads**, so there is no threading headroom — any win has to be algorithmic.

We had comparison *counts* but no *time* split, and the counts alone pointed two
incompatible ways. The k-mer screen is paid on **100%** of comparisons; the
alignment only on the 12–25% that pass it. So a screen that is individually
cheap can still dominate by volume, and the ratio decides the target:

- **screen-dominated** → the lever is the screen itself: an inverted index or
  minimizers, replacing the ESPRIT-era linear scan (never the partition — see
  [k-mer screen provenance](kmer-size-screening.md)).
- **align-dominated** → the lever is the aligner, which is where most prior
  effort already went (#51 WFA, band size, #15 k-mer size).

`--verbose` now reports the split directly, gated on the verbose flag because an
`Instant` pair costs ~40–50 ns and the screen itself is sub-µs.

## The result

Percentages are of `b_compare`'s **busy** time (summed across worker threads);
`ns/comp` divides by the number of units each row is actually paid on.

| | **MiSeq R1** | **MiSeq R2** | **PacBio HiFi** |
|---|---|---|---|
| uniques (`nraw`) | 272,574 | 296,893 | 547,273 |
| clusters | 2,564 | 2,007 | 2,844 |
| comparisons | 482.6 M | 400.9 M | 771.7 M |
| passed the screen | 25.0% | 25.5% | **11.9%** |
| `b_compare` busy | 2,427 s | 1,427 s | 13,426 s |
| map / parallel eff. | 105 s / 97% | 62 s / 95% | 564 s / 99% |
| **k-mer screen** | 405.8 s (**16.7%**) | 327.7 s (**23.0%**) | 2,038.5 s (**15.2%**) |
| — per comparison | 841 ns | 818 ns | 2,642 ns |
| **align total** | 1,877.8 s (**77.4%**) | 994.3 s (**69.7%**) | 10,958.5 s (**81.6%**) |
| — DP kernel | 1,741.1 s (71.7%) | 909.3 s (63.7%) | 10,454.5 s (77.9%) |
| — `al2subs` + quals | 136.8 s (5.6%) | 85.0 s (6.0%) | 504.0 s (3.8%) |
| — per aligned pair | 15,537 ns | 9,713 ns | 119,634 ns |
| other (`compute_lambda`, overhead) | 143.7 s (5.9%) | 105.1 s (7.4%) | 429.1 s (3.2%) |

**Both platforms are align-dominated, and PacBio is the *most* align-dominated
of the three.** That refutes the prediction recorded when #127 was opened: k=7
gives PacBio 16,384-entry k-mer vectors against MiSeq's 1,024, and 88% of its
pairs are shrouded, so PacBio looked like the screen-bound case. The screen does
cost 3.1× more per comparison there (2,642 vs 841 ns) — but the alignment it
guards costs **7.7×** more (119.6 vs 15.5 µs), and it fires on 2.1× fewer pairs.
The two effects move in opposite directions and the alignment wins.

So there is one lever, not two. Rolling the DP kernel up to whole-run shares:

| | MiSeq R1 | MiSeq R2 | PacBio |
|---|---|---|---|
| DP kernel, % of `run_dada` | **43.7%** | **29.7%** | **62.8%** |
| DP kernel, % of pooled wall | ~41% | ~28% | ~58% |

**The banded Needleman-Wunsch DP kernel is the single hottest thing in pooled
`dada` — by a wide margin, on both platforms.** Every phase #124 examined is a
rounding error beside it.

## The screen is not the problem, and it is worth keeping

Worth stating plainly, because "16.7% of compare" invites the wrong conclusion.
Without the screen, every comparison would align: MiSeq R1's 482.6 M
comparisons × 15.5 µs ≈ **7,450 s** against the 2,427 s actually spent. On
PacBio it is 771.7 M × 119.6 µs ≈ **25.6 hours** of busy time against 3.7. The
screen returns roughly **3× on MiSeq and 7× on PacBio**.

A *perfect, free* index — one that costs nothing and shrouds exactly what the
current screen shrouds — is therefore capped at 15–23% of compare. Worth having
eventually; not the thing to build first, and not the thing #127 was opened to
find.

## Where the DP kernel's time goes

The per-aligned-pair costs are large enough to be worth decomposing. The
vectorized aligner walks anti-diagonals of a banded matrix, so cells computed
≈ `nrow × (band+1)` with `nrow = len1 + len2 + 1`:

| | cells/align | ns/align | **ns/cell** |
|---|---|---|---|
| MiSeq R1 (240 bp, band 16) | 8,177 | 14,405 | **1.76** |
| MiSeq R2 (160 bp, band 16) | 5,457 | 8,882 | **1.63** |
| PacBio (~1,455 bp, band 32) | 96,063 | 114,132 | **1.19** |

Two things are reassuring here: cost scales with cells (MiSeq R2/R1 = 1.62×
measured vs 1.50× from the length ratio), and the vectorized path is confirmed
in use — HiFi reads sit under the 3,500 bp `i16` guard and take
`homo_gap_p == gap_p`, so neither fallback fires.

But **1.2–1.8 ns/cell is 5–10× off what an 8-lane `i16` SIMD kernel should
achieve**, and the two platforms appear to miss for different reasons.

### It is not memory traffic — a falsified hypothesis

The first explanation tried was memory bandwidth, and it was wrong. Recording it
because the refutation is the useful part.

`align_vectorized_with_buf` allocated the **full** compressed matrix — `d16`
(`i16` scores) *and* `p8` (traceback pointers), both `ncol × nrow` — and
zero-filled both on every call. Each cell was therefore touched twice (zeroed,
then written) across 3 bytes: **6 bytes of write traffic per cell**, independent
of platform. That produced an apparent signature:

| | matrix (`d16` + `p8`) | traffic/align | GB/s/thread | × 24 threads |
|---|---|---|---|---|
| MiSeq R1 | 18 KB + 9 KB | 54 KB | 3.8 | 91 GB/s |
| MiSeq R2 | 12 KB + 6 KB | 36 KB | 4.1 | 99 GB/s |
| PacBio | 199 KB + 99 KB | 597 KB | 5.4 | **129 GB/s** |

All three land in the same 91–129 GB/s band despite a 12× spread in matrix size,
which reads like a bandwidth wall. **It is not.** The three figures agree because
traffic-per-cell is *fixed at 6 bytes by construction*, so dividing it by a
roughly constant ns/cell necessarily gives a roughly constant GB/s. The
convergence restates the ratio; it is not evidence about the bottleneck.

The structural observation was nonetheless sound: **`d16` is write-only outside a
3-row window.** The DP reads rows `row`, `row-1`, `row-2`; the traceback reads
`p8` alone and never touches `d16`. So the experiment was run — `d16` reduced to
three rotating one-row slots and its zero-fill dropped (199 KB → 210 B per
alignment on PacBio), holding `p8` full and byte-identity intact, verified by a
5,808-case fuzz against the retained pre-change implementation. If traffic were
binding, per-cell traffic falling 6 bytes → ~1 should have moved PacBio most.

Measured on the benchmark node, per-thread ns, full-matrix → rolling:

| threads | MiSeq R1 (240 bp, band 16) | PacBio (1455 bp, band 32) |
|---|---|---|
| 1 | 32,025 → 30,382 (1.05×) | 204,542 → 201,632 (1.01×) |
| 4 | 34,840 → 32,522 (1.07×) | 222,708 → 219,277 (1.02×) |
| 12 | 34,765 → 33,253 (1.05×) | 221,698 → 219,012 (1.01×) |
| 24 | 70,990 → 68,348 (1.04×) | 416,960 → 404,454 (1.03×) |

**Flat at every thread count, and PacBio — the 199 KB case, the one that does not
fit L2 — benefits *least*.** A bandwidth mechanism predicts the ratio rising with
contention and the large matrix gaining most; neither happened.

The full pooled A/B agrees, and supplies its own control. Because the change
touches only the DP, the **k-mer screen is an untouched channel in the same
runs** — anything it does is pure run-to-run drift:

| | prior → rolling | Δ | screen (untouched) | Δ |
|---|---|---|---|---|
| MiSeq R1 dp ns/comp | 14,377 → 14,654 | +1.9% | 843 → 898 | **+6.5%** |
| MiSeq R2 dp ns/comp | 8,903 → 9,092 | +2.1% | 757 → 797 | **+5.3%** |
| PacBio dp ns/comp | 113,834 → 110,667 | −2.8% | 2,659 → 2,700 | +1.5% |

`run_dada` moves +1.1% / +1.5% / −2.4% respectively. **The drift channel moves as
much as or more than the signal channel in every run**, so the noise floor on
this benchmark is ~2–6% — wider than the ~1.04× the microbenchmark predicted.
MiSeq's apparent regression and PacBio's apparent gain are both inside it. The
change is not measurable in production, in either direction.

Worth keeping as a technique: when a change touches one phase, the untimed
phases in the same run are a free null channel, and they calibrate the noise
floor far better than repeat runs would. Corroborating evidence
from the other direction: an Apple-silicon laptop, with an entirely different
memory system, reproduces the production ns/cell within 10% (1.75 vs 1.76 on
MiSeq R1) — which a DRAM-bound kernel would not do.

The kernel is **execution-bound**. `dploop` does vectorise (`sqadd.8h` is present
in the emitted assembly), so the remaining candidates are work *per cell* and
cells *per alignment*, not layout:

- The **diag-fill is a second full pass** over every cell of the anti-diagonal,
  walking `s1` in reverse — a pattern much less likely to vectorise than
  `dploop` itself.
- **Short anti-diagonals**: at band 16 the inner loop runs ~17 cells, two `i16`
  vectors plus a scalar remainder, so prologue/epilogue is amortised over almost
  nothing. Band 32 (~33 cells, four vectors) amortises better — consistent with
  PacBio's *lower* ns/cell despite its far larger matrix.

Byte-identity was confirmed on the production data, not just the fixtures: all
97 PacBio outputs identical, and 726/728 MiSeq, where the two exceptions are the
`seqtab` files whose row order is per-process `HashMap` nondeterminism — the ASV
sets are identical and the per-sample counts match under the column permutation
across all 362 samples in both reads.

The rolling-`d16` change is retained on the `perf/dp-rolling-d16-127` branch,
**unmerged**: byte-identical and harmless, but below the measurement floor.

### Threading saturates at 12, not 24

The thread ladder above carries a second result, visible in aggregate throughput
(threads ÷ per-thread ns) rather than per-thread latency:

| threads | 1 | 4 | 8 | 12 | 16 | 24 |
|---|---|---|---|---|---|---|
| MiSeq R1 | 0.031 | 0.115 | 0.230 | **0.345** | **0.347** | **0.338** |
| PacBio | 0.0049 | 0.0180 | 0.0362 | **0.0541** | **0.0544** | **0.0576** |

Scaling is near-perfect to 12 threads (per-thread cost 32.0 → 34.8 µs — note no
bandwidth degradation whatsoever), then **aggregate throughput stops dead** while
per-thread cost doubles by 24.

The node itself is not the constraint: it has 72 physical cores (144 SMT
threads); the newer node in the queue has 128 cores / 256 threads. The runs are
SLURM jobs asking for 24 threads, and the allocation turns out to be **12
physical cores plus their 12 SMT siblings** — confirmed from inside a job on
both nodes (`0-11,72-83` on the 72-core node; `0-11,128-139` on the 128-core
node that ran the production benchmarks):

```console
$ grep Cpus_allowed_list /proc/self/status
Cpus_allowed_list:      0-11,72-83
```

Linux enumerates each core's first SMT sibling as CPU 0–71 and its second as
72–143, so CPU *n* and CPU *n+72* are the same physical core. `0-11,72-83` is
precisely the pairs (0,72) … (11,83): 24 logical CPUs on 12 cores. Half the
threads share execution units with the other half, which is exactly the
doubling-at-constant-throughput the ladder shows. To check the sibling mapping
directly:

```console
$ lscpu -p=CPU,CORE | grep -v '^#' | awk -F, '$1==0 || $1==72'
0,0
72,0
```

Both CPUs report `CORE 0`, so the pairing is confirmed outright. The newer
128-core node follows the same convention with the offset at 128 (CPU 0 and
CPU 128 are both `CORE 0`).

Two cautions when repeating this. The check is only meaningful **inside a job**:
run from a plain shell on the node it reports the whole machine
(`Cpus_allowed_list: 0-255`) and says nothing about what an allocation would
grant. And where a job is *not* CPU-restricted, the Linux scheduler will
generally spread threads across distinct physical cores first — so the same
`--threads 24` can behave quite differently depending on whether the allocation
confines it to SMT-paired cores.

Both allocations have now been checked and both are 12 cores, so the ladder and
the production benchmarks share this topology. It is also worth carrying back to
earlier results —
absolute wall times here have been attributed to node generation
([#124's runs](shuffle-build-scan.md) were ~1.5–2× slower on an older node) — but
allocation topology is a competing explanation for part of that gap. Both are
untested, and the same in-job command settles them.

The cause is the job request rather than the code. `-n 24` asks for 24 *tasks*,
which the scheduler is free to satisfy with 12 SMT-paired cores. For a
single-process threaded program the shape to ask for is cores, and it works:

```console
$ srun -N 1 -n 1 -c 24 --threads-per-core=1 -p hpcbioamd --pty bash
$ grep Cpus_allowed_list /proc/self/status
Cpus_allowed_list:    0-23,128-151
```

`0-23,128-151` is the pairs (0,128) … (23,151) — **24 distinct physical cores**,
both siblings of each exposed — against 12 cores for the same nominal "24"
before. As a batch script:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=1      # or: --hint=nomultithread
...
dada2-rs dada-pooled --threads "${SLURM_CPUS_PER_TASK}" ...
```

`$SLURM_NTASKS` is the wrong variable to drive `--threads`: it counts tasks
(MPI ranks), not the cores one process may use. `$SLURM_CPUS_PER_TASK` is the
threading budget.

!!! tip "Check this from `srun`, not `salloc`"

    `salloc` grants the allocation but leaves the shell on the *login* node
    unless the cluster configures `SallocDefaultCommand` to step onto the
    compute node — so `Cpus_allowed_list` read after `salloc` reports the login
    node's cpuset and looks fine no matter what the job got. Read it from
    `srun --pty bash`, and confirm the hostname changed before trusting it.

**Expected effect: ~2× on the parallel portion of `run_dada`** — which is
63–88% of it — for a job-script change and no code change. On nodes of 72 and
128 cores there is a great deal more headroom above 24.

!!! warning "Every absolute constant on this page was measured on 12 cores"

    Confirmed on the node that ran the production benchmarks (`-N 1 -n 24`,
    128-core node):

    ```console
    $ grep Cpus_allowed_list /proc/self/status
    Cpus_allowed_list:    0-11,128-139
    $ echo $SLURM_NTASKS $SLURM_CPUS_ON_NODE $SLURM_TASKS_PER_NODE
    24 24 24
    ```

    CPU *n* and *n+128* are the same physical core, so that allocation is
    **cores 0–11: 12 physical cores presented as 24 logical CPUs**. Every
    per-unit cost quoted here is therefore roughly **2× what a dedicated core
    would show**, and the pooled runs were using about a tenth of a 128-core
    node.

    Note that all three SLURM variables report `24`, and none of them can
    distinguish cores from threads — which is how a recommendation to drive
    `--threads` from `$SLURM_NTASKS` produces this. Nor could
    `map parallel efficiency`: it is `busy ÷ (map × nthreads)`, so threads each
    running at half speed report as fully efficient. It measures whether the
    threads were *busy*, never whether they had cores to be busy on.

    The **ratios** on this page are unaffected — every run shared these
    conditions — and they carry all of its conclusions.

## Design constants

For costing any future change without re-running (24 threads, pinned node;
`ns/comp` and `ns/align` are per-thread busy time):

| constant | MiSeq (k=5, band 16) | PacBio (k=7, band 32) |
|---|---|---|
| k-mer screen | 818–841 ns/comp | 2,642 ns/comp |
| screen pass rate | 25.0–25.5% | 11.9% |
| DP kernel | 1.63–1.76 ns/cell | 1.19 ns/cell |
| `al2subs` + qual mapping | 830–1,132 ns/align | 5,502 ns/align |
| `compute_lambda` + overhead | 5.9–7.4% of compare | 3.2% of compare |
| comparisons per unique | 1,350–1,771 | 1,410 |

Measured at `--threads 24` on an allocation confirmed to be 12 physical cores
(see the warning above), so as absolutes these run ~2× a dedicated core; they
are reliable for modelling *relative* changes. DP write traffic is a fixed 6 bytes/cell on
both platforms and is **not** a useful cost term — see the falsification above.

Note the last row: `b_compare` is the `O(nraw × nclusters)` term, and at ~1,400
comparisons per unique on both platforms it is the reason pooled runs scale the
way they do. Nothing in this page changes that shape — it identifies which
constant in front of it is worth attacking.

## What this dictates

1. **`b_compare` is align-bound on both platforms.** One lever, not the two
   platform-specific ones #87 ended with.
2. **The k-mer screen is closed as a perf target.** It returns 3–7× on what it
   costs, and a perfect free replacement caps out at 15–23% of compare. Revisit
   only if the aligner gets fast enough to change the ratio.
3. **The DP kernel is the target** — 30–63% of `run_dada`, the largest single
   share this project has measured.
4. **The kernel is execution-bound, so memory layout is not the lever.** Tested
   and falsified: eliminating 5/6 of the DP's memory traffic bought ~1.04×, and
   *least* on the platform with the largest matrix. What remains is work per
   cell — the diag-fill's second pass over every cell, with its reverse `s1`
   walk — and cells per alignment, i.e. band width or a different algorithm.
5. **`al2subs` is not the problem** (3.8–6.0%), which retires the hypothesis
   that post-processing explained the per-pair cost.
6. **Check the SLURM allocation's topology before any further perf work.**
   Kernel throughput saturates at 12 of the 24 allocated threads on a 72-core
   node, which means every constant here is ~2× a dedicated core's and that
   `map parallel efficiency` cannot detect the ceiling. If the allocation is
   12 cores × 2 SMT siblings, fixing the request is worth more than anything
   else on this page — and it is a job-script change, not a code change.
