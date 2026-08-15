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

### Memory traffic: falsified at 24 threads, then confirmed at scale

This section recorded a falsified hypothesis for most of #127's life. It is
kept, because the sequence is the useful part — the hypothesis was wrong about
the *mechanism at 24 threads*, and the measurement that "refuted" it was wrong
about everything.

`align_vectorized_with_buf` allocated the **full** compressed matrix — `d16`
(`i16` scores) *and* `p8` (traceback pointers), both `ncol × nrow` — and
zero-filled both on every call. Each cell was touched twice across 3 bytes:
**6 bytes of write traffic per cell**, independent of platform. That produced an
apparent signature:

| | matrix (`d16` + `p8`) | traffic/align | GB/s/thread | × 24 threads |
|---|---|---|---|---|
| MiSeq R1 | 18 KB + 9 KB | 54 KB | 3.8 | 91 GB/s |
| MiSeq R2 | 12 KB + 6 KB | 36 KB | 4.1 | 99 GB/s |
| PacBio | 199 KB + 99 KB | 597 KB | 5.4 | **129 GB/s** |

**That convergence is an artifact and proves nothing.** Traffic per cell is
fixed at 6 bytes *by construction*, so dividing by a roughly constant ns/cell
necessarily gives a roughly constant GB/s. It restates the ratio.

The structural observation was sound, though: **`d16` is write-only outside a
3-row window.** The DP reads rows `row`, `row-1`, `row-2`; the traceback reads
`p8` alone and never touches `d16`. So the experiment was run — `d16` reduced to
three rotating one-row slots and its zero-fill dropped (199 KB → 210 B per
alignment on PacBio), `p8` left full, byte-identity verified by a 5,808-case
fuzz against the retained pre-change implementation.

**Measured unpinned, it looked flat**, and was twice recommended for closure:
1.04× in a microbenchmark, "within noise" in a 24-thread production A/B, and a
"null" from a 128-thread A/B/A/B whose two paired estimates disagreed by 6
percentage points.

**Measured with both arms pinned to one NUMA domain, it is unambiguous:**

| threads | rep1 | rep2 | screen (untouched control) |
|---|---|---|---|
| 24 | **−3.71%** | **−3.39%** | −1.1% / −0.9% |
| 64 | **−4.72%** | **−4.73%** | −1.1% / −0.5% |
| 128 (SMT) | **−5.80%** | **−6.18%** | −1.1% / −0.1% |

`run_dada` is −2.7 to −4.2%, the parallel map −3.0 to −5.1%, output byte-identical
(95/95 PacBio in every pair). Replicates agree to 0.4 points; `main` reproduces
itself to 0.04–0.14%. The effect **scales monotonically with thread count**,
which is what a memory-pressure mechanism should do and what no amount of drift
would produce.

So the honest statement of the mechanism is narrower than either the original
hypothesis or its refutation:

- At **24 threads** the DP kernel is execution-bound — memory traffic is not the
  binding constraint — and yet removing 5/6 of it still buys 3.5%. Not free,
  just not binding.
- By **64–128 threads** it is bandwidth-bound, and the same change buys 4.7–6%.

The independent confirmation is placement rather than thread count: at a *fixed*
64 threads, interleaving pages across both NUMA domains instead of binding them
to one moves the DP kernel **20%** (149,339 → 119,980 ns/comp). Same code, same
thread count, only page placement — that is a bandwidth effect and nothing else.

The measurement story is the more transferable result and has its own page:
[Measuring on a NUMA node](measuring-on-numa.md).

### A false alarm: SMT, and what the topology reading actually described

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

Two cautions when repeating this, and the second one caught us out.

The check is only meaningful **inside a job** — from a plain shell on the node
it reports the whole machine and says nothing about an allocation. And it must
be read in the *same context* as the run being explained: an interactive
`srun --pty bash` binds its shell to the allocation, whereas a batch script
invoking the binary directly typically runs **unconfined**, because a SLURM
allocation is only enforced by affinity or a cpuset cgroup when the site sets
`ConstrainCores` or the step is launched under `srun`.

That distinction is the whole story here. The 12-core readings above come from
interactive `srun` sessions. **The production benchmarks are batch jobs and were
never confined**, which a four-way A/B makes unambiguous (362-sample MiSeq, R1):

| job request | `--threads` | `run_dada` | map | DP ns/comp | screen ns/comp |
|---|---|---|---|---|---|
| `-n 24` | 24 | 165.0 s | 104.0 s | 14,377 | 843 |
| `-c 24 --threads-per-core=1` | 24 | 165.5 s | 104.3 s | 14,373 | 849 |
| `-n 48` | 48 | 120.7 s | 58.3 s | 14,601 | 1,201 |
| `-c 48 --threads-per-core=1` | 48 | **118.6 s** | **58.2 s** | 14,521 | 1,200 |

**`--ntasks` and `--cpus-per-task` are indistinguishable**, because with nothing
enforced both merely set a thread count. And 24 threads show no sign of SMT
contention: per-comparison DP cost is flat from 24 to 48 threads, and the map
scales 1.79×, neither of which is possible on 12 shared cores. Outputs are
bit-identical across thread counts (362/362 per-sample files).

So the constants elsewhere on this page stand as measured — an earlier version
of this page claimed they were ~2× inflated by a 12-core allocation, and that
claim is **withdrawn**. The lesson is narrower and more portable: *a topology
reading describes the context it was taken in*, and an interactive shell is not
the context your batch jobs run in.

What remains true is that the two phases scale differently, which is a result in
its own right:

- **DP kernel: flat** (14,373 → 14,521 ns/comp from 24 to 48 threads).
- **k-mer screen: +41%** (849 → 1,200 ns/comp).

The screen streams k-mer vectors and contends for bandwidth; the DP is
execution-bound and does not. This is a direct experimental separation of the
two, and it independently corroborates the falsification above. It also means
**the screen's share of `b_compare` grows with thread count**, so a screen-side
optimisation is worth more on wider machines.

This holds **up to 48 threads only** — past that the DP contends as well
(+34% by 128). See
[Thread scaling](#thread-scaling-and-the-ranking-it-overturns).

Practical upshot for these benchmarks: **48 threads is worth −20 to −28% of
`run_dada` over 24**, free and output-identical. There turns out to be little
room above that — see
[Thread scaling](#thread-scaling-and-the-ranking-it-overturns). Whether to run wider than the allocation is a policy question
for the site, not a code one — `--verbose` now reports the allocation, the
reachable core count, and flags when they disagree.

!!! tip "Check this from `srun`, not `salloc`"

    `salloc` grants the allocation but leaves the shell on the *login* node
    unless the cluster configures `SallocDefaultCommand` to step onto the
    compute node — so `Cpus_allowed_list` read after `salloc` reports the login
    node's cpuset and looks fine no matter what the job got. Read it from
    `srun --pty bash`, and confirm the hostname changed before trusting it.

!!! warning "`map parallel efficiency` cannot detect any of this"

    The statistic is `busy ÷ (map × nthreads)`, so threads each running at half
    speed on shared SMT siblings still report as fully efficient. It measures
    whether the threads were *busy*, never whether they had cores to be busy
    on. #127 was opened on the reading that 99% efficiency meant zero threading
    headroom and therefore that any win had to be algorithmic — which is not
    something that statistic can establish. It happens to have been true here,
    but only by luck.

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

Measured at `--threads 24` on an unconfined batch job (see above): the threads
had real cores, so these are not SMT-deflated. They are **thread-count
specific** as well as machine-specific — the DP is +34% and the screen +180% at
128 threads — so recalibrate at whatever thread count a design will run at. See
[Thread scaling](#thread-scaling-and-the-ranking-it-overturns) and
[How far these numbers travel](#how-far-these-numbers-travel). DP write traffic is a fixed 6 bytes/cell on
both platforms and is **not** a useful cost term — see the falsification above.

Note the last row: `b_compare` is the `O(nraw × nclusters)` term, and at ~1,400
comparisons per unique on both platforms it is the reason pooled runs scale the
way they do. Nothing in this page changes that shape — it identifies which
constant in front of it is worth attacking.

## Thread scaling, and the ranking it overturns

Four-point sweeps on **both** platforms (128-core node), run after the
measurements above because the 24-vs-48 comparison raised the question of where
this stops.

**362-sample MiSeq (R1):**

| threads | `run_dada` | map (parallel) | **serial remainder** | `shuffle` | `store` | DP ns/comp | screen ns/comp |
|---|---|---|---|---|---|---|---|
| 24 | 165.5 s | 104.3 s | **61.2 s** | 31.7 s | 13.0 s | 14,373 | 849 |
| 48 | 118.6 s | 58.1 s | **60.4 s** | 31.9 s | 12.1 s | 14,521 | 1,200 |
| 96 | 104.7 s | 42.3 s | **62.4 s** | 33.7 s | 11.9 s | 18,996 | 1,917 |
| 128 | 101.0 s | 37.2 s | **63.8 s** | 34.1 s | 12.2 s | 19,468 | 2,379 |

**95-sample PacBio HiFi:**

| threads | `run_dada` | map (parallel) | **serial remainder** | `shuffle` | `store` | DP ns/comp | screen ns/comp |
|---|---|---|---|---|---|---|---|
| 24 | 703.5 s | 577.4 s | **126.1 s** | 58.7 s | 21.7 s | 116,087 | 2,789 |
| 48 | 414.8 s | 297.4 s | **117.4 s** | 53.7 s | 20.2 s | 116,199 | 3,207 |
| 96 | 337.0 s | 207.8 s | **129.2 s** | 61.6 s | 18.9 s | 144,845 | 6,372 |
| 128 | 320.8 s | 185.0 s | **135.8 s** | 62.5 s | 21.9 s | 165,622 | 7,967 |

### The serial floor is real

The serial remainder is flat on both platforms across a 5.3× thread range —
**61–64 s on MiSeq, 117–136 s on PacBio** — while the parallel map falls by
2.8–3.1×. `store` barely moves; `shuffle` drifts up ~7% on MiSeq, consistent
with a wider `compare` map evicting more of the cache its scattered raw-major
access depends on, but far too small to matter.

So `run_dada` has an Amdahl floor near **62 s (MiSeq)** and **~127 s (PacBio)**,
against 101 s and 321 s at 128 threads. Everything remaining in thread count is
worth at most ~1.6× on MiSeq and ~2.5× on PacBio.

That the floor appears on **PacBio too** is the surprise. PacBio is much more
parallel — `compare` is 88% of `run_dada` at 24 threads against MiSeq's 77% —
and the expectation was that its crossover would sit above any reachable thread
count. It doesn't: being more parallel moves the crossover (the serial share at
128 threads is 42% on PacBio against 63% on MiSeq) but does not remove the
floor.

### Returns fall off a cliff after 48

| step | MiSeq gain (map eff.) | PacBio gain (map eff.) |
|---|---|---|
| 24 → 48 | **−28.4%** (90%) | **−41.0%** (97%) |
| 48 → 96 | −11.7% (62%) | −18.8% (69%) |
| 96 → 128 | −3.5% (53%) | −4.8% (59%) |

**48 threads is the sweet spot on both platforms**, and it is where these
benchmarks should run. PacBio rewards width more throughout — it is the more
parallel workload — so 96 remains defensible there (1.9× the core-seconds for
19% less wall) where the same step on MiSeq buys only 12%. Neither platform
justifies 128.

In core-seconds, MiSeq costs 3,972 at 24 threads, 5,693 at 48 and 12,928 at 128
— **3.3× the compute for 39% less wall**. Where core-hours are charged, 24
threads is the most efficient and 48 is the reasonable turnaround/cost
tradeoff.

### The DP kernel is only execution-bound up to ~48 threads

The claim earlier on this page that the DP is flat while the screen degrades
holds to 48 threads and **not beyond** — reproducibly, on both platforms:

| | MiSeq 24 → 48 | MiSeq 48 → 128 | PacBio 24 → 48 | PacBio 48 → 128 |
|---|---|---|---|---|
| DP ns/comp | +1.0% (flat) | **+34%** | +0.1% (flat) | **+43%** |
| screen ns/comp | +41% | +98% | +15% | +148% |

The mechanism is **memory bandwidth**, established independently of thread count
by varying only page placement: at a fixed 64 threads, binding to one NUMA
domain instead of interleaving across both costs the DP kernel **24%**
(119,776 → 149,747 ns/comp). Same code, same thread count, same data — only how
many memory controllers are serving it. Bound-64 also lands almost exactly on
the 96-thread default-placement value, i.e. **one domain at 64 threads
reproduces the pressure of two domains at 96**.

So the picture is: below ~48 threads the DP is execution-bound and the traffic
it moves is not the constraint; above that, bandwidth is. Both facts are needed
to read the
[rolling-`d16` result](#memory-traffic-falsified-at-24-threads-then-confirmed-at-scale),
which gains 3.5% in the first regime and 6% in the second.

See [Measuring on a NUMA node](measuring-on-numa.md) for the placement data and
for why none of this was visible until the benchmark rig was fixed.

### What this overturns

`b_shuffle` is serial, so its share of `run_dada` grows with every thread added:

| threads | 24 | 48 | 96 | 128 |
|---|---|---|---|---|
| `b_shuffle` share of `run_dada` | 19.2% | 26.9% | 32.2% | **33.8%** |

[#124 closed `b_shuffle`](shuffle-build-scan.md) partly on the grounds that no
remaining phase in it was worth more than ~3% of `run_dada`. That was measured
at 24 threads and **does not survive running wider**: at 48 threads the serial
block is already over half of `run_dada` and `shuffle` is two-thirds of it. The
move pass — the one phase with a favourable cost shape, sized at 1.8–2.9% and
deliberately not built — is worth roughly double that at 48 threads.

So the ranking has flipped. **At the thread counts these runs will actually use,
the serial phases outrank the DP kernel**, and #124's closure should be re-read
as conditional on thread count rather than final.

!!! success "Acted on, and it paid"

    [#132](https://github.com/HPCBio/dada2-rs/issues/132) took the one lever
    #124 had sized and declined — the move-pass dirty set — and built it against
    these numbers rather than the 24-thread ones. Measured at 64 threads: the
    move pass falls **59–66%**, and `run_dada` **4.3% / 5.6% / 1.8%**
    (MiSeq R1 / R2 / PacBio), byte-identical across 1,638 production files.

    That is a larger return than anything the DP-kernel work produced on MiSeq,
    from a phase this page had already established was 27% of `run_dada` at 48
    threads. The re-ranking was the finding; the code change was the
    consequence.

!!! warning "Do not tune for 48 threads at the expense of fewer"

    `dada2-rs` is also run on laptops and small allocations. Any change
    motivated by this table must be checked across the ladder — 4, 8, 24 as
    well as 48 — because a serial-phase optimisation that helps at 48 threads
    can trade against per-core work that dominates at 4. The thread count is a
    property of the deployment, not of the algorithm.

## How far these numbers travel

Everything above was measured on one machine class (AMD, 128 physical cores,
24–48 threads). `dada2-rs` also runs on laptops and other schedulers, so it is
worth being explicit about which of these results are properties of the
*algorithm*, which are properties of *this hardware*, and which depend on both.

**Invariant — determined by the data and parameters alone.** Portable
anywhere, and safe to design against:

- comparison counts (482.6 M / 400.9 M / 771.7 M) and comparisons per unique
  (~1,400 on both platforms)
- screen pass rates (25.0% / 25.5% / 11.9%) and shrouded fractions
- DP cells per alignment (`nrow × (band+1)`) and k-mer vector size (4^k bytes)
- the outputs themselves: bit-identical across thread counts (verified
  362/362 per-sample files at 24 vs 48 threads)

**Machine-dependent — rates.** Every `ns/comp`, `ns/cell` and GB/s figure
scales with clock, vector width and memory system. Use them for *relative*
modelling on comparable hardware, and recalibrate before trusting them
elsewhere — `bench_rolling_d16` (`--ignored`) re-derives the DP constant on any
machine in seconds. This repo has been burned once already by pricing a design
against laptop numbers ([#124](shuffle-build-scan.md)).

**Conditional — depends on machine balance *and* concurrency.** The important
category, because these are conclusions rather than measurements:

- *"`b_compare` is align-dominated."* The margin is wide (70–82% align vs
  15–23% screen) so the ordering is unlikely to flip, but the split is not
  fixed: the screen degrades faster than the DP with thread count (+148% vs
  +43% on PacBio, 24 → 128), so its share grows on wider machines and on any
  machine with less bandwidth per core. A screen-side optimisation is worth
  *more* at high thread counts.
- *"The DP kernel is bandwidth-bound."* True above ~48 threads on this hardware,
  and false below it — at 24 threads the kernel is execution-bound and the same
  memory traffic is not the constraint. Which regime a given machine sits in
  depends on its cores-to-bandwidth ratio, so the
  [rolling-`d16` change](#memory-traffic-falsified-at-24-threads-then-confirmed-at-scale)
  is worth 3.5% in the first regime and 6% in the second — a positive change
  either way here, but the *size* is not transferable.
- The best thread count. 48 beat 24 by 28% (MiSeq) and 41% (PacBio) here; the
  right number elsewhere depends on cores, bandwidth and what else shares the
  node.
- **Measurement conditions, not just hardware.** Every number on this page
  assumes fixed NUMA page placement. Under default placement the same rig
  produced a 25% noise floor that hid a 5% effect — see
  [Measuring on a NUMA node](measuring-on-numa.md). Re-deriving these constants
  elsewhere means controlling that too, not only matching the hardware.

**Site-specific.** All SLURM flags and confinement behaviour above describe one
cluster's configuration. The portable part is the *question* — how many
physical cores can this process actually use, and does that match what the
scheduler granted — which is why `--verbose` now reports both.

## What this dictates

1. **`b_compare` is align-bound on both platforms.** One lever, not the two
   platform-specific ones #87 ended with.
2. **The k-mer screen is closed as a perf target.** It returns 3–7× on what it
   costs, and a perfect free replacement caps out at 15–23% of compare. Revisit
   only if the aligner gets fast enough to change the ratio.
3. **The DP kernel is the target** — 30–63% of `run_dada`, the largest single
   share this project has measured.
4. **Memory layout *is* a lever above ~48 threads — take it.** Eliminating 5/6
   of the DP's write traffic is worth 3.5% at 24 threads and 4.7–6% at 64–128,
   byte-identical (the rolling-`d16` change). Below that the kernel is
   execution-bound and the remaining levers are work per cell — the diag-fill's
   second pass, with its reverse `s1` walk — and cells per alignment.
5. **`al2subs` is not the problem** (3.8–6.0%), which retires the hypothesis
   that post-processing explained the per-pair cost.
6. **Run at 48 threads.** Worth −28% of `run_dada` against 24,
   output-identical — more than any code change on this page returned. Returns
   collapse after that (−11.7% for 2× threads at 96, −3.5% at 128), and
   `run_dada` has an Amdahl floor near 62 s regardless. Confirm what the job
   actually gets from `--verbose` rather than the scheduler's variables, none
   of which distinguish cores from threads.
7. **The serial phases now outrank the DP kernel.** At 48 threads the serial
   remainder is over half of `run_dada` and `b_shuffle` is two-thirds of it
   (26.9% of `run_dada`, rising to 33.8% at 128). #124's closure was measured
   at 24 threads and should be re-read as conditional on thread count; its
   declined move-pass lever is worth roughly double there.
8. **Re-test the rolling-`d16` change at 128 threads.** It was falsified at 24,
   where the DP is genuinely execution-bound; by 128 the DP is +34% and
   contending. Cheap, and the branch already exists.
9. **Weight the next target by where it will run, and do not tune for 48
   threads alone.** The ranking here and every ns constant are properties of
   this hardware and this thread count — see
   [How far these numbers travel](#how-far-these-numbers-travel), and check any
   change across the ladder (4, 8, 24 as well as 48) before adopting it.
