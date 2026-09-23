# The pooled derep load: format-bound, not I/O-bound

**The serial load front of a pooled run was 87% JSON parsing.** The filesystem —
NFS, which had been the standing suspicion — ran at **868 MB/s**, indistinguishable
from local disk with a warm page cache. Fixing the parse took the phase from
**58.3 s to 21.7 s (−63%)**, its share of pooled wall from 14% to 6%, and total
pooled wall down **8.6%**, byte-identical.

Three things on this page are worth more than the speedup:

- **The measurement was built to choose between mutually exclusive fixes**, so
  the wrong one could not be built by accident.
- **The predicted cause was not the found cause.** The double-parse was expected
  to cost time; it cost 36% of peak RSS and *no measurable time*. The time was
  somewhere adjacent.
- **A lever was deliberately not built**, and the arithmetic for that decision is
  recorded so it does not get re-proposed.

## Why look at the load at all

[#41](https://github.com/HPCBio/dada2-rs/issues/41) made the pooled load
serial on purpose — streaming one sample at a time is what
[collapsed the memory peak](pooled-memory.md) — and left a comment justifying it:

> *the load (derep) is a tiny fraction of pooled wall, so serializing the load is
> cheap*

That was true when written. It stopped being true without anyone changing the
load, because the phase is flat in absolute terms while everything around it
shrinks:

| PacBio, 95 samples | 24 threads | 128 threads |
|---|---|---|
| `derep` | 50.7 s (**7%**) | 48.0 s (**12.5%**) |
| `run_dada` | 693.4 s (92%) | 320.8 s (83%) |

At 128 threads the whole serial block is 52% of pooled wall and `derep` is its
second-largest piece. This is the **third** time the same pattern has appeared in
this area — after [#124's closure](shuffle-build-scan.md) and the "screen is
closed" verdict in [#127](compare-screen-vs-align.md) — a correct decision made
against a share that has since moved. The rule it keeps teaching:
**re-measure the premise, not the design.**

## An instrument that picks the branch

The candidate fixes were mutually exclusive, which is what made the measurement
worth doing before any code:

- **~130 MB/s** → filesystem-bound. These benchmarks run on NFS rather than local
  scratch, so the answer is staging inputs locally and *no code change helps*.
- **~1 GB/s or more** → gzip and JSON. The answer is parallelism or a cheaper
  format.

A throughput number selects one branch and closes the other. That is a better
instrument than one which merely reports a total, and it is cheap: a
`--verbose` split of the load into read, parse and build.

Production, 95-sample PacBio on NFS:

```text
[dada-pooled] derep split (of 58.27s over 95 sample(s), 4951 MB uncompressed):
[dada-pooled]   read+gunzip     5.71s ( 9.8%)   868 MB/s
[dada-pooled]   parse          50.89s (87.3%)   (serde_json)
[dada-pooled]   build           0.56s ( 1.0%)
[dada-pooled]   per-sample   min 300ms  median 587ms  max 1196ms
```

**868 MB/s on the network filesystem**, against 838 MB/s measured locally with a
warm cache. Two things close at once: staging inputs to local scratch is
pointless, and the standing worry — raised when
[the NUMA work](measuring-on-numa.md) started — that NFS might be confounding
these benchmarks is retired *for this phase*.

## The fix, and the cause that was not predicted

`read_tagged_json` parsed every input **twice**: once into a `serde_json::Value`
to read the `dada2_rs_command` tag, then again via `from_value` into the target
type. The prediction was that the second parse was costing time.

Measured on a 339 MB local set, warm cache, so the format cost is isolated from
the filesystem question entirely:

| approach | parse | peak RSS after derep+merge |
|---|---|---|
| `Value` + `from_value` | 1.34 s | 1360 MB |
| tag-only, then direct parse | 1.36 s | **874 MB** |
| single pass | **0.90 s** | **871 MB** |

**Building a `Value` to read one string field cost 36% of peak RSS and no extra
time at all.** It allocates the whole document as a node tree; that is a memory
cost, not a CPU one. The *time* was the second scan itself — both passes read
every byte, so validating the tag separately genuinely doubles the work.

Both are fixed by carrying `dada2_rs_command` on the target struct and validating
it from the same parse that produces the data. The memory half benefits **every**
`read_tagged_json` caller, not just the pooled load.

Two effects that looked like one, separated only because the middle row of that
table was measured. Without it, the RSS win would have been silently credited to
the single-pass change.

## The fixture under-stated the effect

The local prediction was −33% on parse. Production delivered **−70%**.

The reason is in the throughputs: `main` parsed at 253 MB/s locally but only
**97 MB/s** at production scale, while the fixed path held up (377 → 323 MB/s).
The `Value`-tree approach degrades under sustained load in a way a six-sample
fixture never shows.

This is worth flagging because it runs against the grain. Small fixtures in this
repo have repeatedly **over**-stated effects — that caution appears on several
pages here. **This is the first case of one under-stating.** The asymmetry to
remember: a fixture that fits in cache and finishes in seconds cannot exhibit
allocator or memory-pressure behaviour, so it flatters whichever arm is the
greedier one.

## The lever deliberately not built

Bounded prefetch — loading *k* samples ahead on the thread pool while the fold
stays serial — was the obvious next step, and was sized against a phase that was
12.5% of pooled wall. After the parse fix that phase is 6%, and only 28% of *that*
is the read it would overlap. Optimistically it is worth **~1.5% of pooled wall**.

Against that: it has to preserve the sample-index fold order exactly, and bound
its own memory — the two invariants #41 was protecting when it serialized the
load in the first place. **Recommendation: do not build it** unless a workload
turns up with a much larger read column.

The remaining parse is still 71% of the reduced 21.7 s, so what is left belongs to
a cheaper format — [#1](https://github.com/HPCBio/dada2-rs/issues/1) — and **must
be costed against the new baseline, not the old one.** That is the same trap this
page opens with, one level down.

## What this dictates

- **Re-measure a premise before optimising against it.** Three phases in this
  project have now had a decision expire because a share moved while the absolute
  number stood still. A phase-share justification has a shelf life, and thread
  count is what expires it.
- **Prefer an instrument that selects between branches.** "How fast is the
  filesystem, really" answered a question that would otherwise have been settled
  by assumption — and NFS was the plausible assumption.
- **Cost models expire.** Prefetch was reasonable at 12.5% and is not at 6%. The
  format work inherits the same obligation.
- **`read_tagged_json` is a shared path.** Its memory behaviour affects every
  caller, so measure it there rather than in the one command that surfaced it.
- **A small fixture can under-state as well as over-state.** Where an arm's cost
  is allocation-heavy, only production scale will show it.
