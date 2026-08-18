//! Decompose what the parallel `b_compare` map loses that is not compute
//! (issue #147 follow-up).
//!
//! ## What this exists to settle
//!
//! Two questions came out of the #147 A/B runs, and the phase accounting is
//! not fine-grained enough to separate them.
//!
//! **1. Where does the +12% `user` time come from?** Removing
//! `b_compare_parallel`'s per-call allocation — either by shrinking the result
//! element under glibc's 32 MB `mmap` cap (lever A) or by reusing one buffer
//! (lever D) — collapses `sys` by 70-94% and *raises* `user` by 11-16%. Both
//! levers do it, at different element widths, so it tracks reusing the
//! allocation rather than its size.
//!
//! A rayon spin-wait explanation was **falsified by accounting** before this
//! benchmark was written: on soil ITS2 R2 the whole run has only 838
//! core-seconds of CPU outside the map's summed per-item busy time (8,056 total
//! vs 7,218 busy), while lever D's `user` rise alone is 1,211. There is no room
//! for it, and 48.7% overall occupancy says workers sleep through the serial
//! phases rather than burning CPU in them. So the rise is *inside* busy: the
//! workers' own stores got more expensive.
//!
//! The mechanism this tests: a freshly-`mmap`ped page is zeroed by the kernel
//! immediately before a worker writes it, so the line is already warm in the
//! local cache and the store is nearly free. A reused buffer's lines are cold
//! and dirty, so every store pays a read-for-ownership from DRAM plus a
//! writeback of what it evicts. Identical work, moved out of kernel zeroing
//! (`sys`) and into user-mode memory stalls (`user`) — and *more* of it, since
//! the kernel zeroes with non-temporal stores and the collect does not.
//!
//! If that is right, `collect_fresh` and `collect_reuse` differ in the same
//! direction and rough proportion as the production `user` figures, and
//! `collect_reuse_nt`-style mitigation becomes the question rather than which
//! lever to merge.
//!
//! **2. What does the map lose to load imbalance?** Production reports 86-90%
//! map parallel efficiency, which on soil ITS2 R2 is 1,166 core-seconds lost
//! *inside* the parallel region across 4,017 calls — 4.6 ms per thread per
//! call. That is far too long to be rayon's spin-then-sleep (microseconds), so
//! it is threads idling at the tail of each collect waiting on stragglers.
//! `join_uniform` and `join_skewed` do the *same total work* per call and differ
//! only in how it is distributed, so the gap between them is the imbalance cost
//! rather than the framework cost. Pass `--base 0` for the fork/join/wake floor
//! with no work at all.
//!
//! This is the largest unexamined parallel loss in the project and it is
//! independent of the allocation question, so both are measured here.
//!
//! ## Reading the results
//!
//! Absolute numbers are not the point; the *arms are matched pairs* and only
//! the differences carry the argument:
//!
//! - `collect_fresh` vs `collect_reuse` — the allocation question. Same work,
//!   same element, differing only in whether the destination is newly mapped.
//! - `join_uniform` vs `join_skewed` — the imbalance question. Same call count,
//!   thread count, **and total work**, differing only in how that work is
//!   distributed across items. An arm that simply does less work would measure
//!   work rather than imbalance.
//!
//! Per-call figures matter more than totals here, because production's cost
//! scales with call count (4,017 on ITS2, 11,283 on 16S), not with wall time.
//!
//! Run under `numactl --interleave=all` on an otherwise idle node, at the
//! thread count production uses. Page placement re-rolls per run and the
//! serial scattered phases have swung 25% between replicates of one binary
//! without it — see docs/findings/measuring-on-numa.md.
//!
//! ## Variants
//!
//! | arm | what it does per round | isolates |
//! |---|---|---|
//! | `collect_fresh` | parallel collect into a newly allocated `Vec` | today's `b_compare_parallel` |
//! | `collect_reuse` | parallel collect into one reused `Vec` | lever D |
//! | `collect_fresh16` | as `collect_fresh`, 16 B element | lever A, allocation kept |
//! | `collect_reuse16` | as `collect_reuse`, 16 B element | levers A+D |
//! | `join_uniform` | parallel map, cost spread evenly over items | work-matched control |
//! | `join_skewed` | same total work, concentrated 1 item in 32 | + load imbalance |
//!
//! ## Usage
//!
//! ```text
//! cargo build --profile release-native --example parallel_overhead
//! numactl --interleave=all ./target/release-native/examples/parallel_overhead \
//!     --raws 939532 --rounds 200 --threads 64
//! ```
//!
//! `--raws` defaults to the soil ITS2 R2 pool (939,532) so the 48 B result
//! vector is 45.1 MB — above glibc's 32 MB cap, which is the regime under
//! test. Drop it below ~700,000 and the fresh/reuse arms should converge,
//! because the allocation stops being `mmap`ped; that is a useful positive
//! control on the mechanism.

use rayon::prelude::*;
use std::hint::black_box;
use std::time::Instant;

/// Mirrors `b_compare_parallel`'s result element: `(f64, u32, bool, CompCost)`.
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct Item48 {
    lambda: f64,
    hamming: u32,
    skipped: bool,
    cost: [u64; 4],
}

/// The same without `CompCost` (lever A's production element).
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct Item16 {
    lambda: f64,
    hamming: u32,
    skipped: bool,
}

fn arg<T: std::str::FromStr>(name: &str, default: T) -> T {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Stand-in for per-item alignment cost. Production's raws are abundance-sorted
/// and a few percent trigger full alignment while the rest are k-mer screened
/// out, so per-task cost is heavily skewed — which is what `with_max_len`
/// exists to rebalance. `skew = 0` gives uniform cost for the `join_empty`
/// floor.
#[inline]
fn work(index: usize, base: usize, skew: usize) -> f64 {
    // `base` iterations for every item; the skewed arm multiplies that by
    // `skew` for one item in 32, standing in for the few percent of
    // comparisons that clear the k-mer screen and pay a full alignment.
    let iters = if skew > 1 && index.is_multiple_of(32) {
        base * skew
    } else {
        base
    };
    let mut acc = index as f64;
    for _ in 0..iters {
        acc = black_box(acc * 1.000_001 + 1.0);
    }
    acc
}

fn main() {
    // Soil ITS2 R2: 939,532 uniques, so the 48 B vector is 45.1 MB.
    let nraw: usize = arg("--raws", 939_532);
    let rounds: usize = arg("--rounds", 200);
    let threads: usize = arg("--threads", 0);
    // Production's alignment pass rate is 0.8% (ITS2) to 4.5% (16S); the skew
    // multiplier stands in for how much more a passing comparison costs.
    let skew: usize = arg("--skew", 64);
    // Per-item work, in iterations of a dependent FLOP. The default is set so
    // `join_skewed` lands near production's map cost on soil ITS2 R2 —
    // 131 s over 4,017 calls at 64 threads is 2,086 core-ms/call, or ~2.2 us
    // per raw. An under-costed stand-in makes the collect look dominant and
    // understates imbalance, which is the thing being measured. Pass `--base 0`
    // to strip the compute out entirely and read the fork/join/wake floor.
    let base: usize = arg("--base", 160);

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("failed to set rayon thread count");
    }
    let nthreads = rayon::current_num_threads();

    let bytes48 = nraw * std::mem::size_of::<Item48>();
    let bytes16 = nraw * std::mem::size_of::<Item16>();
    println!(
        "raws={nraw} rounds={rounds} threads={nthreads} base={base} skew={skew}\n\
         result vector: {:.1} MB at 48 B, {:.1} MB at 16 B \
         (glibc mmap cap is 32 MB)\n\
         production reference: soil ITS2 R2 map = 2086 core-ms/call \
         (131 s / 4017 calls x 64 threads)\n",
        bytes48 as f64 / 1e6,
        bytes16 as f64 / 1e6,
    );

    let grain = 32;
    let report = |name: &str, dur: std::time::Duration, bytes: usize| {
        let per_call = dur.as_secs_f64() / rounds as f64;
        let core_s = per_call * nthreads as f64;
        println!(
            "{name:<16} {:>8.3}s total  {:>8.3} ms/call  {:>8.1} ns/raw  \
             {:>7.1} GB/s  {:>8.1} core-ms/call",
            dur.as_secs_f64(),
            per_call * 1e3,
            per_call * 1e9 / nraw as f64,
            if bytes > 0 {
                bytes as f64 / per_call / 1e9
            } else {
                0.0
            },
            core_s * 1e3,
        );
    };

    // --- allocation: fresh mapping vs reused buffer -----------------------
    //
    // Identical work and identical element; the only difference is whether the
    // destination pages were just handed over zeroed by the kernel.

    let t = Instant::now();
    for _ in 0..rounds {
        let v: Vec<Item48> = (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| Item48 {
                lambda: work(i, base, skew),
                hamming: i as u32,
                skipped: false,
                cost: [0; 4],
            })
            .collect();
        black_box(&v);
    }
    report("collect_fresh", t.elapsed(), bytes48);

    let mut buf48: Vec<Item48> = Vec::with_capacity(nraw);
    // Fault the pages in once so the first round is not charged for it.
    buf48.resize(nraw, Item48::default());
    let t = Instant::now();
    for _ in 0..rounds {
        (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| Item48 {
                lambda: work(i, base, skew),
                hamming: i as u32,
                skipped: false,
                cost: [0; 4],
            })
            .collect_into_vec(&mut buf48);
        black_box(&buf48);
    }
    report("collect_reuse", t.elapsed(), bytes48);

    let t = Instant::now();
    for _ in 0..rounds {
        let v: Vec<Item16> = (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| Item16 {
                lambda: work(i, base, skew),
                hamming: i as u32,
                skipped: false,
            })
            .collect();
        black_box(&v);
    }
    report("collect_fresh16", t.elapsed(), bytes16);

    let mut buf16: Vec<Item16> = Vec::with_capacity(nraw);
    buf16.resize(nraw, Item16::default());
    let t = Instant::now();
    for _ in 0..rounds {
        (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| Item16 {
                lambda: work(i, base, skew),
                hamming: i as u32,
                skipped: false,
            })
            .collect_into_vec(&mut buf16);
        black_box(&buf16);
    }
    report("collect_reuse16", t.elapsed(), bytes16);

    // --- fork/join floor vs load imbalance --------------------------------
    //
    // No collect at all: a sum reduction, so the only memory traffic is the
    // index range. `join_empty` gives uniform per-item cost, `join_skewed`
    // gives production's. The gap is imbalance, not framework overhead.

    // Work-matched against `join_skewed`: same mean iterations per item, spread
    // uniformly instead of concentrated in one item per 32. Matching the *work*
    // is the whole point — an arm that simply does less is measuring work, not
    // imbalance, and the difference between these two is then the cost of the
    // distribution alone.
    let uniform = base + (base * (skew - 1)) / 32;
    let t = Instant::now();
    for _ in 0..rounds {
        let s: f64 = (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| work(i, uniform, 1))
            .sum();
        black_box(s);
    }
    report("join_uniform", t.elapsed(), 0);

    let t = Instant::now();
    for _ in 0..rounds {
        let s: f64 = (0..nraw)
            .into_par_iter()
            .with_max_len(grain)
            .map(|i| work(i, base, skew))
            .sum();
        black_box(s);
    }
    report("join_skewed", t.elapsed(), 0);

    println!(
        "\nfresh vs reuse is the allocation question (#147 levers A and D);\n\
         uniform vs skewed is the load-imbalance question behind the 86-90%\n\
         map parallel efficiency production reports."
    );
}
