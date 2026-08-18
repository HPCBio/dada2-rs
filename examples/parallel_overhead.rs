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
//! ## Running it on a node you do not have to yourself
//!
//! Prefer an idle node. When one is not available this benchmark is built to
//! survive a neighbour, by two deliberate choices:
//!
//! - **Arms run round-robin, one round each**, never as sequential blocks.
//!   Contention drifts over the life of a run, so blocked arms sample different
//!   node states and the between-arm difference picks up whatever the
//!   neighbouring job was doing during each block. Interleaved, every arm sees
//!   the same drift, and because only the matched pairs carry the argument the
//!   common component cancels.
//! - **The estimator is the per-round minimum**, not the mean. The
//!   least-disturbed round is the one closest to the uncontended truth; a mean
//!   is pulled around by the neighbour. Median and max print alongside it, and
//!   the `med/min` line is the contention readout — near 1.00 is a quiet node,
//!   and a wide spread is the signal to distrust the run rather than to reason
//!   about small differences in it.
//!
//! Neither trick rescues a badly oversubscribed node. If the neighbour is using
//! cores this run also wants, `join_skewed` is measuring the neighbour's
//! stragglers as well as its own and the imbalance number is not meaningful.
//! Leave headroom with `--threads` rather than oversubscribing, and say what
//! the node was doing when reporting the result.
//!
//! Run under `numactl --interleave=all`, at the thread count production uses.
//! Page placement re-rolls per run and the serial scattered phases have swung
//! 25% between replicates of one binary without it — see
//! docs/findings/measuring-on-numa.md.
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
    let uniform = base + (base * (skew - 1)) / 32;

    // Reused destinations. Faulted in once so the first round is not charged
    // for page faults the later rounds do not pay.
    let mut buf48: Vec<Item48> = vec![Item48::default(); nraw];
    let mut buf16: Vec<Item16> = vec![Item16::default(); nraw];

    // Arms are run **round-robin, one round each**, not as sequential blocks.
    //
    // On a shared node this is the difference between a usable measurement and
    // a worthless one. Contention drifts over the life of the run, so blocked
    // arms sample different node states and the between-arm difference picks up
    // whatever the neighbouring job happened to be doing during each block.
    // Interleaved, every arm sees the same drift, and since only the matched
    // pairs carry the argument the common component cancels.
    //
    // Per-round times are kept rather than summed so the reporting can use the
    // **minimum** as the estimator: the least-disturbed round is the one
    // closest to the uncontended truth, while the mean is pulled around by
    // whatever else is on the node. Median and max are printed alongside so
    // contention is visible instead of silent -- a wide min-to-median spread is
    // the signal to distrust the run.
    type Arm<'a> = (&'a str, Box<dyn FnMut() + 'a>, usize);
    let mut arms: Vec<Arm> = vec![
        (
            "collect_fresh",
            Box::new(|| {
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
            }),
            bytes48,
        ),
        (
            "collect_reuse",
            Box::new(|| {
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
            }),
            bytes48,
        ),
        (
            "collect_fresh16",
            Box::new(|| {
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
            }),
            bytes16,
        ),
        (
            "collect_reuse16",
            Box::new(|| {
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
            }),
            bytes16,
        ),
        (
            "join_uniform",
            Box::new(|| {
                let s: f64 = (0..nraw)
                    .into_par_iter()
                    .with_max_len(grain)
                    .map(|i| work(i, uniform, 1))
                    .sum();
                black_box(s);
            }),
            0,
        ),
        (
            "join_skewed",
            Box::new(|| {
                let s: f64 = (0..nraw)
                    .into_par_iter()
                    .with_max_len(grain)
                    .map(|i| work(i, base, skew))
                    .sum();
                black_box(s);
            }),
            0,
        ),
    ];

    let mut times: Vec<Vec<f64>> = vec![Vec::with_capacity(rounds); arms.len()];
    // One untimed round per arm: first touch, allocator warm-up, and branch
    // predictor state should not land on the first measured round.
    for (_, run, _) in arms.iter_mut() {
        run();
    }
    for _ in 0..rounds {
        for (a, (_, run, _)) in arms.iter_mut().enumerate() {
            let t = Instant::now();
            run();
            times[a].push(t.elapsed().as_secs_f64());
        }
    }

    println!(
        "{:<16} {:>10} {:>10} {:>10} {:>10} {:>9} {:>12}",
        "arm", "min ms", "med ms", "max ms", "ns/raw", "GB/s", "core-ms"
    );
    let mut mins = Vec::with_capacity(arms.len());
    for (a, (name, _, bytes)) in arms.iter().enumerate() {
        let mut v = times[a].clone();
        v.sort_by(|x, y| x.partial_cmp(y).expect("no NaN in timings"));
        let (min, med, max) = (v[0], v[v.len() / 2], v[v.len() - 1]);
        mins.push(min);
        println!(
            "{name:<16} {:>10.3} {:>10.3} {:>10.3} {:>10.1} {:>9.1} {:>12.1}",
            min * 1e3,
            med * 1e3,
            max * 1e3,
            min * 1e9 / nraw as f64,
            if *bytes > 0 {
                *bytes as f64 / min / 1e9
            } else {
                0.0
            },
            min * nthreads as f64 * 1e3,
        );
    }

    // The two comparisons this benchmark exists for, on the min estimator.
    let pct = |a: usize, b: usize| (mins[b] - mins[a]) / mins[a] * 100.0;
    println!(
        "\nallocation  48 B: reuse is {:+.1}% vs fresh   16 B: {:+.1}%",
        pct(0, 1),
        pct(2, 3),
    );
    println!(
        "imbalance        : skewed is {:+.1}% vs work-matched uniform",
        pct(4, 5),
    );
    // Contention check: on an idle node med/min sits near 1.00. Well above that
    // and the arms were not sampling the same machine.
    let spread: Vec<String> = arms
        .iter()
        .enumerate()
        .map(|(a, (name, _, _))| {
            let mut v = times[a].clone();
            v.sort_by(|x, y| x.partial_cmp(y).expect("no NaN in timings"));
            format!("{name} {:.2}", v[v.len() / 2] / v[0])
        })
        .collect();
    println!("med/min (1.00 = quiet node): {}", spread.join("  "));

    println!(
        "\nfresh vs reuse is the allocation question (#147 levers A and D);\n\
         uniform vs skewed is the load-imbalance question behind the 86-90%\n\
         map parallel efficiency production reports."
    );
}
