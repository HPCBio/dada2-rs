//! Price `b_p_update`'s two candidate costs against each other (issue #154).
//!
//! ## The question
//!
//! `b_p_update` is 10-12% of `run_dada` on soil ITS2 and runs at ~56 ns per
//! repricing. #154 proposed parallelising it, on the reasoning that `get_pA`
//! calls `calc_pA`, which evaluates a regularised incomplete gamma — expensive
//! arithmetic that would scale with threads and, unlike the map, is not
//! competing for already-saturated bandwidth.
//!
//! Then the path counters went in, and on a small fixture **94% of repricings
//! never reach `calc_pA` at all** — they are singletons, which return `1.0`
//! after two comparisons. If that holds on real pools, the arithmetic cannot be
//! what the 56 ns is buying, and parallelising it would target a few percent of
//! the work.
//!
//! The alternative is the pattern this project keeps rediscovering. The loop
//! walks `Bi::raws`, a scattered `Vec<usize>` of raw indices, and reads four
//! fields — `reads`, `prior`, `comp.lambda`, `comp.hamming` — out of a
//! **160-byte `Raw`**. That is a full cache line pulled per raw to use ~20
//! bytes of it, in an order that defeats the prefetcher. It is exactly the
//! shape that made the serial store 83-87% strided read in #147, where the fix
//! was a dense parallel array and not parallelism.
//!
//! The two hypotheses want opposite work, so pricing them first is the whole
//! point:
//!
//! - **arithmetic-bound** -> parallelise `calc_pA` (#154 as filed)
//! - **memory-bound** -> hoist the hot fields into dense arrays (#147's fix)
//!
//! ## Arms
//!
//! | arm | what it does per item | isolates |
//! |---|---|---|
//! | `scattered` | reads 4 fields from a 160 B `Raw` via a shuffled index | today's traversal |
//! | `dense` | reads the same 4 fields from packed 16 B records | the layout fix |
//! | `calc_only` | `calc_pA` on resident values, no traversal | the arithmetic |
//! | `scattered_full` | scattered traversal + `calc_pA` at `--full-frac` | today, end to end |
//! | `dense_full` | dense traversal + `calc_pA` at `--full-frac` | the fix, end to end |
//!
//! `--full-frac` is the measured fraction of repricings that reach `calc_pA`;
//! take it from the `full calc_pA` line of a verbose run on **your** pool rather
//! than the default, which is a placeholder from a 896-raw fixture and is
//! certainly wrong for a real one.
//!
//! Every read goes through `black_box`. The first draft of `store_bandwidth`
//! reported 574 GB/s because single-array arms auto-vectorised and paired ones
//! could not, which measures SIMD width rather than access pattern.
//!
//! ## Usage
//!
//! ```text
//! cargo build --profile release-native --example pval_cost
//! numactl --interleave=all ./target/release-native/examples/pval_cost \
//!     --raws 825214 --rounds 50 --full-frac 0.05
//! ```
//!
//! ## It under-models, and by how much
//!
//! On a laptop this reports ~7 ns/repricing against production's ~56 ns. The
//! touched set here is `--visit x 160 B` (~16 MB at the ITS2 default), which
//! fits caches that the real traversal misses, and the loop omits work the real
//! one does: the scattered write-back of `p`, the bud-candidate bookkeeping, and
//! **two `Vec<usize>` clones of the member list per dirty cluster per round**
//! (`src/pval.rs:113` and `:169`) — allocation churn of the kind #147 found to
//! matter.
//!
//! So do not read absolute ns off this. What it can establish is an *ordering*
//! between the two candidate mechanisms, which is what decides the lever. Check
//! `scattered_full` against your own production ns/repricing (`p_update`
//! seconds / `raws repriced`, both on the `p-update churn` line); if the gap is
//! large, this benchmark is not modelling your bottleneck and the in-run path
//! counters are the instrument to trust instead.

use dada2_rs::pval::calc_pA;
use std::hint::black_box;
use std::time::Instant;

/// Same footprint as `Raw` (160 bytes), with the four fields `b_p_update`
/// actually reads placed as they are in the real struct.
#[derive(Clone, Copy)]
#[repr(C)]
struct RawLike {
    reads: u32,
    hamming: u32,
    lambda: f64,
    prior: bool,
    _pad: [u8; 143],
}

/// The same four fields, packed. What a dense parallel array would cost.
#[derive(Clone, Copy)]
#[repr(C)]
struct Packed {
    reads: u32,
    hamming: u32,
    lambda: f64,
}

fn arg<T: std::str::FromStr>(name: &str, default: T) -> T {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn main() {
    // Soil ITS2 R1: 825,214 uniques, 100,735 repriced per round.
    let nraw: usize = arg("--raws", 825_214);
    let visit: usize = arg("--visit", 100_735);
    let rounds: usize = arg("--rounds", 50);
    // Fraction of repricings reaching calc_pA. TAKE THIS FROM YOUR OWN RUN.
    let full_frac: f64 = arg("--full-frac", 0.05);

    println!(
        "raws={nraw} visited/round={visit} rounds={rounds} full_frac={full_frac}\n\
         Raw is {} B; packed record is {} B\n",
        std::mem::size_of::<RawLike>(),
        std::mem::size_of::<Packed>(),
    );

    let raws: Vec<RawLike> = (0..nraw)
        .map(|i| RawLike {
            reads: (i % 50 + 1) as u32,
            hamming: (i % 4) as u32,
            lambda: 1e-6 * ((i % 100) + 1) as f64,
            prior: false,
            _pad: [0; 143],
        })
        .collect();
    let dense: Vec<Packed> = raws
        .iter()
        .map(|r| Packed {
            reads: r.reads,
            hamming: r.hamming,
            lambda: r.lambda,
        })
        .collect();

    // Cluster membership is neither sequential nor random: raws are
    // abundance-sorted and members accumulate by cluster, so the traversal is
    // scattered but correlated. A deterministic stride with a large odd step
    // reproduces that better than either extreme, and keeps the arm reproducible.
    let step = 7919usize;
    let members: Vec<usize> = (0..visit).map(|i| (i * step) % nraw).collect();
    let bi_reads = 10_000u32;
    let period = if full_frac > 0.0 {
        (1.0 / full_frac).round() as usize
    } else {
        usize::MAX
    };

    let report = |name: &str, dur: std::time::Duration| {
        let per = dur.as_secs_f64() / (rounds * visit) as f64;
        println!(
            "{name:<18} {:>8.3}s total  {:>7.1} ns/repricing",
            dur.as_secs_f64(),
            per * 1e9
        );
    };

    // --- traversal only ---------------------------------------------------
    let t = Instant::now();
    for _ in 0..rounds {
        let mut acc = 0f64;
        for (n, &m) in members.iter().enumerate() {
            let r = &raws[m];
            acc += black_box(r.lambda)
                + black_box(r.reads) as f64
                + black_box(r.hamming) as f64
                + black_box(r.prior) as u8 as f64
                + n as f64;
        }
        black_box(acc);
    }
    report("scattered", t.elapsed());

    let t = Instant::now();
    for _ in 0..rounds {
        let mut acc = 0f64;
        for (n, &m) in members.iter().enumerate() {
            let r = &dense[m];
            acc += black_box(r.lambda)
                + black_box(r.reads) as f64
                + black_box(r.hamming) as f64
                + n as f64;
        }
        black_box(acc);
    }
    report("dense", t.elapsed());

    // --- arithmetic only --------------------------------------------------
    let t = Instant::now();
    for _ in 0..rounds {
        let mut acc = 0f64;
        for n in 0..visit {
            if n % period == 0 {
                acc += black_box(calc_pA(black_box(3), black_box(0.01), false));
            }
        }
        black_box(acc);
    }
    report("calc_only", t.elapsed());

    // --- end to end -------------------------------------------------------
    let t = Instant::now();
    for _ in 0..rounds {
        let mut acc = 0f64;
        for (n, &m) in members.iter().enumerate() {
            let r = &raws[m];
            let (reads, hamming, lambda) = (
                black_box(r.reads),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            acc += if n % period == 0 && hamming != 0 && lambda != 0.0 {
                calc_pA(reads.max(1), lambda * bi_reads as f64, false)
            } else {
                1.0
            };
        }
        black_box(acc);
    }
    report("scattered_full", t.elapsed());

    let t = Instant::now();
    for _ in 0..rounds {
        let mut acc = 0f64;
        for (n, &m) in members.iter().enumerate() {
            let r = &dense[m];
            let (reads, hamming, lambda) = (
                black_box(r.reads),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            acc += if n % period == 0 && hamming != 0 && lambda != 0.0 {
                calc_pA(reads.max(1), lambda * bi_reads as f64, false)
            } else {
                1.0
            };
        }
        black_box(acc);
    }
    report("dense_full", t.elapsed());

    println!(
        "\nscattered vs dense = the layout question (#147's fix applied here).\n\
         calc_only          = the arithmetic #154 proposed parallelising.\n\
         Whichever dominates scattered_full decides which lever is worth building."
    );
}
