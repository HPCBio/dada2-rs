//! Price the layout candidates for `b_p_update`'s hot loop (issue #154).
//!
//! ## What this decides
//!
//! #154's premise — that `calc_pA` is the cost — was falsified: it is reached by
//! 4.6-7.9% of repricings across four arms, while the phase runs at a nearly
//! flat 56.6-65.4 ns per repricing across pools differing 1.5x in `nraw`. A
//! per-unit cost that does not move with the working set is a cache-line
//! signature, so the target is the scattered gather, not the arithmetic.
//!
//! What is *not* settled is which layout change to make, and the obvious one
//! turns out not to work:
//!
//! | candidate | why it might win | why it might not |
//! |---|---|---|
//! | reorder `Raw`'s fields | free, exact, no duplication | the hot set is **already** contiguous (offsets 104..146); at 152 B and align 8, that 42-byte span straddles a 64-byte line ~2/3 of the time wherever it is put |
//! | `repr(C, align(64))` on `Raw` | guarantees one line per repricing | pads `Raw` 152 -> 192 B, **+26% on the `Raw` array** (125 MB -> 158 MB on soil ITS2) |
//! | dense parallel arrays | ~32 B/raw touched instead of 152, and the hot arrays may become L3-resident | *three* arrays gathered instead of one struct, so potentially three misses instead of two |
//!
//! The third is what #147 did for `e_minmax`, and it is the one whose payoff is
//! hardest to predict: it trades fewer bytes per raw against more independent
//! streams. On the 7713 the deciding question is whether ~26 MB of hot arrays
//! stays resident in a 32 MB per-CCD L3 while a 125 MB `Raw` array cannot.
//!
//! ## Arms
//!
//! | arm | models |
//! |---|---|
//! | `raw_current` | today: 152 B, align 8, hot fields at offset 104 |
//! | `raw_aligned` | `repr(C, align(64))`, hot fields first, 192 B |
//! | `dense_split` | `p` and `comp` in dense arrays; `reads`/`prior` still in `Raw` |
//! | `dense_all` | one packed 32 B record — the ceiling, if every hot field could move |
//! | `prefetch` | today's layout, with a software prefetch issued `--distance` iterations ahead |
//!
//! `dense_all` is deliberately **not implementable as written**: `reads` and
//! `prior` are read all over the codebase and moving them is a much larger
//! refactor than this phase justifies. It is here as an upper bound — if
//! `dense_split` lands close to it, the cheap version captures the win; if
//! `dense_all` is far better, the refactor has a price worth quoting.
//!
//! ## The size prediction, which is the thing to test
//!
//! With a working instrument the win tracks **whether the hot array fits L3**,
//! not how many bytes each access reads. A serial phase runs on one core, so it
//! sees one CCD's 32 MB, not the 7713's 256 MB total. That makes the payoff
//! pool-dependent:
//!
//! | pool | `nraw` | packed (32 B/raw) | fits 32 MB L3? |
//! |---|---|---|---|
//! | soil ITS2 | 825,214 | 26.4 MB | **yes** |
//! | soil 16S | 1,225,523 | 39.2 MB | **no** |
//!
//! So `dense_all` should win on ITS2 and win **less, or not at all, on 16S** —
//! and if it wins equally on both, the residency explanation is wrong and the
//! benefit is simply fewer bytes touched. Run both before building anything.
//! `--raws` also sweeps this directly: shrink it until even `raw_current` fits
//! and every arm should converge.
//!
//! ## Why `prefetch` is the interesting arm
//!
//! Every layout arm changes *where the data lives*. The prefetch arm changes
//! only *when it is asked for*: the member list is fully known before the loop
//! runs, so the address needed `D` iterations from now is available now. That
//! makes it the one candidate which
//!
//! - needs no duplication and no second source of truth,
//! - is **byte-identical by construction** — the loop body is untouched, only
//!   the memory system is warned, and
//! - does not depend on the hot array fitting L3, so it should not decay as
//!   pools grow. `dense_all` wins 2.1x on soil ITS2 (26.4 MB, fits) and loses
//!   5% on soil 16S (39.2 MB, does not), which makes it a cache-capacity result
//!   with an expiry date.
//!
//! Tune `--distance` (default 8): too short and the line has not arrived; too
//! long and it is evicted before use, or the prefetches themselves saturate the
//! load/store queues.
//!
//! ## Read this before believing any number
//!
//! `examples/pval_cost.rs` reported ~7 ns/repricing against production's ~57,
//! because its touched set fit caches the real traversal misses. **This
//! benchmark has the same failure mode on a laptop.** Run it on the benchmark
//! node, at the real `--raws`, and check `raw_current` against the production
//! figure (`p_update` seconds / `raws repriced`) before reading the others. If
//! `raw_current` is far below production, the whole run is modelling the wrong
//! machine and the arm ordering is all it can offer.
//!
//! ```text
//! cargo build --profile release-native --example pval_layout
//! numactl --interleave=all ./target/release-native/examples/pval_layout \
//!     --raws 825214 --visit 100735 --rounds 40
//! ```

use std::hint::black_box;
use std::time::Instant;

/// Today's `Raw`: 152 bytes, align 8, hot fields at offset 104.
#[derive(Clone, Copy)]
#[repr(C)]
struct RawCurrent {
    _cold: [u8; 104],
    p: f64,
    lambda: f64,
    hamming: u32,
    reads: u32,
    prior: bool,
    _tail: [u8; 7],
}

/// Hot fields first, whole struct on a 64-byte boundary.
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct RawAligned {
    p: f64,
    lambda: f64,
    hamming: u32,
    reads: u32,
    prior: bool,
    _cold: [u8; 143],
}

/// What `Raw` still holds under `dense_split`.
#[derive(Clone, Copy)]
#[repr(C)]
struct RawSplit {
    _cold: [u8; 136],
    reads: u32,
    prior: bool,
    _tail: [u8; 11],
}

/// `Comparison`, dense.
#[derive(Clone, Copy)]
#[repr(C)]
struct Comp {
    i: u32,
    index: u32,
    lambda: f64,
    hamming: u32,
    _pad: u32,
}

/// Every hot field in one packed record — the ceiling.
#[derive(Clone, Copy)]
#[repr(C)]
struct Packed {
    p: f64,
    lambda: f64,
    hamming: u32,
    reads: u32,
    prior: bool,
    _pad: [u8; 7],
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
    let nraw: usize = arg("--raws", 825_214);
    let visit: usize = arg("--visit", 100_735);
    let rounds: usize = arg("--rounds", 40);

    // Member lists are arbitrary raw indices in insertion order — effectively
    // random with respect to memory — and, crucially, **a different subset each
    // round**: which clusters are dirty changes as the partition evolves.
    //
    // Two instrument failures got us here, both of which made the arms measure
    // cache hits instead of the gather they exist to price:
    //
    // 1. Members were `(i * 7919) % nraw`, described in a comment as "scattered
    //    but correlated". A large odd stride is a CONSTANT stride — the most
    //    prefetchable pattern there is.
    // 2. With that fixed, one fixed member list was reused every round. Visiting
    //    the same 100,735 of 825,214 raws touches ~6.4 MB of distinct lines,
    //    which is L3-resident after the first round; the 112 MB array never
    //    mattered. `raw_current` came back *faster* than `stride_control`, which
    //    is impossible if misses are being paid — and is exactly what the
    //    control exists to say.
    //
    // So: a fresh random sample per round, precomputed so generation cost stays
    // out of the timed region.
    let mut rng = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        // xorshift64*, reproducible without a dependency.
        rng ^= rng >> 12;
        rng ^= rng << 25;
        rng ^= rng >> 27;
        rng.wrapping_mul(0x2545_F491_4F6C_DD1D)
    };
    let rounds_members: Vec<Vec<usize>> = (0..rounds)
        .map(|_| {
            (0..visit)
                .map(|_| (next() % nraw as u64) as usize)
                .collect()
        })
        .collect();

    // Positive control: the same arm over a constant stride, i.e. fully
    // prefetchable. A null is evidence only once the instrument is shown able to
    // return something else, so if `stride_control` is not markedly FASTER than
    // `raw_current`, this machine is not exhibiting miss costs at this
    // working-set size and no layout arm here means anything.
    let stride: Vec<usize> = (0..visit).map(|i| (i * 7919) % nraw).collect();

    let cur = vec![
        RawCurrent {
            _cold: [0; 104],
            p: 0.0,
            lambda: 1e-6,
            hamming: 3,
            reads: 5,
            prior: false,
            _tail: [0; 7]
        };
        nraw
    ];
    let ali = vec![
        RawAligned {
            p: 0.0,
            lambda: 1e-6,
            hamming: 3,
            reads: 5,
            prior: false,
            _cold: [0; 143]
        };
        nraw
    ];
    let spl = vec![
        RawSplit {
            _cold: [0; 136],
            reads: 5,
            prior: false,
            _tail: [0; 11]
        };
        nraw
    ];
    let mut sp_p = vec![0.0f64; nraw];
    let sp_c = vec![
        Comp {
            i: 0,
            index: 0,
            lambda: 1e-6,
            hamming: 3,
            _pad: 0
        };
        nraw
    ];
    let mut pk = vec![
        Packed {
            p: 0.0,
            lambda: 1e-6,
            hamming: 3,
            reads: 5,
            prior: false,
            _pad: [0; 7]
        };
        nraw
    ];

    let mb = |b: usize| b as f64 / 1e6;
    println!(
        "raws={nraw} visited/round={visit} rounds={rounds}\n\
         resident hot state:  current {:.1} MB   aligned {:.1} MB   \
         split {:.1}+{:.1}+{:.1} MB   packed {:.1} MB\n\
         (per-CCD L3 on an EPYC 7713 is 32 MB — whether the hot state fits is \
         the question)\n\
         a fresh random sample of {visit} raws per round, so lines are not \
         reused across rounds\n",
        mb(nraw * size_of::<RawCurrent>()),
        mb(nraw * size_of::<RawAligned>()),
        mb(nraw * size_of::<RawSplit>()),
        mb(nraw * 8),
        mb(nraw * size_of::<Comp>()),
        mb(nraw * size_of::<Packed>()),
    );

    let report = |name: &str, d: std::time::Duration| {
        println!(
            "{name:<14} {:>8.3}s  {:>7.1} ns/repricing",
            d.as_secs_f64(),
            d.as_secs_f64() / (rounds * visit) as f64 * 1e9
        );
    };

    // Each arm mirrors the real loop: read the four fields, take the singleton
    // early exit ~94% of the time, write `p` back. Reads go through `black_box`
    // so nothing vectorises — the first draft of `store_bandwidth` reported an
    // impossible 574 GB/s for want of exactly that.
    let mut cur_p = vec![0.0f64; nraw];
    let t = Instant::now();
    for r in 0..rounds {
        for &m in &rounds_members[r] {
            let r = &cur[m];
            let (reads, prior, h, l) = (
                black_box(r.reads),
                black_box(r.prior),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            cur_p[m] = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&cur_p);
    }
    report("raw_current", t.elapsed());

    let mut ctl_p = vec![0.0f64; nraw];
    let t = Instant::now();
    for _ in 0..rounds {
        for &m in &stride {
            let r = &cur[m];
            let (reads, prior, h, l) = (
                black_box(r.reads),
                black_box(r.prior),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            ctl_p[m] = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&ctl_p);
    }
    report("stride_control", t.elapsed());

    let mut ali_p = vec![0.0f64; nraw];
    let t = Instant::now();
    for r in 0..rounds {
        for &m in &rounds_members[r] {
            let r = &ali[m];
            let (reads, prior, h, l) = (
                black_box(r.reads),
                black_box(r.prior),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            ali_p[m] = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&ali_p);
    }
    report("raw_aligned", t.elapsed());

    let t = Instant::now();
    for r in 0..rounds {
        for &m in &rounds_members[r] {
            let (reads, prior) = (black_box(spl[m].reads), black_box(spl[m].prior));
            let c = &sp_c[m];
            let (h, l) = (black_box(c.hamming), black_box(c.lambda));
            sp_p[m] = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&sp_p);
    }
    report("dense_split", t.elapsed());

    let t = Instant::now();
    for r in 0..rounds {
        for &m in &rounds_members[r] {
            let r = &pk[m];
            let (reads, prior, h, l) = (
                black_box(r.reads),
                black_box(r.prior),
                black_box(r.hamming),
                black_box(r.lambda),
            );
            pk[m].p = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&pk);
    }
    report("dense_all", t.elapsed());

    // Software prefetch, D iterations ahead. `Raw`'s hot fields span offsets
    // 104..146 — two lines — so both are requested; prefetching only the first
    // would leave half the stall in place.
    let dist: usize = arg("--distance", 8);
    // The prefetch is x86_64-only. On any other target the arm degrades to
    // `raw_current` plus a bounds check and would report a convincing null —
    // the fourth instrument in this issue that could not exhibit the effect it
    // was built to measure. Say so instead.
    let pf_active = cfg!(target_arch = "x86_64");
    if !pf_active {
        println!(
            "  NOTE: prefetch arm is INACTIVE on this target \
                  (x86_64 only); its number is raw_current, not a result"
        );
    }
    let mut pf_p = vec![0.0f64; nraw];
    let t = Instant::now();
    for r in 0..rounds {
        let ms = &rounds_members[r];
        for (i, &m) in ms.iter().enumerate() {
            #[cfg(target_arch = "x86_64")]
            if let Some(&ahead) = ms.get(i + dist) {
                // SAFETY: prefetch has no architectural effect — an out-of-range
                // or unmapped address is a no-op, never a fault. The pointer is
                // in-bounds regardless, since `ahead` indexes `cur`.
                unsafe {
                    let p = cur.as_ptr().add(ahead) as *const i8;
                    std::arch::x86_64::_mm_prefetch(p.add(104), std::arch::x86_64::_MM_HINT_T0);
                    std::arch::x86_64::_mm_prefetch(p.add(140), std::arch::x86_64::_MM_HINT_T0);
                }
            }
            #[cfg(not(target_arch = "x86_64"))]
            let _ = ms.get(i + dist);
            let rr = &cur[m];
            let (reads, prior, h, l) = (
                black_box(rr.reads),
                black_box(rr.prior),
                black_box(rr.hamming),
                black_box(rr.lambda),
            );
            pf_p[m] = if reads == 1 && !prior {
                1.0
            } else {
                l + h as f64
            };
        }
        black_box(&pf_p);
    }
    report(
        &format!(
            "prefetch(d={dist}){}",
            if pf_active { "" } else { " [INACTIVE]" }
        ),
        t.elapsed(),
    );

    println!(
        "\nstride_control is the SAME arm over a constant stride, i.e. fully\n\
         prefetchable. raw_current must be markedly slower than it; if the two\n\
         are close, this machine is not exhibiting miss costs at this working-set\n\
         size and no layout arm here means anything.\n\
         Then check raw_current against production: p_update seconds / raws\n\
         repriced was 56.6-65.4 ns on the four measured arms.\n\
         prefetch is compared against raw_current, not against dense_all: it is\n\
         the same layout, so any gain is latency hiding and carries no refactor\n\
         and no size ceiling. Sweep --distance before concluding it does nothing."
    );
}
