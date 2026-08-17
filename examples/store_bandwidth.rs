//! Decompose `b_compare`'s serial store scan into its two memory accesses
//! (issue #147).
//!
//! ## What this exists to settle
//!
//! After #146 folded the serial reduction away, the store is the entire serial
//! block inside `b_compare` — 23-30% of `run_dada` — at a flat **18.4-19.6
//! ns/raw-visit** across a 5x range in pool size. Its commit path is *not* the
//! cost: stores are 0.102% of raw-visits (12,082,545 of 11,888,798,623), so
//! the loop is ~entirely a scan.
//!
//! The scan touches exactly two things per raw-visit:
//!
//! | | bytes used | access | traffic |
//! |---|---|---|---|
//! | `comps[index]` | 48 | sequential | 48 B |
//! | `b.raws[index].e_minmax` | 8 | **160 B stride** | 64 B (one line) |
//!
//! `Raw` is 160 bytes, so consecutive `e_minmax` reads land 2.5 cache lines
//! apart and every visit pulls a full line to use 8 bytes of it. ~112 B/visit
//! against 19.0 ns implies ~5.9 GB/s single-threaded, which is *low* for a
//! prefetchable constant stride on this hardware. That gap is the whole
//! question, and it decides which lever is worth building:
//!
//! - **bandwidth-bound** -> lever A (drop `CompCost` from the result element,
//!   48 B -> 16 B) is the cheap win and B is incremental.
//! - **latency-bound on the strided read** -> lever B (hoist `e_minmax` into a
//!   dense `Vec<f64>`) is the real win and A is a rounding error.
//!
//! Do not infer this from the production `ns/raw` figure alone. #146 falsified
//! a residency explanation for the *reduction* that looked just as clean.
//!
//! ## Why a synthetic loop is trustworthy here
//!
//! Unlike `screen_bandwidth` — which reproduced ordering but never magnitudes,
//! because it stood in for real alignment work — this loop has no compute to
//! stand in for. It is two array reads, a multiply, and a compare. The
//! variants below differ *only* in element size and access pattern, which is
//! precisely what is under test. Absolute ns should still be checked against
//! the production figure before any of it is believed.
//!
//! Every arm reads through `black_box`, so none of them vectorises. Without
//! that the single-array arms auto-vectorise and the paired arms cannot,
//! which measures SIMD width rather than access pattern — the first draft of
//! this benchmark reported 574 GB/s for `comps48`, which is impossible and is
//! how the mistake was caught. The `black_box` cost is a per-iteration
//! constant, identical across arms, so it inflates absolute ns and leaves the
//! between-arm differences intact.
//!
//! ## Variants
//!
//! | arm | comps | e_minmax | models |
//! |---|---|---|---|
//! | `both48` | 48 B | strided 160 B | **today** |
//! | `both16` | 16 B | strided 160 B | lever A |
//! | `dense48` | 48 B | dense 8 B | lever B alone |
//! | `dense16` | 16 B | dense 8 B | levers A+B |
//! | `comps48` | 48 B | *not read* | isolates the sequential stream |
//! | `comps16` | 16 B | *not read* | ditto, lever-A element |
//! | `stride` | *not read* | strided 160 B | isolates the strided read |
//!
//! `both48 - comps48` and `both48 - stride` bracket how much of the scan each
//! access owns. If `stride` alone is near `both48`, the strided read dominates
//! and lever A cannot help much.
//!
//! ## Invocation
//!
//! ```sh
//! cargo run --release --example store_bandwidth -- --raws 1225523 --rounds 20
//! ```
//!
//! Defaults reproduce soil 16S R1 (`nraw = 1,225,523`, store rate 0.102%).
//! Single-threaded by design: the store is serial, which is the entire reason
//! it matters.

use std::hint::black_box;
use std::time::Instant;

/// Mirrors `containers::Raw`'s 160-byte footprint. Only `e_minmax` is read by
/// the store scan; the rest exists to reproduce the stride, which is the point.
#[repr(C)]
#[derive(Clone)]
struct RawLike {
    e_minmax: f64,
    _pad: [u8; 152],
}

/// The map's result element as it is today: `(f64, u32, bool, CompCost)`.
#[repr(C)]
#[derive(Clone, Copy)]
struct Item48 {
    lambda: f64,
    hamming: u32,
    skipped: bool,
    cost: [u64; 4],
}

/// The same element with `CompCost` gated out (lever A): `(f64, u32, bool)`.
#[repr(C)]
#[derive(Clone, Copy)]
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

fn main() {
    let nraw: usize = arg("--raws", 1_225_523);
    let rounds: usize = arg("--rounds", 20);
    // Production stores on 0.102% of raw-visits. The branch is essentially
    // never taken, so it predicts almost perfectly — reproduce that rather
    // than a 50/50 branch, which would measure mispredicts instead of memory.
    let store_rate: f64 = arg("--store-rate", 0.00102);

    // Deterministic lambdas: `store_rate` of them clear the threshold.
    let mut lambdas = vec![0.0f64; nraw];
    let period = (1.0 / store_rate).round() as usize;
    for (i, l) in lambdas.iter_mut().enumerate() {
        *l = if i % period == 0 { 2.0 } else { 0.5 };
    }

    let items48: Vec<Item48> = lambdas
        .iter()
        .map(|&lambda| Item48 {
            lambda,
            hamming: 3,
            skipped: false,
            cost: [0; 4],
        })
        .collect();
    let items16: Vec<Item16> = lambdas
        .iter()
        .map(|&lambda| Item16 {
            lambda,
            hamming: 3,
            skipped: false,
        })
        .collect();
    let raws: Vec<RawLike> = (0..nraw)
        .map(|_| RawLike {
            e_minmax: 1.0,
            _pad: [0; 152],
        })
        .collect();
    let dense: Vec<f64> = vec![1.0; nraw];

    let mb = |b: usize| b as f64 / 1e6;
    println!(
        "nraw={nraw}  rounds={rounds}  items48={:.1} MB  items16={:.1} MB  \
         raws={:.1} MB  dense={:.1} MB",
        mb(nraw * 48),
        mb(nraw * 16),
        mb(nraw * 160),
        mb(nraw * 8),
    );
    println!(
        "{:<9}{:>12}{:>12}{:>14}",
        "arm", "ns/visit", "GB/s", "traffic B"
    );

    // Each arm returns (accumulator, bytes of traffic per visit) so a wrong
    // model shows up as a GB/s that cannot be right, rather than silently.
    let run = |name: &str, traffic: usize, f: &mut dyn FnMut() -> f64| {
        f(); // warm
        let t = Instant::now();
        let mut acc = 0.0;
        for _ in 0..rounds {
            acc += f();
        }
        let el = t.elapsed().as_secs_f64();
        let visits = (rounds * nraw) as f64;
        black_box(acc);
        println!(
            "{:<9}{:>12.2}{:>12.2}{:>14}",
            name,
            el / visits * 1e9,
            visits * traffic as f64 / el / 1e9,
            traffic
        );
    };

    run("both48", 112, &mut || {
        let mut hits = 0.0;
        for (i, it) in items48.iter().enumerate() {
            if black_box(it.lambda) > black_box(raws[i].e_minmax) {
                hits += 1.0;
            }
        }
        hits
    });
    run("both16", 80, &mut || {
        let mut hits = 0.0;
        for (i, it) in items16.iter().enumerate() {
            if black_box(it.lambda) > black_box(raws[i].e_minmax) {
                hits += 1.0;
            }
        }
        hits
    });
    run("dense48", 56, &mut || {
        let mut hits = 0.0;
        for (i, it) in items48.iter().enumerate() {
            if black_box(it.lambda) > black_box(dense[i]) {
                hits += 1.0;
            }
        }
        hits
    });
    run("dense16", 24, &mut || {
        let mut hits = 0.0;
        for (i, it) in items16.iter().enumerate() {
            if black_box(it.lambda) > black_box(dense[i]) {
                hits += 1.0;
            }
        }
        hits
    });
    run("comps48", 48, &mut || {
        let mut hits = 0.0;
        for it in items48.iter() {
            if black_box(it.lambda) > 1.0 {
                hits += 1.0;
            }
        }
        hits
    });
    run("comps16", 16, &mut || {
        let mut hits = 0.0;
        for it in items16.iter() {
            if black_box(it.lambda) > 1.0 {
                hits += 1.0;
            }
        }
        hits
    });
    run("stride", 64, &mut || {
        let mut hits = 0.0;
        for r in raws.iter() {
            if black_box(r.e_minmax) > 0.5 {
                hits += 1.0;
            }
        }
        hits
    });
}
