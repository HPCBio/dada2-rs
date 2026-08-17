//! Isolate the k-mer screen's throughput as a function of pool size, thread
//! count, and what else is resident alongside it (issue #142).
//!
//! ## What this exists to settle
//!
//! On the #139 archives, three pools run the *identical* screen kernel on the
//! *identical* amount of data — all `k=5`, dense, 1024 B per vector — and
//! differ 3.8x in cost per comparison:
//!
//! | pool     | raws      | screen ns/comp | implied GB/s/thread |
//! |----------|-----------|----------------|---------------------|
//! | soil 16S | 1,225,523 | 495            | 2.07                |
//! | ITS2     | 825,214   | 1890           | 0.54                |
//! | MiSeq    | 272,574   | 1091           | 0.94                |
//!
//! Equal work, unequal time, and the *largest* pool is the *fastest* — which is
//! backwards for both cache-residency and TLB-reach explanations. So the spread
//! is in the memory system, not in the kernel, and it is reproducible without
//! any real sequence data: synthetic vectors with the observed fill statistics
//! stream exactly the same way.
//!
//! Enabling the #139 carry moves 16S from 495 to 1109 ns/comp — essentially
//! MiSeq's rate. The two arms differ in one memory-system respect: the rebuilt
//! arm allocates and frees `compmax`/`emax` (~39 MB) about 9,700 times, while
//! the carry arm allocates once and keeps it resident for the whole run.
//! `--carry-mode` reproduces that difference directly.
//!
//! ## Reading the output
//!
//! `ns/comp` is scaled by thread count (`wall x threads / comps`) so it is
//! comparable to the `[dada] compare split` line in a run log, which divides
//! summed worker-busy time by comparisons. `GB/s` is aggregate across threads.
//!
//! Two outcomes, and they point opposite ways:
//!
//! - **Throughput rises with pool size.** 16S's 495 ns/comp is the normal rate
//!   and the smaller pools are the anomaly. #142 inverts: the question stops
//!   being "recover 16S's loss" and becomes "why do the others never get
//!   there", which is worth more than the -47%.
//! - **Throughput is flat in pool size, and only `--carry-mode` moves it.** The
//!   effect is allocator/page behaviour, and the next stop is THP and
//!   fault-count instrumentation, not the screen.
//!
//! If neither reproduces, the effect is not in the screen kernel at all and the
//! whole-run A/B was measuring something downstream of it.
//!
//! ## Usage
//!
//! ```text
//! cargo build --release --example screen_bandwidth
//! numactl --interleave=all ./target/release/examples/screen_bandwidth \
//!     --raws 272574,825214,1225523 --threads 1,8,24,48,64 --rounds 20
//! ```
//!
//! Always under `numactl --interleave=all` on the EPYC node — see
//! `docs/findings/measuring-on-numa.md`. Output is TSV on stdout for archiving;
//! the run's own configuration is echoed to stderr.

use dada2_rs::kmers::KmerScreen;
use rayon::prelude::*;

/// Grain matching `cluster::par_max_len`'s default, so rayon splits the range
/// the way `b_compare_parallel` does. The access pattern is the thing under
/// test, so this must not drift from the real one.
const PAR_GRAIN: usize = 32;

/// Deterministic xorshift — no rand dependency, and identical vectors across
/// runs so two invocations are comparable.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Build one dense k-mer screen with the fill statistics the real pools show:
/// ~198 distinct k-mers out of 1024 buckets (19.3% of dense 4^k), counts
/// summing to `len - k + 1`.
///
/// Fill matters beyond footprint: `kmer_dist8` branches on `min == 255`, and a
/// vector dense enough to saturate would take the overflow path and measure
/// something else entirely. Amplicon data never does.
fn make_screen(rng: &mut Rng, k: usize, len: usize, distinct: usize) -> KmerScreen {
    let buckets = 4usize.pow(k as u32);
    let mut v = vec![0u8; buckets];
    let total = len - k + 1;
    // Place `distinct` occupied buckets, then distribute the remaining count
    // mass over them — mirrors a real sequence, where most k-mers appear once
    // and a few repeat.
    let mut placed = 0usize;
    while placed < distinct {
        let b = rng.below(buckets);
        if v[b] == 0 {
            v[b] = 1;
            placed += 1;
        }
    }
    let mut remaining = total.saturating_sub(distinct);
    while remaining > 0 {
        let b = rng.below(buckets);
        if v[b] > 0 && v[b] < 200 {
            v[b] += 1;
            remaining -= 1;
        }
    }
    KmerScreen::Dense(v)
}

/// The persistent working set whose residency separates the two #139 arms.
///
/// `compmax`/`emax` are per-raw f64 arrays touched in scattered order by the
/// shuffle. Allocating it is not enough to reproduce the arm — untouched pages
/// are never faulted in — so this writes across it in scattered order.
fn touch_carry(buf: &mut [f64], rng: &mut Rng) {
    let n = buf.len();
    if n == 0 {
        return;
    }
    // Touch *every* element, in scattered order. Coverage has to be total: a
    // partial touch would leave the churn arm holding fewer resident pages than
    // the resident arm, and the resulting difference would be a property of how
    // much each arm faulted in rather than of how long it held the allocation —
    // which is the only thing this is meant to isolate. The real rebuilt arm
    // reads all of `compmax`/`emax` on every build scan, so full coverage is
    // also the faithful choice.
    let stride = 1_000_003usize; // coprime with any realistic n: a full permutation
    let mut i = 0usize;
    for _ in 0..n {
        i = (i + stride) % n;
        buf[i] += 1.0;
    }
    std::hint::black_box(rng.next());
}

#[derive(Clone, Copy, PartialEq)]
enum CarryMode {
    /// No competing working set — the floor.
    None,
    /// Allocated once and kept alive across all rounds (the #139 carry arm).
    Resident,
    /// Allocated and freed every round (the rebuilt arm, ~9,700 cycles).
    Churn,
}

fn parse_list(s: &str) -> Vec<usize> {
    s.split(',')
        .filter(|t| !t.is_empty())
        .map(|t| t.parse().unwrap_or_else(|_| panic!("not a number: {t:?}")))
        .collect()
}

fn main() {
    let mut raws = vec![272_574usize, 825_214, 1_225_523];
    let mut threads = vec![1usize, 8, 24, 48, 64];
    let mut rounds = 20usize;
    let mut k = 5usize;
    let mut len = 250usize;
    let mut distinct = 198usize;
    let mut carry_mb = 39usize;
    let mut carry_mode = CarryMode::None;
    // Heap spacing between consecutive k-mer vectors. In a real run each raw
    // also owns `seq`, `qual` and `kord`, so the screen's stream touches 1024
    // useful bytes out of a much longer stride and wastes about half of every
    // cache line it pulls. The 16S log reports 1728.3 MB of k-mer vectors and
    // 540.9 MB of seq+qual over 1,225,523 raws — 1411 + 441 B/raw, of which
    // 1024 is the dense screen, leaving ~830 B of interleaved other data.
    // Allocating that padding reproduces the spacing; without it the benchmark
    // measures a best-case contiguous stream that the run never gets.
    let mut pad = 830usize;

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let val = |i: usize| -> String {
            args.get(i + 1)
                .unwrap_or_else(|| panic!("{} needs a value", args[i]))
                .clone()
        };
        match args[i].as_str() {
            "--raws" => raws = parse_list(&val(i)),
            "--threads" => threads = parse_list(&val(i)),
            "--rounds" => rounds = val(i).parse().expect("--rounds"),
            "--k" => k = val(i).parse().expect("--k"),
            "--len" => len = val(i).parse().expect("--len"),
            "--distinct" => distinct = val(i).parse().expect("--distinct"),
            "--carry-mb" => carry_mb = val(i).parse().expect("--carry-mb"),
            "--pad" => pad = val(i).parse().expect("--pad"),
            "--carry-mode" => {
                carry_mode = match val(i).as_str() {
                    "none" => CarryMode::None,
                    "resident" => CarryMode::Resident,
                    "churn" => CarryMode::Churn,
                    other => panic!("--carry-mode: none|resident|churn, got {other:?}"),
                }
            }
            "--help" | "-h" => {
                eprintln!(
                    "screen_bandwidth --raws N,N --threads T,T --rounds R \\\n  \
                     [--k 5] [--len 250] [--distinct 198] \\\n  \
                     [--pad 830] [--carry-mode none|resident|churn] [--carry-mb 39]"
                );
                return;
            }
            other => panic!("unknown flag {other:?} (try --help)"),
        }
        i += 2;
    }

    let bytes_per_vec = 4usize.pow(k as u32);
    eprintln!(
        "[screen] k={k} len={len} distinct={distinct}/{bytes_per_vec} buckets, {bytes_per_vec} B/vector, {rounds} rounds"
    );
    eprintln!(
        "[screen] heap stride: {} B per raw ({} B screen + {pad} B interleaved filler)",
        bytes_per_vec + pad,
        bytes_per_vec,
    );
    eprintln!(
        "[screen] carry mode: {}",
        match carry_mode {
            CarryMode::None => "none".to_string(),
            CarryMode::Resident => format!("resident, {carry_mb} MB held across all rounds"),
            CarryMode::Churn => format!("churn, {carry_mb} MB allocated+freed per round"),
        }
    );

    println!(
        "raws\tthreads\tMB_resident\tcomps\twall_s\tns_per_comp\tGB_per_s\tGB_per_s_per_thread"
    );

    for &nraw in &raws {
        // Each screen is its own heap allocation, matching `Raw::kmer8`
        // (`Option<KmerScreen>` per raw). A single contiguous slab would be a
        // different access pattern and would not reproduce the run.
        let mut rng = Rng(0x2545F4914F6CDD1D);
        // Interleave a per-raw filler allocation between the screens so their
        // heap spacing matches a real run's (see `pad` above). The fillers are
        // held for the whole measurement — freeing them would let the allocator
        // reuse the gaps and collapse the stride back to contiguous.
        let mut fillers: Vec<Vec<u8>> = Vec::with_capacity(if pad > 0 { nraw } else { 0 });
        let mut screens: Vec<KmerScreen> = Vec::with_capacity(nraw);
        for _ in 0..nraw {
            screens.push(make_screen(&mut rng, k, len, distinct));
            if pad > 0 {
                fillers.push(vec![0u8; pad]);
            }
        }
        std::hint::black_box(&fillers);
        let resident_mb = (nraw * bytes_per_vec) as f64 / (1024.0 * 1024.0);

        // A fixed center, as in a real compare call: it stays hot in L1 while
        // every other vector is streamed once.
        let center = &screens[0];
        let carry_len = carry_mb * 1024 * 1024 / std::mem::size_of::<f64>();

        for &nthreads in &threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build()
                .expect("rayon pool");

            let mut held: Option<Vec<f64>> = if carry_mode == CarryMode::Resident {
                let mut v = vec![0.0f64; carry_len];
                touch_carry(&mut v, &mut rng);
                Some(v)
            } else {
                None
            };

            // One untimed round first: the measurement is of steady-state
            // streaming, not of first-touch page faults on the vectors.
            pool.install(|| {
                let s: f64 = screens
                    .par_iter()
                    .with_max_len(PAR_GRAIN)
                    .map(|s| center.dist8(s, len, len, k))
                    .sum();
                std::hint::black_box(s);
            });

            let mut elapsed = std::time::Duration::ZERO;
            for _ in 0..rounds {
                if carry_mode == CarryMode::Churn {
                    let mut v = vec![0.0f64; carry_len];
                    touch_carry(&mut v, &mut rng);
                    std::hint::black_box(&v);
                    drop(v);
                } else if let Some(v) = held.as_mut() {
                    touch_carry(v, &mut rng);
                }

                let t = std::time::Instant::now();
                let sum = pool.install(|| {
                    screens
                        .par_iter()
                        .with_max_len(PAR_GRAIN)
                        .map(|s| center.dist8(s, len, len, k))
                        .sum::<f64>()
                });
                elapsed += t.elapsed();
                std::hint::black_box(sum);
            }

            let comps = (nraw * rounds) as f64;
            let wall = elapsed.as_secs_f64();
            let bytes = comps * bytes_per_vec as f64;
            // Scaled by threads so it is directly comparable to the run log's
            // `[dada] compare split` ns/comp, which divides summed worker-busy
            // time by comparisons.
            let ns_per_comp = wall * 1e9 * nthreads as f64 / comps;
            let gbs = bytes / wall / 1e9;
            println!(
                "{nraw}\t{nthreads}\t{resident_mb:.0}\t{comps:.0}\t{wall:.3}\t{ns_per_comp:.1}\t{gbs:.2}\t{:.3}",
                gbs / nthreads as f64
            );
            drop(held.take());
        }
    }
}
