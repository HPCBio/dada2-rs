//! #295 gate: how tight is the presence-bitmap bound in front of `kmer_dist8`?
//!
//! For two raws with `D` shared distinct k-mers and `excess = nkmers - distinct`:
//!
//! ```text
//! sum(min(a_i, b_i)) <= D + min(excess_a, excess_b)
//! ```
//!
//! so `1 - (D + min excess) / scale` is a lower bound on the k-mer distance,
//! and a pair whose bound already exceeds the cutoff can be rejected without
//! reading the count vectors. This measures, offline, what share of the pairs
//! the full screen rejects the bound also rejects, and checks that it never
//! rejects a pair the full screen passes (the bound must be exact). Counts are
//! exact (16-bit), so the 8-bit overflow fallback cannot skew the result.
//!
//! As in `batch_nw_bench`, the most abundant uniques act as centres and every
//! unique is screened against each.
//!
//! ```text
//! cargo run --release --example kmer_bound_check -- --k 7 --centers 100 *.fastq.gz
//! ```
//! Flags: --k (5), --cutoff (0.42), --centers (100), --max-uniques (all),
//! --threads (all cores), --time (also time the screen with and without the
//! bound, over dense vectors as production holds them for k <= 7), --reps (3).
//!
//! `--time` runs two arms over every centre x raw pair, interleaved over reps,
//! best rep reported:
//! - `full`: `kmer_dist8` on every pair, as the screen does today.
//! - `bound`: `popcount(bitmap_a & bitmap_b)` first; `kmer_dist8` only where
//!   the bound cannot reject. The bitmap test uses POPCNT, so on x86-64 it runs
//!   an AVX2+POPCNT build when the CPU has both, as production would.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use rayon::prelude::*;

struct Args {
    k: usize,
    cutoff: f64,
    centers: usize,
    max_uniques: usize,
    threads: usize,
    time: bool,
    reps: usize,
    files: Vec<String>,
}

fn parse_args() -> Args {
    let mut a = Args {
        k: 5,
        cutoff: 0.42,
        centers: 100,
        max_uniques: usize::MAX,
        threads: 0,
        time: false,
        reps: 3,
        files: vec![],
    };
    let mut it = std::env::args().skip(1);
    while let Some(x) = it.next() {
        let mut val = || it.next().expect("flag needs a value");
        match x.as_str() {
            "--k" => a.k = val().parse().unwrap(),
            "--cutoff" => a.cutoff = val().parse().unwrap(),
            "--centers" => a.centers = val().parse().unwrap(),
            "--max-uniques" => a.max_uniques = val().parse().unwrap(),
            "--threads" => a.threads = val().parse().unwrap(),
            "--time" => a.time = true,
            "--reps" => a.reps = val().parse().unwrap(),
            f => a.files.push(f.to_string()),
        }
    }
    assert!(!a.files.is_empty(), "give at least one FASTQ");
    a
}

/// Uniques across all files, most abundant first (as `batch_nw_bench`).
fn load_uniques(files: &[String]) -> Vec<Vec<u8>> {
    let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
    for f in files {
        let fh = std::fs::File::open(f).unwrap_or_else(|e| panic!("{f}: {e}"));
        let rd: Box<dyn std::io::Read> = if f.ends_with(".gz") {
            Box::new(flate2::read::MultiGzDecoder::new(fh))
        } else {
            Box::new(fh)
        };
        for (n, line) in BufReader::new(rd).lines().enumerate() {
            let line = line.unwrap();
            if n % 4 != 1 {
                continue;
            }
            let s: Option<Vec<u8>> = line
                .trim_end()
                .bytes()
                .map(|c| match c {
                    b'A' => Some(0),
                    b'C' => Some(1),
                    b'G' => Some(2),
                    b'T' => Some(3),
                    _ => None,
                })
                .collect();
            if let Some(s) = s {
                *counts.entry(s).or_default() += 1;
            }
        }
    }
    let mut v: Vec<_> = counts.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    v.into_iter().map(|(s, _)| s).collect()
}

/// Sorted `(k-mer index, exact count)` for one sequence (bases 0..=3).
struct Profile {
    len: usize,
    kv: Vec<(u32, u16)>,
    /// `nkmers - distinct`.
    excess: u32,
}

fn profile(seq: &[u8], k: usize) -> Profile {
    let mut idx: Vec<u32> = seq
        .windows(k)
        .map(|w| w.iter().fold(0u32, |acc, &b| acc * 4 + b as u32))
        .collect();
    let nkmers = idx.len() as u32;
    idx.sort_unstable();
    let mut kv: Vec<(u32, u16)> = Vec::new();
    for i in idx {
        match kv.last_mut() {
            Some((j, c)) if *j == i => *c += 1,
            _ => kv.push((i, 1)),
        }
    }
    Profile {
        len: seq.len(),
        excess: nkmers - kv.len() as u32,
        kv,
    }
}

/// `(sum of min counts, shared distinct k-mers)`.
fn overlap(a: &Profile, b: &Profile) -> (u32, u32) {
    let (mut i, mut j, mut dot, mut shared) = (0, 0, 0u32, 0u32);
    while i < a.kv.len() && j < b.kv.len() {
        let (ai, ac) = a.kv[i];
        let (bi, bc) = b.kv[j];
        if ai < bi {
            i += 1;
        } else if ai > bi {
            j += 1;
        } else {
            dot += ac.min(bc) as u32;
            shared += 1;
            i += 1;
            j += 1;
        }
    }
    (dot, shared)
}

#[derive(Default, Clone, Copy)]
struct Tally {
    pairs: u64,
    pass: u64,
    reject: u64,
    bound_reject: u64,
    /// Pairs the bound rejects but the full screen passes: must stay zero.
    unsound: u64,
}

impl std::ops::Add for Tally {
    type Output = Tally;
    fn add(self, o: Tally) -> Tally {
        Tally {
            pairs: self.pairs + o.pairs,
            pass: self.pass + o.pass,
            reject: self.reject + o.reject,
            bound_reject: self.bound_reject + o.bound_reject,
            unsound: self.unsound + o.unsound,
        }
    }
}

fn main() {
    let a = parse_args();
    if a.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(a.threads)
            .build_global()
            .unwrap();
    }
    let k = a.k;
    let mut seqs = load_uniques(&a.files);
    let total = seqs.len();
    seqs.truncate(a.max_uniques);
    let profs: Vec<Profile> = seqs.par_iter().map(|s| profile(s, k)).collect();
    let nc = a.centers.min(profs.len());

    let t = (0..nc)
        .into_par_iter()
        .map(|c| {
            let pc = &profs[c];
            let mut t = Tally::default();
            for (r, pr) in profs.iter().enumerate() {
                if r == c {
                    continue;
                }
                let scale = (pc.len.min(pr.len) - k + 1) as f64;
                let (dot, shared) = overlap(pc, pr);
                let dist = 1.0 - dot as f64 / scale;
                let bound = 1.0 - (shared + pc.excess.min(pr.excess)) as f64 / scale;
                let pass = dist <= a.cutoff;
                let bound_rej = bound > a.cutoff;
                t.pairs += 1;
                if pass {
                    t.pass += 1;
                    t.unsound += bound_rej as u64;
                } else {
                    t.reject += 1;
                    t.bound_reject += bound_rej as u64;
                }
            }
            t
        })
        .reduce(Tally::default, |x, y| x + y);

    let mut lens: Vec<usize> = profs.iter().map(|p| p.len).collect();
    lens.sort_unstable();
    let mean_excess = profs.iter().map(|p| p.excess as f64).sum::<f64>() / profs.len() as f64;
    let mean_nk = profs.iter().map(|p| (p.len - k + 1) as f64).sum::<f64>() / profs.len() as f64;
    println!(
        "files={} uniques={total} (used {}) len[min/med/max]={}/{}/{} k={k} cutoff={} centers={nc}",
        a.files.len(),
        profs.len(),
        lens[0],
        lens[lens.len() / 2],
        lens[lens.len() - 1],
        a.cutoff
    );
    println!(
        "mean nkmers={mean_nk:.1} mean excess={mean_excess:.1} ({:.1}% of nkmers)",
        100.0 * mean_excess / mean_nk
    );
    println!(
        "pairs={} pass={} ({:.2}%) reject={} ({:.2}%)",
        t.pairs,
        t.pass,
        100.0 * t.pass as f64 / t.pairs as f64,
        t.reject,
        100.0 * t.reject as f64 / t.pairs as f64
    );
    println!(
        "bound rejects {} of {} screen rejections: tightness {:.3}%; full screen still needed on {:.2}% of pairs",
        t.bound_reject,
        t.reject,
        100.0 * t.bound_reject as f64 / t.reject.max(1) as f64,
        100.0 * (t.pairs - t.bound_reject) as f64 / t.pairs as f64
    );
    println!("unsound (bound rejects a passing pair): {}", t.unsound);
    assert_eq!(
        t.unsound, 0,
        "the bound rejected a pair the full screen passes"
    );

    if a.time {
        time_arms(&seqs, &profs, nc, k, a.cutoff, a.reps, t);
    }
}

/// What the screen would hold per raw: the dense 8-bit counts it reads today,
/// plus the presence bitmap and `excess` the bound needs.
struct Screen {
    len: usize,
    kmer8: Vec<u8>,
    bits: Vec<u64>,
    excess: u32,
}

#[inline(always)]
fn shared_body(a: &[u64], b: &[u64]) -> u32 {
    a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum()
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
fn shared_avx2(a: &[u64], b: &[u64]) -> u32 {
    shared_body(a, b)
}

fn shared(a: &[u64], b: &[u64], fast: bool) -> u32 {
    #[cfg(target_arch = "x86_64")]
    if fast {
        // SAFETY: `fast` is set only when the CPU has AVX2 and POPCNT.
        return unsafe { shared_avx2(a, b) };
    }
    let _ = fast;
    shared_body(a, b)
}

fn time_arms(
    seqs: &[Vec<u8>],
    profs: &[Profile],
    nc: usize,
    k: usize,
    cutoff: f64,
    reps: usize,
    t: Tally,
) {
    use dada2_rs::kmers::kmer_dist8;
    #[cfg(target_arch = "x86_64")]
    let fast = std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("popcnt");
    #[cfg(not(target_arch = "x86_64"))]
    let fast = false;
    let words = (1usize << (2 * k)).div_ceil(64);
    let screens: Vec<Screen> = seqs
        .par_iter()
        .zip(profs)
        .map(|(s, p)| {
            let enc: Vec<u8> = s.iter().map(|&b| b + 1).collect();
            let mut bits = vec![0u64; words];
            for &(i, _) in &p.kv {
                bits[i as usize / 64] |= 1 << (i % 64);
            }
            Screen {
                len: s.len(),
                kmer8: dada2_rs::kmers::assign_kmer8(&enc, k),
                bits,
                excess: p.excess,
            }
        })
        .collect();

    // Returns (passes, bound rejections), so both arms can be checked against
    // the tightness pass above.
    let run = |bound: bool| -> (u64, u64, f64) {
        let t0 = std::time::Instant::now();
        let (pass, brej) = (0..nc)
            .into_par_iter()
            .map(|c| {
                let sc = &screens[c];
                let (mut pass, mut brej) = (0u64, 0u64);
                for (r, sr) in screens.iter().enumerate() {
                    if r == c {
                        continue;
                    }
                    let scale = (sc.len.min(sr.len) - k + 1) as f64;
                    if bound {
                        let d = shared(&sc.bits, &sr.bits, fast);
                        if 1.0 - (d + sc.excess.min(sr.excess)) as f64 / scale > cutoff {
                            brej += 1;
                            continue;
                        }
                    }
                    let d8 = kmer_dist8(&sc.kmer8, sc.len, &sr.kmer8, sr.len, k);
                    // Overflow (-1) falls back to 16-bit in production; count
                    // it as a pass here, as that fallback is rare and exact.
                    pass += (d8 <= cutoff) as u64;
                }
                (pass, brej)
            })
            .reduce(|| (0, 0), |x, y| (x.0 + y.0, x.1 + y.1));
        (pass, brej, t0.elapsed().as_secs_f64())
    };
    let (mut best_full, mut best_bound) = (f64::INFINITY, f64::INFINITY);
    for _ in 0..reps {
        let (p, _, s) = run(false);
        assert_eq!(p, t.pass, "full arm disagrees with the exact pass count");
        best_full = best_full.min(s);
        let (p, b, s) = run(true);
        assert_eq!(p, t.pass, "bound arm disagrees with the exact pass count");
        assert_eq!(b, t.bound_reject, "bound arm rejected a different set");
        best_bound = best_bound.min(s);
    }
    let thr = rayon::current_num_threads() as f64;
    let per = |s: f64| s * thr * 1e9 / t.pairs as f64;
    println!(
        "time (best of {reps}, {} threads, popcnt={}): full {best_full:.3} s ({:.0} ns/pair busy)  bound {best_bound:.3} s ({:.0} ns/pair busy)  ratio {:.3}",
        thr,
        fast,
        per(best_full),
        per(best_bound),
        best_bound / best_full
    );
}
