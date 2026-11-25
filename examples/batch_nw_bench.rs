//! PROTOTYPE bench: single-pair vs inter-sequence batched banded NW on real
//! screen-passing pairs (`src/nwalign_batch.rs`).
//!
//! Mimics `b_compare`: the most abundant uniques act as cluster centers and
//! every unique is screened against each (k-mer distance ≤ cutoff, gapless
//! shortcut excluded). Pairs that pass are aligned by both kernels, so the
//! only difference between arms is the kernel. Batches are formed per center
//! from raws of equal length; a short final group is padded.
//!
//! Also reports the shape statistics that decide whether same-shape batching
//! is viable at all, and the gapless-shortcut share of passing pairs.
//!
//! ```text
//! cargo run --profile release-native --example batch_nw_bench -- \
//!     --k 5 --band 16 --threads 1 filtered/*.fastq.gz
//! ```
//! Flags: --k, --band, --cutoff (0.42), --centers (50), --max-uniques (20000),
//! --threads (1), --global (default ends-free), --reps (3), --arms (all;
//! comma list of single,L8,L16,L32,hybrid16,L16fill,L16alloc; `single` first for the
//! speedup column).

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use dada2_rs::kmers::{assign_kmer_order, assign_kmer8_sparse, kmer_dist8_sparse, kord_dist};
use dada2_rs::nwalign::{AlignBuffers, VectorizedAlignScores, align_vectorized_with_buf};
use dada2_rs::nwalign_batch::{
    BatchBuffers, align_batch_all, align_vectorized_batch, align_vectorized_batch_fill_only,
};
use rayon::prelude::*;

struct Args {
    k: usize,
    band: i32,
    cutoff: f64,
    centers: usize,
    max_uniques: usize,
    threads: usize,
    global: bool,
    reps: usize,
    arms: Option<Vec<String>>,
    files: Vec<String>,
}

fn parse_args() -> Args {
    let mut a = Args {
        k: 5,
        band: 16,
        cutoff: 0.42,
        centers: 50,
        max_uniques: 20_000,
        threads: 1,
        global: false,
        reps: 3,
        arms: None,
        files: vec![],
    };
    let mut it = std::env::args().skip(1);
    while let Some(x) = it.next() {
        let mut val = || it.next().expect("flag needs a value");
        match x.as_str() {
            "--k" => a.k = val().parse().unwrap(),
            "--band" => a.band = val().parse().unwrap(),
            "--cutoff" => a.cutoff = val().parse().unwrap(),
            "--centers" => a.centers = val().parse().unwrap(),
            "--max-uniques" => a.max_uniques = val().parse().unwrap(),
            "--threads" => a.threads = val().parse().unwrap(),
            "--reps" => a.reps = val().parse().unwrap(),
            "--global" => a.global = true,
            "--arms" => a.arms = Some(val().split(',').map(String::from).collect()),
            f => a.files.push(f.to_string()),
        }
    }
    assert!(!a.files.is_empty(), "give at least one FASTQ");
    a
}

fn encode(s: &[u8]) -> Option<Vec<u8>> {
    s.iter()
        .map(|&c| match c {
            b'A' => Some(1),
            b'C' => Some(2),
            b'G' => Some(3),
            b'T' => Some(4),
            _ => None,
        })
        .collect()
}

/// Uniques across all files, most abundant first.
fn load_uniques(files: &[String]) -> Vec<(Vec<u8>, u32)> {
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
            if n % 4 == 1
                && let Some(s) = encode(line.trim_end().as_bytes())
            {
                *counts.entry(s).or_default() += 1;
            }
        }
    }
    let mut v: Vec<_> = counts.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    v
}

struct Uniq {
    seq: Vec<u8>,
    kv: Vec<(u16, u8)>,
    kord: Vec<u16>,
}

/// Screen-passing, non-gapless pairs for one center, grouped by raw length.
type Groups = Vec<Vec<usize>>;

fn cells(l1: usize, l2: usize, band: i32) -> f64 {
    ((l1 + l2 + 1) * (band.max(0) as usize + 1)) as f64
}

fn main() {
    let a = parse_args();
    rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build_global()
        .unwrap();

    let mut uniq = load_uniques(&a.files);
    let total_uniques = uniq.len();
    uniq.truncate(a.max_uniques);
    let k = a.k;
    let raws: Vec<Uniq> = uniq
        .par_iter()
        .map(|(s, _)| Uniq {
            kv: assign_kmer8_sparse(s, k),
            kord: assign_kmer_order(s, k),
            seq: s.clone(),
        })
        .collect();
    let nc = a.centers.min(raws.len());

    // Screen.
    let screened: Vec<(Groups, usize, usize)> = (0..nc)
        .into_par_iter()
        .map(|c| {
            let rc = &raws[c];
            let mut by_len: HashMap<usize, Vec<usize>> = HashMap::new();
            let (mut pass, mut gapless) = (0usize, 0usize);
            for (r, rr) in raws.iter().enumerate() {
                if r == c {
                    continue;
                }
                let d = kmer_dist8_sparse(&rc.kv, rc.seq.len(), &rr.kv, rr.seq.len(), k);
                if d > a.cutoff || d < 0.0 {
                    continue;
                }
                pass += 1;
                let kod = kord_dist(&rc.kord, rc.seq.len(), &rr.kord, rr.seq.len(), k);
                if kod == d {
                    gapless += 1;
                    continue;
                }
                by_len.entry(rr.seq.len()).or_default().push(r);
            }
            let mut g: Groups = by_len.into_values().collect();
            g.sort_by_key(|v| std::cmp::Reverse(v.len()));
            (g, pass, gapless)
        })
        .collect();

    let comparisons = nc * (raws.len() - 1);
    let pass: usize = screened.iter().map(|s| s.1).sum();
    let gapless: usize = screened.iter().map(|s| s.2).sum();
    let aligned: usize = screened.iter().flat_map(|s| &s.0).map(Vec::len).sum();
    let ngroups: usize = screened.iter().map(|s| s.0.len()).sum();
    let mut lens: Vec<usize> = raws.iter().map(|u| u.seq.len()).collect();
    lens.sort_unstable();

    println!(
        "files={} uniques={} (used {}) len[min/med/max]={}/{}/{} centers={nc} k={k} band={} ends={} threads={}",
        a.files.len(),
        total_uniques,
        raws.len(),
        lens[0],
        lens[lens.len() / 2],
        lens[lens.len() - 1],
        a.band,
        if a.global { "global" } else { "free" },
        a.threads
    );
    println!(
        "comparisons={comparisons} pass={pass} ({:.1}%) gapless={gapless} ({:.1}% of pass) aligned={aligned}",
        100.0 * pass as f64 / comparisons as f64,
        100.0 * gapless as f64 / pass.max(1) as f64
    );
    println!(
        "shape groups={ngroups} mean pairs/group={:.1}",
        aligned as f64 / ngroups.max(1) as f64
    );
    for l in [8usize, 16, 32] {
        let full: usize = screened
            .iter()
            .flat_map(|s| &s.0)
            .map(|g| g.len() / l * l)
            .sum();
        let lanes: usize = screened
            .iter()
            .flat_map(|s| &s.0)
            .map(|g| g.len().div_ceil(l) * l)
            .sum();
        println!(
            "  L={l:>2}: pairs in full batches {:.1}%  lane occupancy (padded) {:.1}%",
            100.0 * full as f64 / aligned.max(1) as f64,
            100.0 * aligned as f64 / lanes.max(1) as f64
        );
    }

    let edge = if a.global { -8 } else { 0 };
    let scores = VectorizedAlignScores {
        match_score: 5,
        mismatch: -4,
        gap_p: -8,
        begin_gap_p: edge,
        end_gap_p: edge,
        band: a.band,
    };
    let total_cells: f64 = screened
        .iter()
        .enumerate()
        .flat_map(|(c, s)| s.0.iter().map(move |g| (c, g)))
        .map(|(c, g)| g.len() as f64 * cells(raws[c].seq.len(), raws[g[0]].seq.len(), a.band))
        .sum();

    // One arm = one pass over all centers; per-center work is the parallel unit.
    let run = |arm: &str| -> f64 {
        let t = Instant::now();
        let n: usize = screened
            .par_iter()
            .enumerate()
            .map_init(
                || (AlignBuffers::new(), BatchBuffers::new(), Vec::new()),
                |(sb, bb, out), (c, s)| {
                    let center = &raws[c].seq;
                    let mut n = 0usize;
                    for g in &s.0 {
                        let pairs: Vec<(&[u8], &[u8])> = g
                            .iter()
                            .map(|&r| (center.as_slice(), raws[r].seq.as_slice()))
                            .collect();
                        let mut nb = 0usize;
                        let mut sink = |x: &[u8], _: &[u8]| nb += x.len();
                        match arm {
                            "single" => {
                                for &(x, y) in &pairs {
                                    align_vectorized_with_buf(x, y, &scores, sb);
                                    n += sb.al0.len();
                                }
                            }
                            "L8" => align_batch_all::<8>(&pairs, &scores, bb, &mut sink),
                            "L16" => align_batch_all::<16>(&pairs, &scores, bb, &mut sink),
                            "L32" => align_batch_all::<32>(&pairs, &scores, bb, &mut sink),
                            "L16alloc" => {
                                // The first cluster run's harness: lane buffers moved
                                // out per batch and held for the whole group. Kept as
                                // a control for the allocator-bound collapse.
                                for chunk in pairs.chunks(16) {
                                    let mut lanes = [chunk[0]; 16];
                                    lanes[..chunk.len()].copy_from_slice(chunk);
                                    align_vectorized_batch::<16>(&lanes, &scores, bb);
                                    out.extend(bb.al.drain(..chunk.len()));
                                }
                                n += out
                                    .iter()
                                    .map(|o: &(Vec<u8>, Vec<u8>)| o.0.len())
                                    .sum::<usize>();
                                out.clear();
                            }
                            "L16fill" => {
                                for chunk in pairs.chunks(16) {
                                    let mut lanes = [chunk[0]; 16];
                                    lanes[..chunk.len()].copy_from_slice(chunk);
                                    align_vectorized_batch_fill_only::<16>(&lanes, &scores, bb);
                                }
                            }
                            "hybrid16" => {
                                // Full 16-lane batches batched, the remainder single.
                                let full = pairs.len() / 16 * 16;
                                align_batch_all::<16>(&pairs[..full], &scores, bb, &mut sink);
                                for &(x, y) in &pairs[full..] {
                                    align_vectorized_with_buf(x, y, &scores, sb);
                                    n += sb.al0.len();
                                }
                            }
                            _ => unreachable!(),
                        }
                        n += nb;
                    }
                    n
                },
            )
            .sum();
        std::hint::black_box(n);
        t.elapsed().as_secs_f64()
    };

    // Identity spot-check on real pairs before timing anything.
    {
        let mut sb = AlignBuffers::new();
        let mut bb = BatchBuffers::new();
        let mut out: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        let mut checked = 0usize;
        for (c, s) in screened.iter().enumerate().take(5) {
            for g in &s.0 {
                let pairs: Vec<(&[u8], &[u8])> = g
                    .iter()
                    .map(|&r| (raws[c].seq.as_slice(), raws[r].seq.as_slice()))
                    .collect();
                out.clear();
                align_batch_all::<16>(&pairs, &scores, &mut bb, |x, y| {
                    out.push((x.to_vec(), y.to_vec()))
                });
                for (&(x, y), got) in pairs.iter().zip(&out) {
                    align_vectorized_with_buf(x, y, &scores, &mut sb);
                    assert_eq!((&sb.al0, &sb.al1), (&got.0, &got.1), "identity");
                    checked += 1;
                }
            }
        }
        println!("identity check: {checked} real pairs byte-identical");
    }

    println!(
        "{:<9} {:>9} {:>11} {:>9} {:>8}",
        "arm", "wall_s", "ns/pair*", "ns/cell*", "speedup"
    );
    let all = [
        "single", "L8", "L16", "L32", "hybrid16", "L16fill", "L16alloc",
    ];
    let arms: Vec<&str> = match &a.arms {
        Some(v) => all
            .iter()
            .copied()
            .filter(|x| v.iter().any(|y| y == x))
            .collect(),
        None => all.to_vec(),
    };
    let mut best = vec![f64::INFINITY; arms.len()];
    // Interleave arms across reps so drift hits every arm alike.
    for _ in 0..a.reps {
        for (i, arm) in arms.iter().enumerate() {
            best[i] = best[i].min(run(arm));
        }
    }
    for (i, arm) in arms.iter().enumerate() {
        let busy = best[i] * a.threads as f64 * 1e9;
        println!(
            "{arm:<9} {:>9.3} {:>11.0} {:>9.3} {:>7.2}x",
            best[i],
            busy / aligned as f64,
            busy / total_cells,
            best[0] / best[i]
        );
    }
    println!(
        "* per-thread busy (wall x threads); cells exclude padded lanes; best of {} reps",
        a.reps
    );
}
