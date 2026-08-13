// Cost model for b_shuffle_converge's *build* phase (#124).
//
// Reproduces the build loop's access pattern in isolation and A/Bs data
// layouts. The build is:
//     for each cluster: for each stored comp (ascending raw index):
//         e = lambda * reads;  if e > emax[idx] { emax[idx]=e; compmax[idx]=comp }
// i.e. one sequential stream (the comp vec) plus a *sparse gather* into two
// per-raw arrays. The hypothesis under test is that cost is set by the gather's
// cache-line touches per comp, not by the comparison count.

use std::time::Instant;

#[derive(Clone, Copy, Default)]
struct Comparison {
    i: u32,
    index: u32,
    lambda: f64,
    hamming: u32,
}

// Cheap deterministic RNG so every variant sees identical data.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn f64(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Per-cluster comp vecs, ascending raw index, `density` of all raws present.
fn gen(nraw: usize, nclusters: usize, density: f64, seed: u64, decay: f64, zipf: f64) -> (Vec<Vec<Comparison>>, Vec<u32>) {
    let mut rng = Rng(seed);
    let mut clusters = Vec::with_capacity(nclusters);
    for ci in 0..nclusters {
        // cluster 0 holds every raw (b_compare runs unscreened on it)
        let d = if ci == 0 { 1.0 } else { density };
        let mut v = Vec::with_capacity((nraw as f64 * d) as usize + 16);
        for idx in 0..nraw {
            if ci == 0 || rng.f64() < d {
                v.push(Comparison {
                    i: ci as u32,
                    index: idx as u32,
                    // lambda decays hard with hamming distance; range spans
                    // many orders of magnitude, as in the real error model.
                    lambda: (rng.f64() * -decay).exp(),
                    hamming: (rng.f64() * 8.0) as u32,
                });
            }
        }
        clusters.push(v);
    }
    // Cluster read counts. zipf=0 -> uniform; zipf>0 -> power-law skew
    // (rank^-zipf), which is what real cluster abundances look like and which
    // *weakens* the max_reads bound by inflating its numerator.
    let reads: Vec<u32> = (0..nclusters)
        .map(|ci| {
            if zipf <= 0.0 {
                1 + (rng.next() % 100_000) as u32
            } else {
                let r = (ci + 1) as f64;
                (1.0 + 1e7 * r.powf(-zipf)) as u32
            }
        })
        .collect();
    (clusters, reads)
}

fn main() {
    let nraw: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(285_000);
    let nclusters: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(400);
    let density: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.082);
    let reps: usize = std::env::args().nth(4).and_then(|s| s.parse().ok()).unwrap_or(5);
    let decay: f64 = std::env::args().nth(5).and_then(|s| s.parse().ok()).unwrap_or(30.0);
    let zipf: f64 = std::env::args().nth(6).and_then(|s| s.parse().ok()).unwrap_or(0.0);

    let (clusters, reads) = gen(nraw, nclusters, density, 0x2545F4914F6CDD1D, decay, zipf);
    let ncomp: usize = clusters.iter().map(|c| c.len()).sum();
    println!(
        "nraw={nraw} nclusters={nclusters} density={density} decay={decay} zipf={zipf} comps={:.2}e6  \
         (comp stream {:.1} MB)",
        ncomp as f64 / 1e6,
        (ncomp * 24) as f64 / 1e6,
    );
    println!(
        "  per-raw map footprints: A emax+compmax = {:.1} MB | B fused = {:.1} MB | C emax+best_ci = {:.1} MB",
        (nraw * 32) as f64 / 1e6,
        (nraw * 16) as f64 / 1e6,
        (nraw * 12) as f64 / 1e6,
    );

    // ---- SoA copies of the same data, for the layout variants -------------
    let soa: Vec<(Vec<f64>, Vec<u32>)> = clusters
        .iter()
        .map(|c| (c.iter().map(|x| x.lambda).collect(), c.iter().map(|x| x.index).collect()))
        .collect();

    let mut results: Vec<(&str, f64, u64)> = Vec::new();

    // ---- A: current — AoS 24B comps, emax f64 + compmax Comparison(24B) ----
    {
        let mut best = f64::MAX;
        let mut check = 0u64;
        for _ in 0..reps {
            let mut emax = vec![f64::NEG_INFINITY; nraw];
            let mut compmax = vec![Comparison::default(); nraw];
            let t = Instant::now();
            for (ci, comps) in clusters.iter().enumerate() {
                let r = reads[ci] as f64;
                for comp in comps {
                    let idx = comp.index as usize;
                    let e = comp.lambda * r;
                    if e > emax[idx] {
                        emax[idx] = e;
                        compmax[idx] = *comp;
                    }
                }
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = compmax.iter().map(|c| c.i as u64).sum();
        }
        results.push(("A current: AoS24 comp, emax f64 + compmax 24B", best, check));
    }

    // ---- B: fused per-raw map (f64 emax, u32 best_ci) = 16B, AoS comps ----
    // compmax's other fields (index is the key itself; lambda/hamming) are only
    // needed for raws that actually move (~550 of ~9.6e6 comps), so they can be
    // recovered at move time instead of being written during the scan.
    {
        let mut best = f64::MAX;
        let mut check = 0u64;
        for _ in 0..reps {
            let mut map: Vec<(f64, u32)> = vec![(f64::NEG_INFINITY, 0); nraw];
            let t = Instant::now();
            for (ci, comps) in clusters.iter().enumerate() {
                let r = reads[ci] as f64;
                let cu = ci as u32;
                for comp in comps {
                    let idx = comp.index as usize;
                    let e = comp.lambda * r;
                    let slot = unsafe { map.get_unchecked_mut(idx) };
                    if e > slot.0 {
                        *slot = (e, cu);
                    }
                }
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = map.iter().map(|m| m.1 as u64).sum();
        }
        results.push(("B fused map 16B (f64,u32), AoS24 comp", best, check));
    }

    // ---- C: fused map 16B + SoA comp streams (lambda f64 + index u32 = 12B) ----
    {
        let mut best = f64::MAX;
        let mut check = 0u64;
        for _ in 0..reps {
            let mut map: Vec<(f64, u32)> = vec![(f64::NEG_INFINITY, 0); nraw];
            let t = Instant::now();
            for (ci, (lam, idxs)) in soa.iter().enumerate() {
                let r = reads[ci] as f64;
                let cu = ci as u32;
                for (l, &ix) in lam.iter().zip(idxs.iter()) {
                    let idx = ix as usize;
                    let e = l * r;
                    let slot = unsafe { map.get_unchecked_mut(idx) };
                    if e > slot.0 {
                        *slot = (e, cu);
                    }
                }
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = map.iter().map(|m| m.1 as u64).sum();
        }
        results.push(("C fused map 16B + SoA12 comp", best, check));
    }

    // ---- D: split arrays emax f64 + best_ci u32 (12B/raw), SoA comps -------
    // Two gathers again, but a much smaller second one — separates "fewer
    // bytes" from "fewer cache lines touched per comp".
    {
        let mut best = f64::MAX;
        let mut check = 0u64;
        for _ in 0..reps {
            let mut emax = vec![f64::NEG_INFINITY; nraw];
            let mut best_ci = vec![0u32; nraw];
            let t = Instant::now();
            for (ci, (lam, idxs)) in soa.iter().enumerate() {
                let r = reads[ci] as f64;
                let cu = ci as u32;
                for (l, &ix) in lam.iter().zip(idxs.iter()) {
                    let idx = ix as usize;
                    let e = l * r;
                    if e > unsafe { *emax.get_unchecked(idx) } {
                        unsafe {
                            *emax.get_unchecked_mut(idx) = e;
                            *best_ci.get_unchecked_mut(idx) = cu;
                        }
                    }
                }
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = best_ci.iter().map(|&c| c as u64).sum();
        }
        results.push(("D emax f64 + best_ci u32 (2 gathers), SoA12 comp", best, check));
    }

    // ---- E: control — no gather at all (sequential-only lower bound) -------
    {
        let mut best = f64::MAX;
        for _ in 0..reps {
            let t = Instant::now();
            let mut acc = 0.0f64;
            for (ci, (lam, _)) in soa.iter().enumerate() {
                let r = reads[ci] as f64;
                for l in lam {
                    acc += l * r;
                }
            }
            std::hint::black_box(acc);
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
        }
        results.push(("E control: SoA12 stream, no gather", best, 0));
    }

    // ---- F: raw-major CSR build, lambda-descending, max_reads bound -------
    // Invert the scan: iterate raws ascending over a *flat* CSR candidate index
    // (so the stream is sequential, not the scattered per-raw Vec walk the
    // reconcile pays). Within a raw, candidates are pre-sorted by lambda
    // descending — a static order, since lambda is fixed by b_compare and only
    // reads_j changes. Then reads_j <= max_reads gives the exact bound
    //     lambda * max_reads <= best_e  =>  no later candidate can win
    // and the scan breaks. Exact, not approximate.
    let (csr_lam, csr_ci, csr_off) = {
        let mut per_raw: Vec<Vec<(f64, u32)>> = vec![Vec::new(); nraw];
        for (ci, comps) in clusters.iter().enumerate() {
            for c in comps {
                per_raw[c.index as usize].push((c.lambda, ci as u32));
            }
        }
        let mut lam = Vec::with_capacity(ncomp);
        let mut cis = Vec::with_capacity(ncomp);
        let mut off = Vec::with_capacity(nraw + 1);
        off.push(0u32);
        for v in per_raw.iter_mut() {
            // lambda descending; ties broken by ascending ci so the scan still
            // sees the lowest ci first among equal lambdas.
            v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
            for &(l, c) in v.iter() {
                lam.push(l);
                cis.push(c);
            }
            off.push(lam.len() as u32);
        }
        (lam, cis, off)
    };
    {
        let max_reads = *reads.iter().max().unwrap() as f64;
        let mut best = f64::MAX;
        let mut check = 0u64;
        let mut examined = 0usize;
        for _ in 0..reps {
            let mut emax = vec![f64::NEG_INFINITY; nraw];
            let mut best_ci = vec![0u32; nraw];
            examined = 0;
            let t = Instant::now();
            for raw in 0..nraw {
                let (s, e) = (csr_off[raw] as usize, csr_off[raw + 1] as usize);
                let mut be = f64::NEG_INFINITY;
                let mut bc = u32::MAX;
                let mut k = s;
                while k < e {
                    let l = unsafe { *csr_lam.get_unchecked(k) };
                    // Exact prune: no remaining candidate (lambda is
                    // non-increasing) can beat `be` even at max_reads.
                    if l * max_reads <= be {
                        break;
                    }
                    let ci = unsafe { *csr_ci.get_unchecked(k) };
                    let ev = l * unsafe { *reads.get_unchecked(ci as usize) } as f64;
                    if ev > be || (ev == be && ci < bc) {
                        be = ev;
                        bc = ci;
                    }
                    k += 1;
                }
                examined += k - s;
                emax[raw] = be;
                best_ci[raw] = bc;
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = best_ci.iter().map(|&c| c as u64).sum();
        }
        println!(
            "\n  F prune: examined {:.2}e6 of {:.2}e6 candidates ({:.1}% of the scan, {:.1}x fewer)",
            examined as f64 / 1e6,
            ncomp as f64 / 1e6,
            100.0 * examined as f64 / ncomp as f64,
            ncomp as f64 / examined as f64,
        );
        results.push(("F CSR raw-major, lambda-desc + max_reads bound", best, check));
    }


    // ---- G: two-level bound — hot clusters cluster-major, rest raw-major ----
    // F's bound is lambda * max_reads over ALL clusters, which a power-law
    // abundance distribution makes useless: one huge cluster inflates the
    // bound for every raw. But the huge clusters are *few*. So split the build:
    //   1. the top-T clusters by reads are scanned cluster-major (contiguous,
    //      exactly the access pattern the current build already uses) to seed
    //      emax/best_ci;
    //   2. every other candidate is then bounded by reads_T = the (T+1)-th
    //      largest cluster reads, which is far smaller than max_reads.
    // Selecting top-T costs O(nclusters) per build (hundreds) — negligible.
    // Order-independent: the max uses (e > be) || (e == be && ci < bc), so the
    // two passes compose to the same lowest-ci argmax regardless of visit order.
    for &t_hot in &[8usize, 32, 128] {
        let t_hot = t_hot.min(nclusters);
        let mut order: Vec<u32> = (0..nclusters as u32).collect();
        order.sort_by_key(|&c| std::cmp::Reverse(reads[c as usize]));
        let hot: Vec<u32> = order[..t_hot].to_vec();
        let mut is_hot = vec![false; nclusters];
        for &c in &hot {
            is_hot[c as usize] = true;
        }
        // Bound for everything not scanned in pass 1.
        let reads_t = if t_hot < nclusters { reads[order[t_hot] as usize] as f64 } else { 0.0 };

        let mut best = f64::MAX;
        let mut check = 0u64;
        let mut examined = 0usize;
        for _ in 0..reps {
            let mut emax = vec![f64::NEG_INFINITY; nraw];
            let mut best_ci = vec![u32::MAX; nraw];
            examined = 0;
            let t = Instant::now();
            // pass 1: hot clusters, cluster-major
            for &ci in &hot {
                let (lam, idxs) = &soa[ci as usize];
                let r = reads[ci as usize] as f64;
                for (l, &ix) in lam.iter().zip(idxs.iter()) {
                    let idx = ix as usize;
                    let e = l * r;
                    let em = unsafe { emax.get_unchecked_mut(idx) };
                    if e > *em {
                        *em = e;
                        unsafe { *best_ci.get_unchecked_mut(idx) = ci };
                    } else if e == *em {
                        let b = unsafe { best_ci.get_unchecked_mut(idx) };
                        if ci < *b {
                            *b = ci;
                        }
                    }
                }
                examined += lam.len();
            }
            // pass 2: raw-major CSR over the cold candidates, bounded by reads_t
            for raw in 0..nraw {
                let (s, e) = (csr_off[raw] as usize, csr_off[raw + 1] as usize);
                let mut be = unsafe { *emax.get_unchecked(raw) };
                let mut bc = unsafe { *best_ci.get_unchecked(raw) };
                let mut k = s;
                while k < e {
                    let l = unsafe { *csr_lam.get_unchecked(k) };
                    if l * reads_t <= be {
                        break;
                    }
                    let ci = unsafe { *csr_ci.get_unchecked(k) };
                    if !is_hot[ci as usize] {
                        let ev = l * unsafe { *reads.get_unchecked(ci as usize) } as f64;
                        if ev > be || (ev == be && ci < bc) {
                            be = ev;
                            bc = ci;
                        }
                    }
                    k += 1;
                }
                examined += k - s;
                unsafe {
                    *emax.get_unchecked_mut(raw) = be;
                    *best_ci.get_unchecked_mut(raw) = bc;
                }
            }
            let ns = t.elapsed().as_secs_f64() * 1e9 / ncomp as f64;
            best = best.min(ns);
            check = best_ci.iter().map(|&c| c as u64).sum();
        }
        println!(
            "  G(T={t_hot}) : examined {:.2}e6 of {:.2}e6 ({:.1}% of the scan)  reads_T/max_reads = {:.4}",
            examined as f64 / 1e6,
            ncomp as f64 / 1e6,
            100.0 * examined as f64 / ncomp as f64,
            reads_t / *reads.iter().max().unwrap() as f64,
        );
        let name: &'static str = Box::leak(format!("G two-level bound, T={t_hot} hot clusters").into_boxed_str());
        results.push((name, best, check));
    }

    println!("\n  ns/comp (best of {reps})   vs A     variant");
    let base = results[0].1;
    for (name, ns, check) in &results {
        println!("  {ns:8.2}          {:5.2}x   {name}   [checksum {check}]", base / ns);
    }
}
