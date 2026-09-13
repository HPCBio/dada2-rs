// Gate for the index-build rewrite: current HashMap<u64, Vec<_>> vs a CSR
// (count / prefix-sum / fill) build, at the real pool shapes.
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

#[derive(Default)]
struct Id(u64);
impl Hasher for Id {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, _: &[u8]) { unreachable!() }
    fn write_u64(&mut self, v: u64) { self.0 = v; }
}
type IdMap<V> = HashMap<u64, V, BuildHasherDefault<Id>>;

fn splitmix(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// One sketch per raw: `epr` distinct keys drawn from a `distinct`-key universe.
fn make(nraw: usize, epr: usize, distinct: u64, seed: u64) -> Vec<Vec<(u64, u8)>> {
    let mut st = seed;
    let keys: Vec<u64> = (0..distinct).map(|_| splitmix(&mut st)).collect();
    (0..nraw)
        .map(|_| {
            let mut v: Vec<(u64, u8)> = (0..epr)
                .map(|_| (keys[(splitmix(&mut st) % distinct) as usize], 1u8))
                .collect();
            v.sort_unstable_by_key(|e| e.0);
            v.dedup_by_key(|e| e.0);
            v
        })
        .collect()
}

fn build_current(sk: &[Vec<(u64, u8)>]) -> usize {
    let inc: usize = sk.iter().map(|s| s.len()).sum();
    let mut p: IdMap<Vec<(u32, u8)>> = IdMap::with_capacity_and_hasher(inc / 8 + 16, Default::default());
    for (i, s) in sk.iter().enumerate() {
        for &(h, c) in s {
            p.entry(h).or_default().push((i as u32, c));
        }
    }
    p.values().map(|v| v.len()).sum()
}

fn build_csr(sk: &[Vec<(u64, u8)>]) -> usize {
    let inc: usize = sk.iter().map(|s| s.len()).sum();
    // pass 1: count per key
    let mut slot: IdMap<u32> = IdMap::with_capacity_and_hasher(inc / 8 + 16, Default::default());
    let mut counts: Vec<u32> = Vec::new();
    for s in sk {
        for &(h, _) in s {
            match slot.get(&h) {
                Some(&k) => counts[k as usize] += 1,
                None => { slot.insert(h, counts.len() as u32); counts.push(1); }
            }
        }
    }
    // prefix sum
    let mut off: Vec<u32> = Vec::with_capacity(counts.len() + 1);
    let mut acc = 0u32;
    for &c in &counts { off.push(acc); acc += c; }
    off.push(acc);
    // pass 2: fill
    let mut cursor = off.clone();
    let mut flat: Vec<(u32, u8)> = vec![(0, 0); acc as usize];
    for (i, s) in sk.iter().enumerate() {
        for &(h, c) in s {
            let k = slot[&h] as usize;
            flat[cursor[k] as usize] = (i as u32, c);
            cursor[k] += 1;
        }
    }
    flat.len()
}

fn main() {
    for (name, nraw, epr, distinct) in [
        ("ITS2  k8 w5", 825_214usize, 74usize, 48_497u64),
        ("ITS2  k6 w5", 825_214, 73, 3_211),
        ("ITS2  k5 w1", 412_607, 199, 1_024),
        ("PacBio k8 w5", 273_636, 472, 40_177),
    ] {
        let sk = make(nraw, epr, distinct, 0x1234);
        let inc: usize = sk.iter().map(|s| s.len()).sum();
        let t = Instant::now(); let a = build_current(&sk); let cur = t.elapsed();
        let t = Instant::now(); let b = build_csr(&sk); let csr = t.elapsed();
        assert_eq!(a, b);
        println!(
            "{name:<13} entries {inc:>11}  posting {:>7}  current {:>7.2}s ({:>6.1} ns/e)  csr {:>6.2}s ({:>5.1} ns/e)  speedup {:>5.1}x",
            inc as u64 / distinct,
            cur.as_secs_f64(), cur.as_nanos() as f64 / inc as f64,
            csr.as_secs_f64(), csr.as_nanos() as f64 / inc as f64,
            cur.as_secs_f64() / csr.as_secs_f64()
        );
    }
}
