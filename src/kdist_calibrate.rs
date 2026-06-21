//! k-mer-distance vs true-divergence calibration (analysis-only, issue: KDIST_CUTOFF study).
//!
//! DADA2's k-mer screen skips alignment for pairs with k-mer distance above
//! `KDIST_CUTOFF = 0.42`, said to correspond to ~10% nucleotide divergence —
//! a value calibrated on Illumina 16S (k-mer distance traces to ESPRIT, Sun et
//! al. 2009; the 0.42/10% calibration to DADA2, Callahan et al. 2016). The
//! ESPRIT reference implementation is no longer available, so this re-derives
//! the kdist <-> true-divergence relationship empirically on real input data:
//! for a sampled set of unique-sequence pairs it emits both the k-mer distance
//! (`kmer_dist8`, our faithful port of the ESPRIT metric) and the alignment
//! (`align_endsfree`) edit divergence, so the constant can be checked per
//! dataset / platform / pooling regime.
//!
//! POOLING: the pair *population* the screen sees depends on pooling mode —
//! per-sample (one derep JSON) vs full-pool (a pooled seqtab / concatenated
//! uniques). This tool computes pairs over whatever sequence set it is given,
//! so feed it a single-sample derep for per-sample, or a pooled table for
//! pooled; pseudo's screen population is per-sample (priors change the
//! partition, not which pairs are screened).
//!
//! Run (ignored; needs an input JSON):
//!   DADA2RS_KDIST_INPUT=derep.json DADA2RS_KDIST_OUT=kdist.csv \
//!     cargo test --release --bins kdist_calibrate -- --ignored --nocapture
//! Env: DADA2RS_KDIST_K (default 5, = DADA2 default), DADA2RS_KDIST_MAXPAIRS
//! (default 200000; random-subsample above this to bound the O(n^2) blow-up),
//! DADA2RS_KDIST_CUTOFF (default 0.42, for the screened-in flag/summary).

#[cfg(test)]
mod tests {
    use crate::kmers::{assign_kmer8, kmer_dist8};
    use crate::nwalign::align_endsfree;
    use std::io::Write;

    const GAP: u8 = b'-';

    fn encode(seq: &str) -> Vec<u8> {
        seq.bytes()
            .map(|b| match b {
                b'A' | b'a' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 3,
                b'T' | b't' => 4,
                _ => 5, // N etc. — never matches, never a valid k-mer
            })
            .collect()
    }

    /// Internal edit divergence of an ends-free alignment: trim terminal gap
    /// overhang (length difference, not divergence), then count substitutions
    /// and indel columns in the aligned core. Returns (edits, core_len).
    fn aln_divergence(al: &[Vec<u8>; 2]) -> (usize, usize) {
        let (a, b) = (&al[0], &al[1]);
        let n = a.len();
        let mut lo = 0;
        while lo < n && (a[lo] == GAP || b[lo] == GAP) {
            lo += 1;
        }
        let mut hi = n;
        while hi > lo && (a[hi - 1] == GAP || b[hi - 1] == GAP) {
            hi -= 1;
        }
        let mut edits = 0;
        for k in lo..hi {
            if a[k] == GAP || b[k] == GAP {
                edits += 1; // internal indel column
            } else if a[k] != b[k] {
                edits += 1; // substitution
            }
        }
        (edits, hi - lo)
    }

    /// Parse unique sequences from either a derep JSON (`uniques[].sequence` +
    /// `count`) or a seqtab JSON (`sequences[]`, abundance summed over samples).
    fn load_sequences(path: &str) -> Vec<(String, u64)> {
        let txt = std::fs::read_to_string(path).expect("read input JSON");
        let v: serde_json::Value = serde_json::from_str(&txt).expect("parse JSON");
        if let Some(uniques) = v.get("uniques").and_then(|u| u.as_array()) {
            return uniques
                .iter()
                .filter_map(|e| {
                    let s = e.get("sequence")?.as_str()?.to_string();
                    let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
                    Some((s, c))
                })
                .collect();
        }
        if let Some(seqs) = v.get("sequences").and_then(|s| s.as_array()) {
            let counts = v.get("counts").and_then(|c| c.as_array());
            return seqs
                .iter()
                .enumerate()
                .filter_map(|(j, s)| {
                    let seq = s.as_str()?.to_string();
                    let total = counts
                        .map(|rows| {
                            rows.iter()
                                .filter_map(|r| r.as_array()?.get(j)?.as_u64())
                                .sum::<u64>()
                        })
                        .unwrap_or(1);
                    Some((seq, total))
                })
                .collect();
        }
        panic!("input has neither `uniques` (derep) nor `sequences` (seqtab)");
    }

    #[test]
    #[ignore]
    fn kdist_calibrate() {
        let input = match std::env::var("DADA2RS_KDIST_INPUT") {
            Ok(p) => p,
            Err(_) => {
                eprintln!("set DADA2RS_KDIST_INPUT=derep_or_seqtab.json");
                return;
            }
        };
        let k: usize = std::env::var("DADA2RS_KDIST_K")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(5);
        let max_pairs: usize = std::env::var("DADA2RS_KDIST_MAXPAIRS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(200_000);
        let cutoff: f64 = std::env::var("DADA2RS_KDIST_CUTOFF")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.42);
        // A screened-in pair this divergent is too far to be an amplicon-error
        // copy, so its alignment is "leaked" work. Crude (the true ceiling is
        // abundance-dependent); tune via env.
        let leak_pct: f64 = std::env::var("DADA2RS_KDIST_LEAK_PCT")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(5.0);
        let out = std::env::var("DADA2RS_KDIST_OUT").unwrap_or_else(|_| "kdist.csv".into());

        let seqs = load_sequences(&input);
        let n = seqs.len();
        let enc: Vec<Vec<u8>> = seqs.iter().map(|(s, _)| encode(s)).collect();
        let kmers: Vec<Vec<u8>> = enc.iter().map(|e| assign_kmer8(e, k)).collect();
        let total_pairs = n * (n - 1) / 2;
        eprintln!(
            "{n} uniques, {total_pairs} pairs (k={k}, cutoff={cutoff}); \
             {}",
            if total_pairs > max_pairs {
                format!("random-sampling {max_pairs}")
            } else {
                "all pairs".into()
            }
        );

        let mut f = std::io::BufWriter::new(std::fs::File::create(&out).expect("create out"));
        writeln!(f, "kdist,edits,core_len,pct_div,screened_in,ab_i,ab_j").unwrap();

        // simple LCG for reproducible subsampling
        let mut st: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rnd = |m: usize| {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((st >> 33) as usize) % m
        };

        let (mut emitted, mut screened_in, mut leak) = (0u64, 0u64, 0u64);
        let sample = total_pairs > max_pairs;
        let mut do_pair = |i: usize,
                           j: usize,
                           f: &mut std::io::BufWriter<std::fs::File>,
                           emitted: &mut u64,
                           screened_in: &mut u64,
                           leak: &mut u64| {
            let kd = kmer_dist8(&kmers[i], enc[i].len(), &kmers[j], enc[j].len(), k);
            let al = align_endsfree(&enc[i], &enc[j], 5, -4, -8, -1);
            let (edits, core) = aln_divergence(&al);
            let pct = if core > 0 {
                100.0 * edits as f64 / core as f64
            } else {
                0.0
            };
            let scr = kd < cutoff;
            if scr {
                *screened_in += 1;
                if pct > leak_pct {
                    *leak += 1; // screened-in but too divergent to be an error copy
                }
            }
            writeln!(
                f,
                "{kd:.4},{edits},{core},{pct:.3},{},{},{}",
                scr as u8, seqs[i].1, seqs[j].1
            )
            .unwrap();
            *emitted += 1;
        };

        if sample {
            for _ in 0..max_pairs {
                let i = rnd(n);
                let mut j = rnd(n);
                if i == j {
                    j = (j + 1) % n;
                }
                do_pair(
                    i.min(j),
                    i.max(j),
                    &mut f,
                    &mut emitted,
                    &mut screened_in,
                    &mut leak,
                );
            }
        } else {
            for i in 0..n {
                for j in (i + 1)..n {
                    do_pair(i, j, &mut f, &mut emitted, &mut screened_in, &mut leak);
                }
            }
        }
        eprintln!(
            "wrote {emitted} pairs -> {out}\n  screened-in (kdist<{cutoff}): {screened_in} \
             ({:.1}%); of those, {leak} are >{leak_pct}% divergent — too far to be an \
             error copy (the 'leakage' / wasted-alignment cost)",
            100.0 * screened_in as f64 / emitted as f64
        );
    }
}
