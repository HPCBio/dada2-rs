//! Experimental indel-aware λ (issue #291). Off by default; R DADA2 scores
//! substitutions only.
//!
//! DADA2's λ is the probability that a read of the centre comes out as the
//! unique, a product over the unique's positions of error-model terms in which
//! a gap scores as a match. This model multiplies in one probability per indel
//! event between the two, so an indel-only variant is no longer certain to be
//! absorbed:
//!
//! - insertion of k bases: the per-base insertion rate at the inserted bases'
//!   lowest quality × r_ins^(k-1), divided by the match terms λ already gave
//!   those bases;
//! - deletion of k bases: the per-base deletion rate at the centre run length,
//!   scaled by the deletion enrichment at the lower flanking quality, ×
//!   r_del^(k-1);
//! - an indel that only lengthens or shortens a homopolymer run (largest run
//!   involved >= 2): a fixed ceiling probability, set above the rate of
//!   recurring slippage hotspots. A ceiling of 1 never splits on them; `off`
//!   scores them at the measured rates;
//! - an insertion and a deletion of equal length within `pair_window` bases:
//!   the substitutions they stand for, scored from the error matrix;
//! - terminal gaps and indels within `end_zone` bases of the centre's 3' end:
//!   ignored, as the ends-free alignment already treats read ends.
//!
//! Counted events also count toward `hamming`, so an indel-only unique can
//! reach the abundance test. Rates come from a prior run's
//! `--indel-opportunities` table. The evidence for each choice is on #291.

use std::io;
use std::path::Path;

use crate::containers::Sub;
use crate::indel_events::{Event, Kind, pair_events};
use crate::nwalign::GAP_GLYPH;

/// How homopolymer-length indels are scored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HpCeiling {
    /// A fixed probability per event; 1 means they never count.
    Fixed(f64),
    /// The measured rates, like any other indel.
    Off,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndelSettings {
    pub hp_ceiling: HpCeiling,
    pub pair_window: usize,
    pub end_zone: usize,
}

#[derive(Debug, Clone)]
pub struct IndelModel {
    p_ins: Vec<f64>,
    p_del: Vec<f64>,
    /// Deletion enrichment by lower flanking quality (1 = the mean rate).
    del_fq: Vec<f64>,
    r_ins: f64,
    r_del: f64,
    pub settings: IndelSettings,
}

/// Read-weighted counts per key from the opportunities TSV, summed over samples.
#[derive(Default)]
struct Counts {
    opp: Vec<u64>,
    ev: Vec<u64>,
}

impl Counts {
    fn add(&mut self, k: usize, o: u64, e: u64) {
        if self.opp.len() <= k {
            self.opp.resize(k + 1, 0);
            self.ev.resize(k + 1, 0);
        }
        self.opp[k] += o;
        self.ev[k] += e;
    }

    /// Rate per key, widening a window over neighbouring keys until it holds
    /// `min_events` events, so sparse keys borrow from their neighbours.
    fn pooled(&self, min_events: u64) -> Vec<f64> {
        let n = self.opp.len();
        (0..n)
            .map(|k| {
                let mut w = 0;
                loop {
                    let (lo, hi) = (k.saturating_sub(w), (k + w).min(n - 1));
                    let e: u64 = self.ev[lo..=hi].iter().sum();
                    let o: u64 = self.opp[lo..=hi].iter().sum();
                    if e >= min_events || (lo == 0 && hi == n - 1) {
                        break if o > 0 { e as f64 / o as f64 } else { 0.0 };
                    }
                    w += 1;
                }
            })
            .collect()
    }

    /// r of a geometric length model, from events by length.
    fn length_factor(&self) -> f64 {
        let n: u64 = self.ev.iter().sum();
        let tot: u64 = self.ev.iter().enumerate().map(|(k, &e)| k as u64 * e).sum();
        let mean = if n > 0 { tot as f64 / n as f64 } else { 1.0 };
        (1.0 - 1.0 / mean).max(1e-6)
    }
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl IndelModel {
    /// Load rates from an `--indel-opportunities` TSV.
    pub fn from_tsv(path: &Path, settings: IndelSettings) -> io::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut tables: std::collections::HashMap<String, Counts> = Default::default();
        for (n, line) in text.lines().enumerate() {
            if line.starts_with('#') || line.starts_with("sample\t") || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let parse = |i: usize| -> io::Result<u64> {
                f.get(i).and_then(|v| v.parse().ok()).ok_or_else(|| {
                    invalid(format!("{}: line {}: bad field {i}", path.display(), n + 1))
                })
            };
            if f.len() < 5 {
                return Err(invalid(format!(
                    "{}: line {}: expected 5 columns",
                    path.display(),
                    n + 1
                )));
            }
            tables.entry(f[1].to_string()).or_default().add(
                parse(2)? as usize,
                parse(3)?,
                parse(4)?,
            );
        }
        let need = |t: &str| {
            tables.get(t).filter(|c| !c.opp.is_empty()).ok_or_else(|| {
                invalid(format!(
                    "no `{t}` table; write rates with --indel-opportunities from this version"
                ))
            })
        };
        let ins_q = need("ins_q")?;
        let del_hp = need("del_hp")?;
        let del_fq = need("del_fq")?;
        let (ins_len, del_len) = (need("ins_len")?, need("del_len")?);
        let mean_fq = {
            let e: u64 = del_fq.ev.iter().sum();
            let o: u64 = del_fq.opp.iter().sum();
            if o > 0 { e as f64 / o as f64 } else { 0.0 }
        };
        let fq = del_fq
            .pooled(20)
            .into_iter()
            .map(|r| if mean_fq > 0.0 { r / mean_fq } else { 1.0 })
            .collect();
        Ok(IndelModel {
            p_ins: ins_q.pooled(20),
            p_del: del_hp.pooled(20),
            del_fq: fq,
            r_ins: ins_len.length_factor(),
            r_del: del_len.length_factor(),
            settings,
        })
    }

    fn at(v: &[f64], k: usize) -> f64 {
        v[k.min(v.len() - 1)]
    }

    /// λ factor and extra `hamming` for `raw` (encoded 1-4) against `center`.
    /// `err_mat` is the 16 × `ncol` error matrix λ was computed from.
    pub fn score(
        &self,
        center: &[u8],
        raw: &[u8],
        qual: Option<&[u8]>,
        sub: &Sub,
        err_mat: &[f64],
        ncol: usize,
    ) -> (f64, u32) {
        // Cheap exit: no gap on either side.
        let has_gap = sub.map.contains(&GAP_GLYPH)
            || sub
                .map
                .windows(2)
                .any(|w| w[1] != GAP_GLYPH && w[0] != GAP_GLYPH && w[1] > w[0] + 1)
            || center.len() != raw.len();
        if !has_gap {
            return (1.0, 0);
        }
        let s = &self.settings;
        let mut events: Vec<Event> = pair_events(center, raw, qual, sub, None)
            .into_iter()
            .filter(|e| !e.terminal && center.len() - e.center_pos > s.end_zone)
            .collect();
        let q_at = |j: usize| qual.map_or(0, |q| q[j] as usize).min(ncol - 1);
        let err = |t: usize, q: usize| err_mat[t * ncol + q];
        let nt = |b: u8| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
        let mut factor = 1.0f64;
        let mut hamming = 0u32;

        // Equal-length insertion/deletion pairs: the substitutions they stand for.
        if s.pair_window > 0 {
            let mut used = vec![false; events.len()];
            for i in 0..events.len() {
                if used[i] || events[i].kind != Kind::Insertion {
                    continue;
                }
                let best = (0..events.len())
                    .filter(|&d| {
                        !used[d]
                            && events[d].kind == Kind::Deletion
                            && events[d].bases.len() == events[i].bases.len()
                            && events[d].center_pos.abs_diff(events[i].center_pos) <= s.pair_window
                    })
                    .min_by_key(|&d| events[d].center_pos.abs_diff(events[i].center_pos));
                if let Some(d) = best {
                    used[i] = true;
                    used[d] = true;
                    let j0 = sub.map[events[i].center_pos] as usize + 1;
                    for (k, (&b1, &b0)) in events[i].bases.iter().zip(&events[d].bases).enumerate()
                    {
                        let q = q_at(j0 + k);
                        let (n0, n1) = (nt(b0), nt(b1));
                        let m = err(n1 * 4 + n1, q);
                        if m > 0.0 {
                            factor *= (err(n0 * 4 + n1, q) / m).min(1.0);
                        }
                    }
                    hamming += events[i].bases.len() as u32;
                }
            }
            events = events
                .into_iter()
                .zip(used)
                .filter_map(|(e, u)| (!u).then_some(e))
                .collect();
        }

        for e in &events {
            let k = e.bases.len();
            let involved = match e.kind {
                Kind::Deletion => e.hp_len.max(1),
                Kind::Insertion if e.hp_len > 0 => e.hp_len + 1,
                Kind::Insertion => 1,
            };
            let p = match (involved >= 2, s.hp_ceiling) {
                (true, HpCeiling::Fixed(c)) if c >= 1.0 => continue,
                (true, HpCeiling::Fixed(c)) => c,
                _ => match e.kind {
                    Kind::Insertion => {
                        let j0 = sub.map[e.center_pos] as usize + 1;
                        let qmin = (j0..j0 + k).map(q_at).min().unwrap_or(0);
                        // λ scored the inserted bases as matches; replace those terms.
                        let matched: f64 = (j0..j0 + k)
                            .map(|j| {
                                let n = (raw[j].saturating_sub(1) as usize).min(3);
                                err(n * 4 + n, q_at(j))
                            })
                            .product();
                        let p = Self::at(&self.p_ins, qmin) * self.r_ins.powi(k as i32 - 1);
                        if matched > 0.0 { p / matched } else { p }
                    }
                    Kind::Deletion => {
                        let mut p =
                            Self::at(&self.p_del, e.hp_len.max(1)) * self.r_del.powi(k as i32 - 1);
                        if let (Some(l), Some(r)) = (e.qual_left, e.qual_right) {
                            p *= Self::at(&self.del_fq, l.min(r) as usize);
                        }
                        p
                    }
                },
            };
            factor *= p.min(1.0);
            hamming += 1;
        }
        (factor, hamming)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nwalign::al2subs;

    fn enc(s: &str) -> Vec<u8> {
        s.bytes()
            .map(|c| match c {
                b'A' => 1,
                b'C' => 2,
                b'G' => 3,
                b'T' => 4,
                _ => b'-',
            })
            .collect()
    }

    fn model(settings: IndelSettings) -> IndelModel {
        IndelModel {
            p_ins: vec![1e-3; 41],
            p_del: vec![1e-5, 1e-5, 1e-4, 1e-3],
            del_fq: vec![1.0; 41],
            r_ins: 0.05,
            r_del: 0.01,
            settings,
        }
    }

    const ALL_ON: IndelSettings = IndelSettings {
        hp_ceiling: HpCeiling::Fixed(0.05),
        pair_window: 5,
        end_zone: 0,
    };

    /// A flat error matrix: 1% per substitution, 97% match, at every quality.
    fn err() -> Vec<f64> {
        (0..16)
            .flat_map(|t| std::iter::repeat_n(if t / 4 == t % 4 { 0.97 } else { 0.01 }, 41))
            .collect()
    }

    fn score(m: &IndelModel, al0: &str, al1: &str) -> (f64, u32) {
        let (a0, a1) = (enc(al0), enc(al1));
        let sub = al2subs(&a0, &a1);
        let c: Vec<u8> = a0.iter().copied().filter(|&b| b != b'-').collect();
        let r: Vec<u8> = a1.iter().copied().filter(|&b| b != b'-').collect();
        let q = vec![40u8; r.len()];
        m.score(&c, &r, Some(&q), &sub, &err(), 41)
    }

    #[test]
    fn no_indel_no_change() {
        assert_eq!(score(&model(ALL_ON), "ACGTACGTAC", "ACGAACGTAC"), (1.0, 0));
    }

    #[test]
    fn lone_base_deletion_scored_at_run_length_rate() {
        let (f, h) = score(&model(ALL_ON), "ACGTACGTAC", "ACG-ACGTAC");
        assert_eq!(h, 1);
        assert!((f - 1e-5).abs() < 1e-12, "{f}");
    }

    #[test]
    fn homopolymer_change_uses_ceiling_or_never_counts() {
        let (f, h) = score(&model(ALL_ON), "ACGAAATAC", "ACG-AATAC");
        assert_eq!((f, h), (0.05, 1));
        let hard = IndelSettings {
            hp_ceiling: HpCeiling::Fixed(1.0),
            ..ALL_ON
        };
        assert_eq!(score(&model(hard), "ACGAAATAC", "ACG-AATAC"), (1.0, 0));
        let off = IndelSettings {
            hp_ceiling: HpCeiling::Off,
            ..ALL_ON
        };
        let (f, _) = score(&model(off), "ACGAAATAC", "ACG-AATAC");
        assert!((f - 1e-3).abs() < 1e-12, "run of 3 rate: {f}");
    }

    #[test]
    fn insertion_replaces_the_match_term() {
        // Novel base G between T and A.
        let (f, h) = score(&model(ALL_ON), "ACGT-ACGTAC", "ACGTGACGTAC");
        assert_eq!(h, 1);
        assert!((f - 1e-3 / 0.97).abs() < 1e-12, "{f}");
    }

    #[test]
    fn equal_length_pair_scored_as_substitution() {
        // Deleted C and inserted T within 2 nt: one C->T substitution.
        let (f, h) = score(&model(ALL_ON), "ACGTAC-GTAC", "A-GTACTGTAC");
        assert_eq!(h, 1);
        assert!((f - 0.01 / 0.97).abs() < 1e-12, "{f}");
        let no_pairs = IndelSettings {
            pair_window: 0,
            ..ALL_ON
        };
        assert_eq!(score(&model(no_pairs), "ACGTAC-GTAC", "A-GTACTGTAC").1, 2);
    }

    #[test]
    fn end_zone_drops_indels_near_the_three_prime_end() {
        let near = IndelSettings {
            end_zone: 3,
            ..ALL_ON
        };
        assert_eq!(
            score(&model(near), "ACGTACGTCATG", "ACGTACGTC-TG"),
            (1.0, 0)
        );
        assert_eq!(score(&model(ALL_ON), "ACGTACGTCATG", "ACGTACGTC-TG").1, 1);
    }

    #[test]
    fn longer_indels_pay_the_length_factor() {
        let (f, _) = score(&model(ALL_ON), "ACGTACGTCATG", "AC---CGTCATG");
        assert!((f - 1e-5 * 0.01f64.powi(2)).abs() < 1e-18, "{f}");
    }

    #[test]
    fn loads_rates_and_rejects_old_tables() {
        let dir = std::env::temp_dir().join(format!("dada2rs_indel_model_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("opp.tsv");
        let mut t =
            String::from("# dada2_rs_command: dada\nsample\ttable\tkey\topportunities\tevents\n");
        // Two samples are summed: 40 insertions in 1e6 bases at Q30.
        for s in ["a", "b"] {
            t += &format!("{s}\tins_q\t30\t500000\t20\n{s}\tdel_hp\t1\t500000\t25\n");
            t += &format!("{s}\tdel_fq\t10\t1000\t20\n{s}\tdel_fq\t30\t499000\t5\n");
            t +=
                &format!("{s}\tins_len\t1\t0\t18\n{s}\tins_len\t2\t0\t2\n{s}\tdel_len\t1\t0\t25\n");
        }
        std::fs::write(&path, &t).unwrap();
        let m = IndelModel::from_tsv(&path, ALL_ON).unwrap();
        assert!((m.p_ins[30] - 4e-5).abs() < 1e-12);
        assert!((m.p_del[1] - 5e-5).abs() < 1e-12);
        // Deletions are 40/2000 at flanking Q10 against 50/1e6 overall.
        assert!((m.del_fq[10] - (40.0 / 2000.0) / (50.0 / 1e6)).abs() < 1e-6);
        assert!((m.r_ins - (1.0 - 40.0 / 44.0)).abs() < 1e-12);
        std::fs::write(&path, t.replace("del_fq", "del_old")).unwrap();
        let e = IndelModel::from_tsv(&path, ALL_ON).unwrap_err().to_string();
        assert!(e.contains("no `del_fq` table"), "{e}");
    }
}
