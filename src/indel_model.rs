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
//!   lowest quality (novel bases and run-lengthening bases fitted apart) ×
//!   r_ins^(k-1), divided by the match terms λ already gave those bases;
//! - deletion of k bases: the per-base deletion rate on a joint surface of
//!   centre run length (1-4, 5+) and lower flanking quality, × r_del^(k-1);
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
//! reach the abundance test. Rates come from an `--indel-opportunities` table
//! (fitted here by pooling each key with its neighbours) or a fitted-rates
//! table, such as `learn-errors --indel-rates-out` writes or an external
//! `--indel-fit-cmd` returns. The evidence for each choice is on #291.

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

/// Centre run-length classes for deletions: 1, 2, 3, 4, 5 or more.
const RUN_CLASSES: usize = 5;

#[derive(Debug, Clone)]
pub struct IndelModel {
    /// Insertion rate per base by the inserted base's quality: novel bases
    /// (matching neither neighbour), then bases lengthening a run.
    ins: [Vec<f64>; 2],
    /// Deletion rate per centre base by the lower flanking quality, one curve
    /// per centre run-length class.
    del: [Vec<f64>; RUN_CLASSES],
    /// Deletion rate per run-length class when flanking qualities are unknown.
    del_run: [f64; RUN_CLASSES],
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

    fn rate(&self) -> f64 {
        let o: u64 = self.opp.iter().sum();
        if o > 0 {
            self.ev.iter().sum::<u64>() as f64 / o as f64
        } else {
            0.0
        }
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

/// Names of the fitted-rate tables, as written by [`IndelModel::write_rates_tsv`]
/// and expected from an external fitting command.
fn del_name(c: usize) -> String {
    format!("del_fq_L{}", c + 1)
}
const INS_NAMES: [&str; 2] = ["ins_q_novel", "ins_q_run"];

/// A vector over keys 0..=max from (key, value) pairs, gaps filled from the
/// nearest key that has a value.
fn dense(points: &[(usize, f64)]) -> Vec<f64> {
    let max = points.iter().map(|p| p.0).max().unwrap_or(0);
    (0..=max)
        .map(|k| {
            points
                .iter()
                .min_by_key(|p| p.0.abs_diff(k))
                .map_or(0.0, |p| p.1)
        })
        .collect()
}

impl IndelModel {
    /// Load an `--indel-opportunities` counts table (fitted here by pooling) or
    /// a fitted-rates table (`table`, `key`, `rate`; used as given).
    pub fn from_tsv(path: &Path, settings: IndelSettings) -> io::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut lines = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty());
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        if header == ["table", "key", "rate"] {
            let mut tables: std::collections::HashMap<String, Vec<(usize, f64)>> =
                Default::default();
            for (n, line) in lines.enumerate() {
                let f: Vec<&str> = line.split('\t').collect();
                let (Some(k), Some(r)) = (
                    f.get(1).and_then(|v| v.parse::<usize>().ok()),
                    f.get(2).and_then(|v| v.parse::<f64>().ok()),
                ) else {
                    return Err(invalid(format!(
                        "rates line {}: expected table, key, rate",
                        n + 2
                    )));
                };
                if !r.is_finite() || r < 0.0 {
                    return Err(invalid(format!(
                        "rates line {}: rate {r} is not a probability",
                        n + 2
                    )));
                }
                tables.entry(f[0].to_string()).or_default().push((k, r));
            }
            return Self::from_rates(&tables, settings);
        }
        if header.len() < 5 || header[..5] != ["sample", "table", "key", "opportunities", "events"]
        {
            return Err(invalid(
                "expected an --indel-opportunities table or a fitted-rates table".into(),
            ));
        }
        let mut tables: std::collections::HashMap<String, Counts> = Default::default();
        for (n, line) in lines.enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            let parse = |i: usize| -> io::Result<u64> {
                f.get(i)
                    .and_then(|v| v.parse().ok())
                    .ok_or_else(|| invalid(format!("line {}: bad field {i}", n + 2)))
            };
            tables.entry(f[1].to_string()).or_default().add(
                parse(2)? as usize,
                parse(3)?,
                parse(4)?,
            );
        }
        Self::from_counts(&tables, settings)
    }

    /// Fit rates from in-memory tables, as `learn-errors` does each iteration.
    pub fn from_opportunities(
        opp: &crate::indel_events::Opportunities,
        settings: IndelSettings,
    ) -> io::Result<Self> {
        let mut tables: std::collections::HashMap<String, Counts> = Default::default();
        for (name, t) in opp.named_tables() {
            let c = tables.entry(name).or_default();
            for (k, (&o, &e)) in t.opportunities.iter().zip(&t.events).enumerate() {
                c.add(k, o, e);
            }
        }
        Self::from_counts(&tables, settings)
    }

    /// The native fit: per-key pooling of each table (the best non-GAM fit on
    /// #291's held-out comparisons).
    fn from_counts(
        tables: &std::collections::HashMap<String, Counts>,
        settings: IndelSettings,
    ) -> io::Result<Self> {
        let need =
            |t: &str| {
                tables.get(t).filter(|c| c.opp.iter().any(|&o| o > 0)).ok_or_else(|| {
                invalid(format!(
                    "no `{t}` table; write rates with --indel-opportunities from this version"
                ))
            })
            };
        let ins = [
            need(INS_NAMES[0])?.pooled(20),
            need(INS_NAMES[1])?.pooled(20),
        ];
        // A run-length class with no opportunities (no runs that long) borrows
        // the class below it.
        let mut del: [Vec<f64>; RUN_CLASSES] = Default::default();
        let mut del_run = [0.0; RUN_CLASSES];
        for c in 0..RUN_CLASSES {
            match tables
                .get(&del_name(c))
                .filter(|t| t.opp.iter().any(|&o| o > 0))
            {
                Some(t) => {
                    del[c] = t.pooled(20);
                    del_run[c] = t.rate();
                }
                None if c > 0 => {
                    del[c] = del[c - 1].clone();
                    del_run[c] = del_run[c - 1];
                }
                None => return Err(invalid(format!("no `{}` table", del_name(c)))),
            }
        }
        let r = |t: &str| tables.get(t).map_or(1e-6, |c| c.length_factor());
        Ok(IndelModel {
            ins,
            del,
            del_run,
            r_ins: r("ins_len"),
            r_del: r("del_len"),
            settings,
        })
    }

    /// Rates fitted elsewhere, e.g. by an external `--indel-fit-cmd`.
    fn from_rates(
        tables: &std::collections::HashMap<String, Vec<(usize, f64)>>,
        settings: IndelSettings,
    ) -> io::Result<Self> {
        let need = |t: &str| {
            tables
                .get(t)
                .filter(|v| !v.is_empty())
                .ok_or_else(|| invalid(format!("fitted rates: no `{t}` table")))
        };
        let ins = [dense(need(INS_NAMES[0])?), dense(need(INS_NAMES[1])?)];
        let mut del: [Vec<f64>; RUN_CLASSES] = Default::default();
        let mut del_run = [0.0; RUN_CLASSES];
        let by_run = tables.get("del_L");
        for c in 0..RUN_CLASSES {
            del[c] = dense(need(&del_name(c))?);
            del_run[c] = by_run
                .and_then(|v| v.iter().find(|p| p.0 == c + 1).map(|p| p.1))
                .unwrap_or_else(|| {
                    let mut v = del[c].clone();
                    v.sort_by(f64::total_cmp);
                    v[v.len() / 2]
                });
        }
        let scalar = |t: &str| -> io::Result<f64> { Ok(need(t)?[0].1) };
        Ok(IndelModel {
            ins,
            del,
            del_run,
            r_ins: scalar("ins_r")?,
            r_del: scalar("del_r")?,
            settings,
        })
    }

    /// Fit rates with a user-supplied command, as `--errfun-cmd` fits error
    /// rates: `cmd` is whitespace-split into argv, and the paths of the counts
    /// TSV (the `--indel-opportunities` layout, one sample `learned`) and of the
    /// fitted-rates TSV it must write are appended.
    pub fn fit_external(
        opp: &crate::indel_events::Opportunities,
        settings: IndelSettings,
        cmd: &str,
    ) -> io::Result<Self> {
        let argv: Vec<&str> = cmd.split_whitespace().collect();
        if argv.is_empty() {
            return Err(invalid("--indel-fit-cmd is empty".into()));
        }
        let tmp = crate::error_models::TempDir::new("dada2-rs-indelfit")?;
        let (counts, rates) = (tmp.0.join("counts.tsv"), tmp.0.join("rates.tsv"));
        let report = crate::indel_events::IndelReport {
            uniques: Vec::new(),
            opportunities: opp.clone(),
        };
        crate::indel_events::write_opportunities_tsv(&counts, &[], &[("learned", &report)])?;
        let output = std::process::Command::new(argv[0])
            .args(&argv[1..])
            .arg(&counts)
            .arg(&rates)
            .output()
            .map_err(|e| {
                io::Error::other(format!("--indel-fit-cmd: failed to spawn '{cmd}': {e}"))
            })?;
        if !output.status.success() {
            return Err(io::Error::other(format!(
                "--indel-fit-cmd: command '{cmd}' exited with {}\n--- stderr ---\n{}",
                output.status,
                String::from_utf8_lossy(&output.stderr).trim_end(),
            )));
        }
        let m = Self::from_tsv(&rates, settings)
            .map_err(|e| invalid(format!("--indel-fit-cmd output: {e}")))?;
        Ok(m)
    }

    /// Write the fitted rates (`table`, `key`, `rate`), the format `from_tsv`
    /// reads back as given.
    pub fn write_rates_tsv(&self, path: &Path, meta: crate::indel_events::Meta) -> io::Result<()> {
        use std::io::Write;
        let mut w = io::BufWriter::new(std::fs::File::create(path)?);
        for (k, v) in meta {
            writeln!(w, "# {k}: {v}")?;
        }
        writeln!(w, "table\tkey\trate")?;
        for (name, v) in INS_NAMES.iter().zip(&self.ins) {
            for (k, r) in v.iter().enumerate() {
                writeln!(w, "{name}\t{k}\t{r:e}")?;
            }
        }
        for (c, v) in self.del.iter().enumerate() {
            for (k, r) in v.iter().enumerate() {
                writeln!(w, "{}\t{k}\t{r:e}", del_name(c))?;
            }
        }
        for (c, r) in self.del_run.iter().enumerate() {
            writeln!(w, "del_L\t{}\t{r:e}", c + 1)?;
        }
        writeln!(w, "ins_r\t0\t{:e}", self.r_ins)?;
        writeln!(w, "del_r\t0\t{:e}", self.r_del)?;
        w.flush()
    }

    /// Largest relative change in any fitted rate between two models, for
    /// tracking whether an iterated fit is settling.
    pub fn max_rel_change(&self, other: &IndelModel) -> f64 {
        let rel = |a: f64, b: f64| {
            let m = a.abs().max(b.abs());
            if m > 0.0 { (a - b).abs() / m } else { 0.0 }
        };
        let vecs = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .map(|(&x, &y)| rel(x, y))
                .fold(0.0f64, f64::max)
        };
        let mut m = rel(self.r_ins, other.r_ins).max(rel(self.r_del, other.r_del));
        for (a, b) in self.ins.iter().zip(&other.ins) {
            m = m.max(vecs(a, b));
        }
        for (a, b) in self.del.iter().zip(&other.del) {
            m = m.max(vecs(a, b));
        }
        m.max(vecs(&self.del_run, &other.del_run))
    }

    fn at(v: &[f64], k: usize) -> f64 {
        if v.is_empty() {
            0.0
        } else {
            v[k.min(v.len() - 1)]
        }
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
                        let class = usize::from(e.hp_len > 0);
                        let p = Self::at(&self.ins[class], qmin) * self.r_ins.powi(k as i32 - 1);
                        if matched > 0.0 { p / matched } else { p }
                    }
                    Kind::Deletion => {
                        let c = e.hp_len.clamp(1, RUN_CLASSES) - 1;
                        let base = match (e.qual_left, e.qual_right) {
                            (Some(l), Some(r)) => Self::at(&self.del[c], l.min(r) as usize),
                            _ => self.del_run[c],
                        };
                        base * self.r_del.powi(k as i32 - 1)
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
            ins: [vec![1e-3; 41], vec![1e-3; 41]],
            del: [1e-5, 1e-4, 1e-3, 1e-3, 1e-3].map(|r| vec![r; 41]),
            del_run: [1e-5, 1e-4, 1e-3, 1e-3, 1e-3],
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
    fn run_lengthening_insertion_uses_its_own_curve() {
        // A inserted next to a lone A makes a run of 2: the ceiling while it
        // is on, the run-lengthening curve when it is off.
        let mut m = model(ALL_ON);
        m.ins[1] = vec![2e-3; 41];
        assert_eq!(score(&m, "ACGT-ACGTAC", "ACGTAACGTAC"), (0.05, 1));
        m.settings.hp_ceiling = HpCeiling::Off;
        let (f, _) = score(&m, "ACGT-ACGTAC", "ACGTAACGTAC");
        assert!((f - 2e-3 / 0.97).abs() < 1e-12, "{f}");
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
    fn loads_counts_and_rejects_old_tables() {
        let dir = std::env::temp_dir().join(format!("dada2rs_indel_model_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("opp.tsv");
        let mut t =
            String::from("# dada2_rs_command: dada\nsample\ttable\tkey\topportunities\tevents\n");
        // Two samples are summed: 40 novel insertions in 1e6 bases at Q30.
        for s in ["a", "b"] {
            t += &format!("{s}\tins_q_novel\t30\t500000\t20\n{s}\tins_q_run\t30\t500000\t10\n");
            t += &format!("{s}\tdel_fq_L1\t10\t1000\t20\n{s}\tdel_fq_L1\t30\t499000\t5\n");
            t += &format!("{s}\tdel_fq_L3\t30\t1000\t30\n");
            t +=
                &format!("{s}\tins_len\t1\t0\t18\n{s}\tins_len\t2\t0\t2\n{s}\tdel_len\t1\t0\t25\n");
        }
        std::fs::write(&path, &t).unwrap();
        let m = IndelModel::from_tsv(&path, ALL_ON).unwrap();
        assert!((m.ins[0][30] - 4e-5).abs() < 1e-12);
        assert!((m.ins[1][30] - 2e-5).abs() < 1e-12);
        // Run 1: 40/2000 at flanking Q10; Q30's 10 events pool out to Q10,
        // giving 50/1e6, the same as the run's overall rate.
        assert!((m.del[0][10] - 0.02).abs() < 1e-12);
        assert!((m.del[0][30] - 5e-5).abs() < 1e-12);
        assert!((m.del_run[0] - 5e-5).abs() < 1e-12);
        // Run 2 has no opportunities and borrows run 1; run 3 pools to 60/2000.
        assert_eq!(m.del[1], m.del[0]);
        assert!((m.del_run[2] - 0.03).abs() < 1e-12);
        assert!((m.r_ins - (1.0 - 40.0 / 44.0)).abs() < 1e-12);
        std::fs::write(&path, t.replace("ins_q_run", "ins_q_old")).unwrap();
        let e = IndelModel::from_tsv(&path, ALL_ON).unwrap_err().to_string();
        assert!(e.contains("no `ins_q_run` table"), "{e}");
    }

    #[test]
    fn fitted_rates_round_trip_and_are_used_as_given() {
        let dir = std::env::temp_dir().join(format!("dada2rs_indel_rates_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("rates.tsv");
        let a = model(ALL_ON);
        a.write_rates_tsv(&path, &[("dada2_rs_command", "learn-errors".into())])
            .unwrap();
        let b = IndelModel::from_tsv(&path, ALL_ON).unwrap();
        assert_eq!(a.max_rel_change(&b), 0.0);
        // Sparse keys from an external fit fill from the nearest key.
        let t = "table\tkey\trate\nins_q_novel\t10\t1e-3\nins_q_novel\t30\t1e-5\n\
                 ins_q_run\t20\t1e-4\n\
                 del_fq_L1\t20\t1e-5\ndel_fq_L2\t20\t1e-5\ndel_fq_L3\t20\t1e-4\n\
                 del_fq_L4\t20\t1e-3\ndel_fq_L5\t20\t1e-2\n\
                 ins_r\t0\t0.1\ndel_r\t0\t0.2\n";
        std::fs::write(&path, t).unwrap();
        let m = IndelModel::from_tsv(&path, ALL_ON).unwrap();
        let at = |q| IndelModel::at(&m.ins[0], q);
        assert_eq!((at(12), at(25), at(40)), (1e-3, 1e-5, 1e-5));
        assert_eq!((m.del_run[4], m.r_ins, m.r_del), (1e-2, 0.1, 0.2));
        std::fs::write(&path, t.replace("del_r\t0\t0.2", "del_r\t0\t-1")).unwrap();
        assert!(IndelModel::from_tsv(&path, ALL_ON).is_err());
    }

    #[test]
    fn tables_in_memory_and_through_tsv_fit_the_same_model() {
        use crate::indel_events::{
            IndelReport, Opportunities, pair_events, write_opportunities_tsv,
        };
        let (a0, a1) = (enc("ACGTAC-GTACGAAATCC"), enc("A-GTACTGTACG-AATCC"));
        let sub = al2subs(&a0, &a1);
        let c: Vec<u8> = a0.iter().copied().filter(|&b| b != b'-').collect();
        let r: Vec<u8> = a1.iter().copied().filter(|&b| b != b'-').collect();
        let mut opp = Opportunities::default();
        pair_events(
            &c,
            &r,
            Some(&vec![35u8; r.len()]),
            &sub,
            Some((&mut opp, 50)),
        );
        let report = IndelReport {
            uniques: Vec::new(),
            opportunities: opp.clone(),
        };
        let dir =
            std::env::temp_dir().join(format!("dada2rs_indel_model_rt_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("opp.tsv");
        write_opportunities_tsv(&path, &[], &[("s", &report)]).unwrap();
        let a = IndelModel::from_opportunities(&opp, ALL_ON).unwrap();
        let b = IndelModel::from_tsv(&path, ALL_ON).unwrap();
        assert_eq!(a.max_rel_change(&b), 0.0);
        assert_eq!(
            (a.ins[0].len(), a.del[0].len()),
            (b.ins[0].len(), b.del[0].len())
        );
    }
}
