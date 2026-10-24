//! Indel events between each unique and its final cluster centre (issue #291).
//!
//! DADA2's error model scores substitutions only; positions inside an alignment
//! gap count as matches, and `hamming` counts substitutions. This diagnostic
//! tabulates what that model ignores, from the same final read-to-centre
//! alignments that give `trans` (R `b_make_transition_by_quality_matrix`), so
//! indel rates and substitution rates share a denominator.
//!
//! Two outputs:
//! - **events**: one row per indel event per unique, all cluster members
//!   (`correct` or not), so indel-only uniques absorbed today can be counted.
//! - **opportunities**: read-weighted event and opportunity counts over
//!   `correct` members only, as `trans` counts, for per-quality and
//!   per-homopolymer-length rates.
//!
//! Gaps outside the first..last aligned centre positions are `terminal`: they
//! come from read-length differences and are free under the ends-free end model.
//! They are reported but excluded from opportunities.

use std::io::{self, Write};
use std::path::Path;

use crate::containers::{B, Sub};
use crate::misc::nt_decode;
use crate::nwalign::GAP_GLYPH;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Insertion,
    Deletion,
}

/// One insertion or deletion of the unique relative to its centre.
#[derive(Debug, Clone, PartialEq)]
pub struct Event {
    pub kind: Kind,
    /// Deletion: first deleted centre position. Insertion: the centre position
    /// the bases follow (0 for a leading terminal insertion).
    pub center_pos: usize,
    /// Inserted (from the unique) or deleted (from the centre) bases, ASCII.
    pub bases: Vec<u8>,
    /// Centre length of the homopolymer run this event lengthens or shortens;
    /// 0 when the bases are mixed or match neither flanking run.
    pub hp_len: usize,
    pub terminal: bool,
    /// Insertion only: min and mean quality of the inserted bases.
    pub qual_min: Option<u8>,
    pub qual_mean: Option<f64>,
    /// Quality of the aligned unique bases either side of the event.
    pub qual_left: Option<u8>,
    pub qual_right: Option<u8>,
}

/// Read-weighted counts for one rate table: `opportunities[k]` and `events[k]`.
#[derive(Debug, Default, Clone)]
pub struct Table {
    pub opportunities: Vec<u64>,
    pub events: Vec<u64>,
}

impl Table {
    fn add(&mut self, key: usize, opp: u64, ev: u64) {
        if self.opportunities.len() <= key {
            self.opportunities.resize(key + 1, 0);
            self.events.resize(key + 1, 0);
        }
        self.opportunities[key] += opp;
        self.events[key] += ev;
    }
}

/// Rate tables, keyed by quality (`*_q`) or centre run length (`*_hp`).
#[derive(Debug, Default, Clone)]
pub struct Opportunities {
    /// Aligned unique bases at quality q; events = substitutions.
    pub sub_q: Table,
    /// Unique bases at quality q; events = inserted bases.
    pub ins_q: Table,
    /// Centre bases in runs of length L; events = deleted bases.
    pub del_hp: Table,
    /// Centre runs of length L; events = insertions lengthening that run.
    pub ins_hp: Table,
}

/// A unique with at least one indel, and its events.
pub struct UniqueEvents {
    pub cluster: usize,
    pub center_len: usize,
    /// Input index of the unique.
    pub unique: u32,
    pub reads: u32,
    pub correct: bool,
    pub p: f64,
    /// λ against the centre and the cluster's reads, the inputs of `p`
    /// (`calc_pA(reads, lambda * cluster_reads)`).
    pub lambda: f64,
    pub cluster_reads: u32,
    pub nsubs: usize,
    pub events: Vec<Event>,
}

pub struct IndelReport {
    pub uniques: Vec<UniqueEvents>,
    pub opportunities: Opportunities,
}

/// Length of the homopolymer run containing each position.
fn run_lengths(seq: &[u8]) -> Vec<usize> {
    let mut out = vec![0; seq.len()];
    let mut start = 0;
    for i in 1..=seq.len() {
        if i == seq.len() || seq[i] != seq[start] {
            out[start..i].fill(i - start);
            start = i;
        }
    }
    out
}

fn same_base(bases: &[u8]) -> Option<u8> {
    let b = *bases.first()?;
    bases.iter().all(|&x| x == b).then_some(b)
}

/// Indel events of `raw` (encoded 1-4) against `center`, from their `Sub`, and
/// this pair's contribution to `opp` at weight `reads` (`None` skips it).
pub fn pair_events(
    center: &[u8],
    raw: &[u8],
    qual: Option<&[u8]>,
    sub: &Sub,
    opp: Option<(&mut Opportunities, u64)>,
) -> Vec<Event> {
    let map = &sub.map;
    let aligned: Vec<usize> = (0..center.len()).filter(|&i| map[i] != GAP_GLYPH).collect();
    let (Some(&first), Some(&last)) = (aligned.first(), aligned.last()) else {
        return Vec::new();
    };
    let runs = run_lengths(center);
    let q = |j: usize| qual.map(|q| q[j]);
    let decode = |s: &[u8]| s.iter().map(|&b| nt_decode(b)).collect::<Vec<u8>>();
    let mut events = Vec::new();

    let insertion = |ins: &[u8], ins_q: Option<&[u8]>, after: Option<usize>, terminal: bool| {
        // An inserted run lengthens a flanking centre run of the same base.
        let flanks = match after {
            Some(a) => vec![a, a + 1],
            None => vec![0],
        };
        let hp_len = same_base(ins).map_or(0, |n| {
            flanks
                .into_iter()
                .find(|&i| i < center.len() && center[i] == n)
                .map_or(0, |i| runs[i])
        });
        Event {
            kind: Kind::Insertion,
            center_pos: after.unwrap_or(0),
            bases: decode(ins),
            hp_len,
            terminal,
            qual_min: ins_q.and_then(|v| v.iter().copied().min()),
            qual_mean: ins_q.map(|v| v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64),
            qual_left: None,
            qual_right: None,
        }
    };
    let deletion = |s: usize, e: usize, terminal: bool| Event {
        kind: Kind::Deletion,
        center_pos: s,
        bases: decode(&center[s..e]),
        hp_len: if same_base(&center[s..e]).is_some() {
            runs[s]
        } else {
            0
        },
        terminal,
        qual_min: None,
        qual_mean: None,
        qual_left: None,
        qual_right: None,
    };
    let slice_q = |a: usize, b: usize| qual.map(|q| &q[a..b]);

    // Leading terminal gaps.
    if first > 0 {
        events.push(deletion(0, first, true));
    }
    let q0 = map[first] as usize;
    if q0 > 0 {
        events.push(insertion(&raw[..q0], slice_q(0, q0), None, true));
    }

    // Internal: walk consecutive aligned centre positions a < b.
    for w in aligned.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (ja, jb) = (map[a] as usize, map[b] as usize);
        if b > a + 1 {
            let mut ev = deletion(a + 1, b, false);
            ev.qual_left = q(ja);
            ev.qual_right = q(jb);
            events.push(ev);
        }
        if jb > ja + 1 {
            let mut ev = insertion(&raw[ja + 1..jb], slice_q(ja + 1, jb), Some(a), false);
            ev.qual_left = q(ja);
            ev.qual_right = q(jb);
            events.push(ev);
        }
    }

    // Trailing terminal gaps.
    if last + 1 < center.len() {
        events.push(deletion(last + 1, center.len(), true));
    }
    let ql = map[last] as usize;
    if ql + 1 < raw.len() {
        events.push(insertion(
            &raw[ql + 1..],
            slice_q(ql + 1, raw.len()),
            Some(last),
            true,
        ));
    }

    if let Some((opp, w)) = opp {
        let qk = |j: usize| qual.map_or(0, |q| q[j] as usize);
        let mut is_sub = vec![false; center.len()];
        for &p in &sub.pos {
            is_sub[p as usize] = true;
        }
        for &i in &aligned {
            opp.sub_q
                .add(qk(map[i] as usize), w, if is_sub[i] { w } else { 0 });
        }
        let mut inserted = vec![false; raw.len()];
        for ev in events
            .iter()
            .filter(|e| e.kind == Kind::Insertion && !e.terminal)
        {
            let j0 = map[ev.center_pos] as usize + 1;
            inserted[j0..j0 + ev.bases.len()].fill(true);
        }
        for (j, &ins) in inserted.iter().enumerate().take(ql + 1).skip(q0) {
            opp.ins_q.add(qk(j), w, if ins { w } else { 0 });
        }
        for i in first..=last {
            let deleted = map[i] == GAP_GLYPH;
            opp.del_hp.add(runs[i], w, if deleted { w } else { 0 });
            if i == 0 || center[i] != center[i - 1] {
                opp.ins_hp.add(runs[i], w, 0);
            }
        }
        for ev in events
            .iter()
            .filter(|e| e.kind == Kind::Insertion && !e.terminal)
        {
            if ev.hp_len > 0 {
                opp.ins_hp.add(ev.hp_len, 0, w);
            }
        }
    }
    events
}

/// Tabulate indel events over the final partition, from `subs` indexed by
/// `raw.index` (the final-subs pass, as for `trans`).
pub fn tabulate(b: &B, subs: &[Option<Sub>]) -> IndelReport {
    let mut uniques = Vec::new();
    let mut opp = Opportunities::default();
    for (ci, bi) in b.clusters.iter().enumerate() {
        let Some(center_idx) = bi.center else {
            continue;
        };
        let center = &b.raws[center_idx];
        for &ri in &bi.raws {
            let raw = &b.raws[ri];
            let Some(sub) = subs.get(raw.index as usize).and_then(|s| s.as_ref()) else {
                continue;
            };
            let weight = raw.correct.then_some((&mut opp, raw.reads as u64));
            let events = pair_events(&center.seq, &raw.seq, raw.qual.as_deref(), sub, weight);
            if !events.is_empty() {
                uniques.push(UniqueEvents {
                    cluster: ci,
                    center_len: center.seq.len(),
                    unique: raw.index,
                    reads: raw.reads,
                    correct: raw.correct,
                    p: raw.p,
                    lambda: raw.comp.lambda,
                    cluster_reads: bi.reads,
                    nsubs: sub.nsubs(),
                    events,
                });
            }
        }
    }
    uniques.sort_by_key(|u| (u.cluster, u.unique));
    IndelReport {
        uniques,
        opportunities: opp,
    }
}

fn create(path: &Path) -> io::Result<io::BufWriter<std::fs::File>> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    Ok(io::BufWriter::new(std::fs::File::create(path)?))
}

fn opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| x.to_string())
}

/// Attributes a unique to the samples it occurs in, with that sample's reads:
/// one entry per sample for `dada`, one per containing sample for `dada-pooled`.
pub type Attribution<'a> = &'a dyn Fn(&UniqueEvents) -> Vec<(String, u32)>;

/// Run identity written as `# key: value` lines above the column header.
pub type Meta<'a> = &'a [(&'static str, String)];

fn write_meta(w: &mut impl Write, meta: Meta) -> io::Result<()> {
    for (k, v) in meta {
        writeln!(w, "# {k}: {v}")?;
    }
    Ok(())
}

/// Write the events TSV for one or more reports; returns the row count.
pub fn write_events_tsv(
    path: &Path,
    meta: Meta,
    parts: &[(&IndelReport, Attribution)],
) -> io::Result<usize> {
    let mut w = create(path)?;
    write_meta(&mut w, meta)?;
    writeln!(
        w,
        "sample\tcluster\tcenter_len\tunique\treads\tcorrect\tp\tlambda\tcluster_reads\tnsubs\tn_internal\t\
         type\tcenter_pos\tlen\tbases\thp_len\tterminal\tqual_min\tqual_mean\tqual_left\tqual_right"
    )?;
    let mut n = 0;
    for (report, samples) in parts {
        for u in &report.uniques {
            let n_internal = u.events.iter().filter(|e| !e.terminal).count();
            for (sample, reads) in samples(u) {
                for e in &u.events {
                    writeln!(
                        w,
                        "{sample}\t{}\t{}\t{}\t{reads}\t{}\t{:e}\t{:e}\t{}\t{}\t{n_internal}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        u.cluster,
                        u.center_len,
                        u.unique,
                        u.correct as u8,
                        u.p,
                        u.lambda,
                        u.cluster_reads,
                        u.nsubs,
                        if e.kind == Kind::Insertion { "I" } else { "D" },
                        e.center_pos,
                        e.bases.len(),
                        String::from_utf8_lossy(&e.bases),
                        e.hp_len,
                        e.terminal as u8,
                        opt(e.qual_min),
                        opt(e.qual_mean.map(|m| format!("{m:.2}"))),
                        opt(e.qual_left),
                        opt(e.qual_right),
                    )?;
                    n += 1;
                }
            }
        }
    }
    w.flush()?;
    Ok(n)
}

/// Write the opportunities TSV (long format: one row per sample, table and key).
pub fn write_opportunities_tsv(
    path: &Path,
    meta: Meta,
    parts: &[(&str, &IndelReport)],
) -> io::Result<usize> {
    let mut w = create(path)?;
    write_meta(&mut w, meta)?;
    writeln!(w, "sample\ttable\tkey\topportunities\tevents")?;
    let mut n = 0;
    for (sample, report) in parts {
        let o = &report.opportunities;
        for (name, t) in [
            ("sub_q", &o.sub_q),
            ("ins_q", &o.ins_q),
            ("del_hp", &o.del_hp),
            ("ins_hp", &o.ins_hp),
        ] {
            for (k, (&opp, &ev)) in t.opportunities.iter().zip(&t.events).enumerate() {
                if opp > 0 || ev > 0 {
                    writeln!(w, "{sample}\t{name}\t{k}\t{opp}\t{ev}")?;
                    n += 1;
                }
            }
        }
    }
    w.flush()?;
    Ok(n)
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

    /// Events for an alignment given as two gapped strings (centre, unique).
    fn run(al0: &str, al1: &str, qual: Option<&[u8]>) -> (Vec<Event>, Opportunities) {
        let (a0, a1) = (enc(al0), enc(al1));
        let sub = al2subs(&a0, &a1);
        let center: Vec<u8> = a0.iter().copied().filter(|&b| b != b'-').collect();
        let raw: Vec<u8> = a1.iter().copied().filter(|&b| b != b'-').collect();
        let mut opp = Opportunities::default();
        let ev = pair_events(&center, &raw, qual, &sub, Some((&mut opp, 10)));
        (ev, opp)
    }

    #[test]
    fn no_indels_no_events() {
        let (ev, opp) = run("ACGTACGT", "ACGAACGT", None);
        assert!(ev.is_empty());
        assert_eq!(opp.sub_q.opportunities[0], 80);
        assert_eq!(opp.sub_q.events[0], 10);
    }

    #[test]
    fn internal_homopolymer_deletion() {
        // Centre run AAAA (len 4) loses one A.
        let (ev, opp) = run("CGAAAATC", "CG-AAATC", Some(&[30, 31, 32, 33, 34, 35, 36]));
        assert_eq!(ev.len(), 1);
        let e = &ev[0];
        assert_eq!(
            (e.kind, e.center_pos, e.hp_len, e.terminal),
            (Kind::Deletion, 2, 4, false)
        );
        assert_eq!(e.bases, b"A");
        assert_eq!((e.qual_left, e.qual_right), (Some(31), Some(32)));
        assert_eq!(opp.del_hp.events[4], 10);
        assert_eq!(opp.del_hp.opportunities[4], 40);
    }

    #[test]
    fn internal_insertion_extends_flanking_run() {
        // Unique has an extra T after the centre's TT run.
        let qual = [20, 21, 22, 23, 9, 25, 26];
        let (ev, opp) = run("CGTT-AC", "CGTTTAC", Some(&qual));
        assert_eq!(ev.len(), 1);
        let e = &ev[0];
        assert_eq!((e.kind, e.center_pos, e.hp_len), (Kind::Insertion, 3, 2));
        assert_eq!(
            (e.qual_min, e.qual_left, e.qual_right),
            (Some(9), Some(23), Some(25))
        );
        assert_eq!(opp.ins_q.events[9], 10);
        assert_eq!(opp.ins_hp.events[2], 10);
    }

    #[test]
    fn novel_and_mixed_insertions_are_not_homopolymer() {
        let (ev, _) = run("CGTA--C", "CGTAGCC", None);
        assert_eq!(ev[0].hp_len, 0); // GC: mixed bases
        let (ev, _) = run("CT-AC", "CTGAC", None);
        assert_eq!(ev[0].hp_len, 0); // G next to T and A
    }

    #[test]
    fn terminal_gaps_flagged_and_excluded() {
        let (ev, opp) = run("--ACGTAC-", "GGACG-ACT", None);
        let kinds: Vec<(Kind, bool)> = ev.iter().map(|e| (e.kind, e.terminal)).collect();
        assert_eq!(
            kinds,
            vec![
                (Kind::Insertion, true),
                (Kind::Deletion, false),
                (Kind::Insertion, true)
            ]
        );
        // Internal unique bases only: ACGAC.
        assert_eq!(opp.ins_q.opportunities[0], 50);
    }

    #[test]
    fn run_lengths_by_position() {
        assert_eq!(run_lengths(&enc("AACGGG")), vec![2, 2, 1, 3, 3, 3]);
    }
}
