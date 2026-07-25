//! `reference-eval` subcommand: NW alignment-based comparison of a query ASV
//! set against a reference (truth) set.
//!
//! v1 scope (issue #91): a single query set vs a single reference. Classifies
//! every query ASV (TP / near / FP) by ends-free global alignment, reports
//! reference recovery (FN) and precision/recall, and joins the classification
//! to the per-ASV p-value already carried in a `dada`/`dada-pooled` JSON. Pure
//! diagnostic — it only *reads* existing fields, never re-infers anything.
//!
//! Deliberately NOT in v1 (kept as forward-compatible hooks, see #91):
//! - Multi-sample / control inputs and the TN question — needs a negative
//!   frame (blanks / closed community); that is a future multi-sample mode.
//!   Hook: every output row is keyed by the ASV *sequence* so independent
//!   single-sample runs are joinable downstream.
//! - Abundance / rRNA-copy-number detection-limit modeling — done downstream.
//!   Hook: reference-FASTA header annotations pass through to the per-ref table
//!   verbatim, and observed abundance per recovered reference is emitted.

use std::io::{self, Write};
use std::path::PathBuf;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::evaluate::eval_pair;
use crate::kmers::{assign_kmer, kmer_dist};
use crate::misc::{Tagged, nt_encode, read_fasta_records, read_tagged_json};
use crate::nwalign::{AlignBuffers, align_endsfree_with_buf};

/// Tuning + I/O parameters for a `reference-eval` run.
pub struct Params {
    pub asvs: PathBuf,
    pub reference: PathBuf,
    /// Edit distance (mismatches + internal indels) at or below which a query
    /// is a true positive. `0` = exact (alignment-internal) match.
    pub max_diffs: u32,
    /// Upper edit-distance bound for the report-only "near" bucket
    /// (`max_diffs < edit <= near_diffs`). Must be `>= max_diffs`.
    pub near_diffs: u32,
    /// Optional k-mer prefilter cutoff. `None` = screen off (recommended for
    /// small reference sets: dropping a true reference misclassifies the ASV
    /// as an FP, so the screen must be permissive when used at all).
    pub kdist_screen: Option<f64>,
    pub kmer_size: usize,
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_p: i32,
    pub band: i32,
    /// Bonferroni divisor (from a pooled trace's `nraw`). When set, a `Prior`
    /// ASV's `birth_pval` (= uncorrected `p_p`) is reported on the abundance
    /// scale as `p_a = birth_pval * nraw`; an `Abundance` ASV's `birth_pval`
    /// is already `p_a`. `None` leaves `p_a` blank.
    pub nraw: Option<f64>,
    pub per_asv: Option<PathBuf>,
    pub per_ref: Option<PathBuf>,
    pub out: Option<PathBuf>,
    pub threads: usize,
    pub compact: bool,
}

/// A query ASV to classify.
struct Query {
    id: String,
    seq_ascii: Vec<u8>,
    enc: Vec<u8>,
    abundance: Option<u32>,
    birth_type: Option<String>,
    birth_pval: Option<f64>,
}

/// A reference (truth) sequence.
struct Reference {
    id: String,
    /// Full FASTA header after `>` (annotations preserved verbatim for the
    /// per-ref metadata pass-through).
    annotation: String,
    enc: Vec<u8>,
    kvec: Vec<u16>,
    len: usize,
}

/// dada / dada-pooled output, minimally deserialized (unknown fields ignored).
#[derive(Deserialize)]
struct DadaDoc {
    asvs: Vec<DadaAsv>,
}
#[derive(Deserialize)]
struct DadaAsv {
    sequence: String,
    #[serde(default)]
    abundance: Option<u32>,
    #[serde(default)]
    birth_type: Option<String>,
    #[serde(default)]
    birth_pval: Option<f64>,
}

fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| nt_encode(b)).collect()
}

/// Load the query set, auto-detecting FASTA (`>`) vs a tagged dada JSON (`{`).
fn load_queries(path: &std::path::Path) -> io::Result<Vec<Query>> {
    let bytes = crate::misc::read_all_maybe_gz(path)?;
    let first = bytes.iter().find(|b| !b.is_ascii_whitespace()).copied();
    if first == Some(b'{') {
        let doc: DadaDoc = read_tagged_json(path, &["dada", "dada-pooled", "dada-pooled-record"])?;
        Ok(doc
            .asvs
            .into_iter()
            .enumerate()
            .map(|(i, a)| {
                let seq_ascii = a.sequence.into_bytes();
                let enc = encode(&seq_ascii);
                Query {
                    id: format!("asv{}", i + 1),
                    seq_ascii,
                    enc,
                    abundance: a.abundance,
                    birth_type: a.birth_type,
                    birth_pval: a.birth_pval,
                }
            })
            .collect())
    } else {
        Ok(read_fasta_records(path)?
            .into_iter()
            .map(|(header, seq_ascii)| {
                let id = header.split_whitespace().next().unwrap_or("").to_string();
                let enc = encode(&seq_ascii);
                Query {
                    id,
                    seq_ascii,
                    enc,
                    abundance: None,
                    birth_type: None,
                    birth_pval: None,
                }
            })
            .collect())
    }
}

fn load_references(path: &std::path::Path, k: usize) -> io::Result<Vec<Reference>> {
    let recs = read_fasta_records(path)?;
    Ok(recs
        .into_iter()
        .map(|(annotation, seq_ascii)| {
            let id = annotation
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let enc = encode(&seq_ascii);
            let len = enc.len();
            let kvec = if len > k {
                assign_kmer(&enc, k)
            } else {
                Vec::new()
            };
            Reference {
                id,
                annotation,
                enc,
                kvec,
                len,
            }
        })
        .collect())
}

/// Classification bucket for a query ASV.
#[derive(Clone, Copy, PartialEq)]
enum Class {
    Tp,
    Near,
    Fp,
}

/// Per-query alignment result.
struct QueryResult {
    idx: usize,
    best_ref: Option<usize>,
    matches: u32,
    mism: u32,
    indel: u32,
    edit: u32,
    class: Class,
}

/// Align one query against the (optionally screened) reference set and pick
/// the best match by minimum edit distance (tie → higher identity → lower ref
/// index). Uses `buf` to avoid per-call allocation.
fn classify_one(q: &Query, refs: &[Reference], p: &Params, buf: &mut AlignBuffers) -> QueryResult {
    let qkvec = if q.enc.len() > p.kmer_size {
        Some(assign_kmer(&q.enc, p.kmer_size))
    } else {
        None
    };

    let mut best: Option<(u32, u32, u32, u32, usize)> = None; // (edit, matches, mism, indel, ref_idx)
    for (ri, r) in refs.iter().enumerate() {
        // Optional permissive k-mer prefilter.
        if let Some(cut) = p.kdist_screen
            && let (Some(qk), false) = (&qkvec, r.kvec.is_empty())
        {
            let d = kmer_dist(qk, q.enc.len(), &r.kvec, r.len, p.kmer_size);
            if d > cut {
                continue;
            }
        }
        if q.enc.is_empty() || r.enc.is_empty() {
            continue;
        }
        align_endsfree_with_buf(
            &r.enc,
            &q.enc,
            p.match_score,
            p.mismatch,
            p.gap_p,
            p.band,
            buf,
        );
        let (al0, al1) = buf.alignment();
        let (m, mism, indel) = match eval_pair(al0, al1) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let edit = mism + indel;
        let better = match best {
            None => true,
            Some((be, bm, _, _, _)) => edit < be || (edit == be && m > bm),
        };
        if better {
            best = Some((edit, m, mism, indel, ri));
        }
    }

    // `idx` is filled in by the caller (map_init) after return.
    match best {
        None => QueryResult {
            idx: 0,
            best_ref: None,
            matches: 0,
            mism: 0,
            indel: 0,
            edit: u32::MAX,
            class: Class::Fp,
        },
        Some((edit, m, mism, indel, ri)) => {
            let class = if edit <= p.max_diffs {
                Class::Tp
            } else if edit <= p.near_diffs {
                Class::Near
            } else {
                Class::Fp
            };
            QueryResult {
                idx: 0,
                best_ref: Some(ri),
                matches: m,
                mism,
                indel,
                edit,
                class,
            }
        }
    }
}

fn identity(m: u32, mism: u32, indel: u32) -> f64 {
    let denom = (m + mism + indel) as f64;
    if denom == 0.0 { 0.0 } else { m as f64 / denom }
}

/// Abundance-scale p-value for the p-value join. `Prior` births carry the
/// uncorrected `p_p`; multiply by `nraw` (if supplied) to compare against
/// omega_a. `Abundance` births already carry `p_a`.
fn p_a(bt: &Option<String>, bp: Option<f64>, nraw: Option<f64>) -> Option<f64> {
    let bp = bp?;
    match bt.as_deref() {
        Some("Prior") => nraw.map(|n| bp * n),
        _ => Some(bp),
    }
}

#[derive(Serialize)]
struct Summary {
    query_input: String,
    reference_input: String,
    n_asvs: usize,
    n_ref: usize,
    tp: usize,
    tp_exact: usize,
    near: usize,
    fp: usize,
    recovered_refs: usize,
    near_recovered_refs: usize,
    fn_refs: usize,
    /// TP fraction of query ASVs (tp / n_asvs). Precision in the open-reference
    /// sense; TN is undefined here (no finite negative set) so specificity /
    /// accuracy are intentionally omitted — see #91.
    precision: f64,
    /// Fraction of references recovered by at least one TP ASV.
    recall: f64,
    /// Edit-distance histogram over best matches: index 0..=3 are exact edits
    /// 0/1/2/3; index 4 is `> 3` (includes no-match). Chosen to make the
    /// exact-vs-near boundary a data decision.
    edit_hist: [usize; 5],
    max_diffs: u32,
    near_diffs: u32,
    band: i32,
    kdist_screen: Option<f64>,
}

pub fn run(p: &Params) -> io::Result<()> {
    if p.near_diffs < p.max_diffs {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--near-diffs ({}) must be >= --max-diffs ({})",
                p.near_diffs, p.max_diffs
            ),
        ));
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads.max(1))
        .build()
        .map_err(io::Error::other)?;

    let mut queries = load_queries(&p.asvs)?;
    for (i, q) in queries.iter_mut().enumerate() {
        q.id = if q.id.is_empty() {
            format!("q{}", i + 1)
        } else {
            std::mem::take(&mut q.id)
        };
    }
    let refs = load_references(&p.reference, p.kmer_size)?;
    if refs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "reference set is empty",
        ));
    }

    // Parallel classification; each query indexes back into `queries` by order.
    let mut results: Vec<QueryResult> = pool.install(|| {
        queries
            .par_iter()
            .enumerate()
            .map_init(AlignBuffers::new, |buf, (i, q)| {
                let mut r = classify_one(q, &refs, p, buf);
                r.idx = i;
                r
            })
            .collect()
    });
    results.sort_by_key(|r| r.idx);

    // Reference recovery: a ref is recovered if some TP query's best match is
    // it (claimant definition; #91). Track per-ref TP count + observed abundance.
    let n_ref = refs.len();
    let mut recovered = vec![false; n_ref];
    let mut near_recovered = vec![false; n_ref];
    let mut ref_tp_count = vec![0usize; n_ref];
    let mut ref_obs_abund = vec![0u64; n_ref];
    let mut ref_best_edit = vec![u32::MAX; n_ref];
    let mut ref_best_query: Vec<Option<usize>> = vec![None; n_ref];

    let mut edit_hist = [0usize; 5];
    let (mut tp, mut tp_exact, mut near, mut fp) = (0, 0, 0, 0);

    for r in &results {
        let h = (r.edit.min(4)) as usize;
        edit_hist[h] += 1;
        match r.class {
            Class::Tp => {
                tp += 1;
                if r.edit == 0 {
                    tp_exact += 1;
                }
            }
            Class::Near => near += 1,
            Class::Fp => fp += 1,
        }
        if let Some(ri) = r.best_ref {
            if r.class == Class::Tp {
                recovered[ri] = true;
                ref_tp_count[ri] += 1;
                ref_obs_abund[ri] += queries[r.idx].abundance.unwrap_or(0) as u64;
            }
            if r.class == Class::Tp || r.class == Class::Near {
                near_recovered[ri] = true;
            }
            if r.edit < ref_best_edit[ri] {
                ref_best_edit[ri] = r.edit;
                ref_best_query[ri] = Some(r.idx);
            }
        }
    }

    // Reference-centric recovery (correct recall). The claimant pass above
    // marks refs that were some query's *best* match, but with near-identical
    // references (e.g. intragenomic 16S copies) an ASV exact to one allele can
    // be assigned to a sibling on a tie, leaving the true allele unclaimed. A
    // targeted rescue pass over still-unrecovered refs restores the intended
    // meaning: a reference is recovered if ANY query aligns within max_diffs.
    // Cheap — only the unrecovered refs, each scanned against all queries.
    let unrec: Vec<usize> = (0..n_ref).filter(|&ri| !recovered[ri]).collect();
    if !unrec.is_empty() {
        let rescued: Vec<(usize, u32, usize)> = pool.install(|| {
            unrec
                .par_iter()
                .filter_map(|&ri| {
                    let r = &refs[ri];
                    if r.enc.is_empty() {
                        return None;
                    }
                    let mut buf = AlignBuffers::new();
                    let mut best: Option<(u32, usize)> = None; // (edit, query_idx)
                    for (qi, q) in queries.iter().enumerate() {
                        if q.enc.is_empty() {
                            continue;
                        }
                        align_endsfree_with_buf(
                            &r.enc,
                            &q.enc,
                            p.match_score,
                            p.mismatch,
                            p.gap_p,
                            p.band,
                            &mut buf,
                        );
                        let (al0, al1) = buf.alignment();
                        if let Ok((_, mism, indel)) = eval_pair(al0, al1) {
                            let edit = mism + indel;
                            if edit <= p.max_diffs && best.map(|(be, _)| edit < be).unwrap_or(true)
                            {
                                best = Some((edit, qi));
                                if edit == 0 {
                                    break;
                                }
                            }
                        }
                    }
                    best.map(|(edit, qi)| (ri, edit, qi))
                })
                .collect()
        });
        for (ri, edit, qi) in rescued {
            recovered[ri] = true;
            near_recovered[ri] = true;
            if edit < ref_best_edit[ri] {
                ref_best_edit[ri] = edit;
                ref_best_query[ri] = Some(qi);
            }
        }
    }

    let recovered_refs = recovered.iter().filter(|&&b| b).count();
    let near_recovered_refs = near_recovered.iter().filter(|&&b| b).count();
    let fn_refs = n_ref - recovered_refs;
    let n_asvs = queries.len();

    let summary = Summary {
        query_input: p.asvs.display().to_string(),
        reference_input: p.reference.display().to_string(),
        n_asvs,
        n_ref,
        tp,
        tp_exact,
        near,
        fp,
        recovered_refs,
        near_recovered_refs,
        fn_refs,
        precision: if n_asvs > 0 {
            tp as f64 / n_asvs as f64
        } else {
            0.0
        },
        recall: if n_ref > 0 {
            recovered_refs as f64 / n_ref as f64
        } else {
            0.0
        },
        edit_hist,
        max_diffs: p.max_diffs,
        near_diffs: p.near_diffs,
        band: p.band,
        kdist_screen: p.kdist_screen,
    };

    // ---- per-ASV table (sequence-keyed) ----
    if let Some(path) = &p.per_asv {
        let mut w = io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "query_id\tseq_len\tabundance\tbirth_type\tbirth_pval\tp_a\tbest_ref\tedit\tmism\tindel\tidentity\tclass\tsequence"
        )?;
        for r in &results {
            let q = &queries[r.idx];
            let best_ref = r.best_ref.map(|ri| refs[ri].id.as_str()).unwrap_or("-");
            let pa = p_a(&q.birth_type, q.birth_pval, p.nraw);
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{}",
                q.id,
                q.seq_ascii.len(),
                q.abundance.map(|a| a.to_string()).unwrap_or_default(),
                q.birth_type.as_deref().unwrap_or(""),
                q.birth_pval.map(|v| format!("{v:.6e}")).unwrap_or_default(),
                pa.map(|v| format!("{v:.6e}")).unwrap_or_default(),
                best_ref,
                if r.edit == u32::MAX {
                    String::new()
                } else {
                    r.edit.to_string()
                },
                r.mism,
                r.indel,
                identity(r.matches, r.mism, r.indel),
                class_str(r.class),
                String::from_utf8_lossy(&q.seq_ascii),
            )?;
        }
        w.flush()?;
    }

    // ---- per-reference table (metadata pass-through) ----
    if let Some(path) = &p.per_ref {
        let mut w = io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            w,
            "ref_id\trecovered\tbest_query_edit\tn_tp_asvs\tobs_abundance\tbest_query_id\tannotation"
        )?;
        for (ri, r) in refs.iter().enumerate() {
            let best_edit = if ref_best_edit[ri] == u32::MAX {
                String::new()
            } else {
                ref_best_edit[ri].to_string()
            };
            let best_q = ref_best_query[ri]
                .map(|qi| queries[qi].id.as_str())
                .unwrap_or("-");
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.id,
                recovered[ri],
                best_edit,
                ref_tp_count[ri],
                ref_obs_abund[ri],
                best_q,
                r.annotation,
            )?;
        }
        w.flush()?;
    }

    // ---- summary JSON + stderr digest ----
    let tagged = Tagged::new("reference-eval", &summary);
    let json = if p.compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;
    match &p.out {
        Some(path) => std::fs::write(path, json)?,
        None => println!("{json}"),
    }

    eprintln!(
        "[reference-eval] {} ASVs vs {} refs: TP={} (exact {}), near={}, FP={} | recovered {}/{} refs (recall {:.1}%), FN={}",
        n_asvs,
        n_ref,
        tp,
        tp_exact,
        near,
        fp,
        recovered_refs,
        n_ref,
        summary.recall * 100.0,
        fn_refs,
    );
    eprintln!(
        "[reference-eval] edit-distance histogram (best match): 0={} 1={} 2={} 3={} >3={}",
        edit_hist[0], edit_hist[1], edit_hist[2], edit_hist[3], edit_hist[4]
    );

    Ok(())
}

fn class_str(c: Class) -> &'static str {
    match c {
        Class::Tp => "TP",
        Class::Near => "near",
        Class::Fp => "FP",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(max_diffs: u32, near_diffs: u32) -> Params {
        Params {
            asvs: PathBuf::new(),
            reference: PathBuf::new(),
            max_diffs,
            near_diffs,
            kdist_screen: None,
            kmer_size: 5,
            match_score: 5,
            mismatch: -4,
            gap_p: -8,
            band: 16,
            nraw: None,
            per_asv: None,
            per_ref: None,
            out: None,
            threads: 1,
            compact: false,
        }
    }

    fn q(seq: &str) -> Query {
        let a = seq.as_bytes().to_vec();
        Query {
            id: "q".into(),
            enc: encode(&a),
            seq_ascii: a,
            abundance: None,
            birth_type: None,
            birth_pval: None,
        }
    }

    fn r(seq: &str, k: usize) -> Reference {
        let a = seq.as_bytes().to_vec();
        let enc = encode(&a);
        let len = enc.len();
        let kvec = if len > k {
            assign_kmer(&enc, k)
        } else {
            Vec::new()
        };
        Reference {
            id: "r".into(),
            annotation: "r".into(),
            enc,
            kvec,
            len,
        }
    }

    const REF: &str = "ACGTACGTACGTACGTACGTACGTACGTACGT";

    #[test]
    fn exact_match_is_tp_edit0() {
        let p = params(0, 3);
        let refs = vec![r(REF, p.kmer_size)];
        let mut buf = AlignBuffers::new();
        let res = classify_one(&q(REF), &refs, &p, &mut buf);
        assert_eq!(res.edit, 0);
        assert!(matches!(res.class, Class::Tp));
        assert_eq!(res.best_ref, Some(0));
    }

    #[test]
    fn one_substitution_is_near_not_tp_at_max0() {
        // flip a single internal base
        let mut m: Vec<u8> = REF.bytes().collect();
        m[10] = if m[10] == b'A' { b'C' } else { b'A' };
        let mut_seq = String::from_utf8(m).unwrap();
        let p = params(0, 3);
        let refs = vec![r(REF, p.kmer_size)];
        let mut buf = AlignBuffers::new();
        let res = classify_one(&q(&mut_seq), &refs, &p, &mut buf);
        assert_eq!(res.mism, 1);
        assert_eq!(res.edit, 1);
        assert!(matches!(res.class, Class::Near));
    }

    #[test]
    fn one_substitution_is_tp_when_max_diffs_1() {
        let mut m: Vec<u8> = REF.bytes().collect();
        m[10] = if m[10] == b'A' { b'C' } else { b'A' };
        let mut_seq = String::from_utf8(m).unwrap();
        let p = params(1, 3);
        let refs = vec![r(REF, p.kmer_size)];
        let mut buf = AlignBuffers::new();
        let res = classify_one(&q(&mut_seq), &refs, &p, &mut buf);
        assert!(matches!(res.class, Class::Tp));
    }

    #[test]
    fn divergent_is_fp() {
        let p = params(0, 3);
        let refs = vec![r(REF, p.kmer_size)];
        let mut buf = AlignBuffers::new();
        // unrelated sequence
        let res = classify_one(&q("TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA"), &refs, &p, &mut buf);
        assert!(res.edit > p.near_diffs);
        assert!(matches!(res.class, Class::Fp));
    }

    #[test]
    fn endsfree_truncation_of_reference_is_tp() {
        // query is an internal substring of the reference: end gaps are free,
        // so the internal edit distance is 0 -> exact TP.
        let p = params(0, 3);
        let refs = vec![r(REF, p.kmer_size)];
        let mut buf = AlignBuffers::new();
        let sub = &REF[4..28]; // drop 4 bases each end
        let res = classify_one(&q(sub), &refs, &p, &mut buf);
        assert_eq!(res.edit, 0);
        assert!(matches!(res.class, Class::Tp));
    }

    #[test]
    fn best_of_multiple_refs_is_min_edit() {
        let p = params(1, 3);
        let mut m: Vec<u8> = REF.bytes().collect();
        m[5] = if m[5] == b'A' { b'G' } else { b'A' };
        let one_off = String::from_utf8(m).unwrap();
        // refs: an unrelated one and the exact one; query == REF should pick exact.
        let refs = vec![
            r("TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA", p.kmer_size),
            r(&one_off, p.kmer_size),
            r(REF, p.kmer_size),
        ];
        let mut buf = AlignBuffers::new();
        let res = classify_one(&q(REF), &refs, &p, &mut buf);
        assert_eq!(res.best_ref, Some(2));
        assert_eq!(res.edit, 0);
    }

    #[test]
    fn p_a_prior_scales_by_nraw_abundance_does_not() {
        let bt_prior = Some("Prior".to_string());
        let bt_ab = Some("Abundance".to_string());
        // Prior: p_p * nraw
        assert!((p_a(&bt_prior, Some(2.0e-21), Some(1000.0)).unwrap() - 2.0e-18).abs() < 1e-32);
        // Abundance: unchanged
        assert_eq!(p_a(&bt_ab, Some(5.0e-30), Some(1000.0)), Some(5.0e-30));
        // Prior without nraw -> None (can't reconstruct)
        assert_eq!(p_a(&bt_prior, Some(2.0e-21), None), None);
    }
}
