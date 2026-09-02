//! k-mer-distance vs true-divergence calibration (hidden `kdist-calibrate` subcommand).
//!
//! DADA2's k-mer screen skips alignment for pairs with k-mer distance above
//! `KDIST_CUTOFF = 0.42` (nominally ~10% nucleotide divergence, calibrated on
//! Illumina 16S). The k-mer distance traces to ESPRIT (Sun et al. 2009, whose
//! reference implementation is gone); the 0.42/10% calibration to DADA2
//! (Callahan et al. 2016). This re-derives the relationship empirically on real
//! data: for sampled unique-sequence pairs it emits the k-mer distance
//! (`kmer_dist8`, our port of the ESPRIT metric) alongside the true UNBANDED
//! `align_endsfree` divergence, so the constant can be checked per dataset /
//! platform / k / pooling regime.
//!
//! Alignment is unbanded by default (`--band -1`): the curve must measure the
//! true divergence of distant pairs, which a band would truncate. It is the
//! expensive part, so it is parallelised across `--threads`.
//!
//! POOLING: with multiple derep JSONs, `--per-sample` computes pairs WITHIN each
//! sample (the per-sample / independent regime); the default pools all uniques
//! into one set and computes pairs across the union (the full-pool regime).
//! Pseudo's screen population is per-sample (priors change the partition, not
//! which pairs are screened), so model it with `--per-sample`.

use crate::kmers::{assign_kmer8, kmer_dist8};
use crate::minimizers::{self, MinimizerSketch};
use crate::misc::WithPath as _;
use crate::nwalign::ScreenBackend;
use crate::nwalign::{AlignBuffers, align_endsfree_with_buf};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const GAP: u8 = b'-';

/// The pre-alignment screen under calibration, in either backend's
/// representation (issue: minimizer screen evaluation).
///
/// The whole point of this subcommand is to re-derive `KDIST_CUTOFF` empirically
/// rather than inherit it, and a cutoff turns out **not** to transfer between
/// screens that merely share a distance formula: on the MiSeq SOP, 0.42 passes
/// 27.6% of pairs on the frequency vector and 9.0% on the minimizer sketch, with
/// the matching operating point near 0.64. So the minimizer backend needs its
/// own calibration curve, produced by the same instrument, against the same
/// unbanded true divergence. See `docs/findings/minimizer-screening.md`.
enum Screens {
    Kmer(Vec<Vec<u8>>),
    Minimizer(Vec<MinimizerSketch>),
}

impl Screens {
    fn build(enc: &[Vec<u8>], p: &Params) -> Screens {
        match p.screen_backend {
            ScreenBackend::Kmer => {
                Screens::Kmer(enc.iter().map(|e| assign_kmer8(e, p.k)).collect())
            }
            ScreenBackend::Minimizer => Screens::Minimizer(
                enc.iter()
                    .map(|e| minimizers::sketch(e, p.minimizer_k, p.minimizer_w))
                    .collect(),
            ),
        }
    }

    fn extend(&mut self, other: Screens) {
        match (self, other) {
            (Screens::Kmer(a), Screens::Kmer(b)) => a.extend(b),
            (Screens::Minimizer(a), Screens::Minimizer(b)) => a.extend(b),
            _ => unreachable!("mixed screen backends (uniform per run)"),
        }
    }

    /// Screen distance between element `i` of this set and element `j` of
    /// `other`. Mirrors what `raw_align_with_buf` computes for the same pair,
    /// including the minimizer path's fail-open rule, so the curve describes the
    /// screen that actually runs.
    #[inline]
    fn dist(
        &self,
        i: usize,
        other: &Screens,
        j: usize,
        len_i: usize,
        len_j: usize,
        k: usize,
    ) -> f64 {
        match (self, other) {
            (Screens::Kmer(a), Screens::Kmer(b)) => kmer_dist8(&a[i], len_i, &b[j], len_j, k),
            (Screens::Minimizer(a), Screens::Minimizer(b)) => minimizers::screen_dist(
                minimizers::shared_count(&a[i], &b[j]),
                Some(&a[i]),
                Some(&b[j]),
            ),
            _ => unreachable!("mixed screen backends (uniform per run)"),
        }
    }
}

/// Parameters for [`run`] (mirrors the CLI flags).
#[derive(Clone)]
pub struct Params {
    pub k: usize,
    /// Which screen to calibrate. `Kmer` reproduces the historical behaviour.
    pub screen_backend: ScreenBackend,
    /// Minimizer sketch parameters; ignored unless `screen_backend` is
    /// `Minimizer`.
    pub minimizer_k: usize,
    pub minimizer_w: usize,
    pub cutoff: f64,
    pub leak_pct: f64,
    pub band: i32,
    pub max_pairs: usize,
    pub max_uniques: usize,
    pub per_sample: bool,
    pub nearest_parent: bool,
    /// Post-inference mode: positional inputs are `dada` output JSONs; each is
    /// paired with its derep JSON (from `derep_dir`) so every input unique can be
    /// labelled by what denoising did to it (center / member / failed).
    pub from_dada: bool,
    /// Pooled post-inference mode: positional inputs are the pooled record(s)
    /// written by `dada-pooled --pooled-record <path>`. Self-contained (merged
    /// uniques with pooled abundance + global map + global ASVs), so no
    /// `derep_dir` is needed and no per-sample projection is re-aggregated — the
    /// pool is one population.
    pub from_dada_pooled: bool,
    /// Directory holding the derep JSONs that fed `dada`. Matched to each dada
    /// output by sample name: an exact `{sample}.json[.gz]` first, else a
    /// `{sample}.*.json[.gz]` file (pipeline-renamed dereps, issue #72).
    /// Required with `from_dada`.
    pub derep_dir: Option<PathBuf>,
    pub threads: usize,
    pub seed: u64,
    pub output: Option<PathBuf>,
    pub verbose: bool,
    /// Derive-only mode: report the minimizer cutoff that matches the k-mer
    /// screen's pass rate, and nothing else. Skips alignment entirely, so it runs
    /// in seconds where the full curve takes hours.
    pub derive_cutoff: bool,
    /// Sample pairs uniformly instead of abundance-weighted. Uniform is what a
    /// calibration CURVE wants (it describes the metric); weighted is what a PASS
    /// RATE wants (it describes the comparisons `b_compare` performs). Kept for
    /// comparison, since the published curves are uniform.
    pub derive_uniform_pairs: bool,
}

fn encode(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|b| match b {
            b'A' | b'a' => 1,
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            b'T' | b't' => 4,
            _ => 5, // N etc. — never a valid k-mer, never matches
        })
        .collect()
}

/// Stats of an ends-free alignment. Returns (edits, core_len, band_req):
/// - edits/core_len: substitution+indel columns over the aligned core (terminal
///   gap overhang trimmed — that's length difference, not divergence).
/// - band_req: the MINIMUM diagonal band that would reproduce this alignment =
///   max over the path of |#seq1-bases − #seq2-bases consumed so far|. A banded
///   aligner with band < band_req cannot find this alignment (it gets truncated
///   to a worse score), so band_req is the cost of correctly aligning the pair.
fn aln_divergence(a: &[u8], b: &[u8]) -> (usize, usize, usize) {
    let n = a.len();
    let mut lo = 0;
    while lo < n && (a[lo] == GAP || b[lo] == GAP) {
        lo += 1;
    }
    let mut hi = n;
    while hi > lo && (a[hi - 1] == GAP || b[hi - 1] == GAP) {
        hi -= 1;
    }
    // band_req over the full path (the band applies to the whole DP matrix):
    // a gap in seq2 advances only seq1 (offset +1); a gap in seq1, only seq2
    // (offset -1); a match/mismatch advances both (no change).
    let (mut off, mut band_req) = (0i32, 0i32);
    for k in 0..n {
        if b[k] == GAP {
            off += 1;
        } else if a[k] == GAP {
            off -= 1;
        }
        band_req = band_req.max(off.abs());
    }
    let mut edits = 0;
    for k in lo..hi {
        if a[k] == GAP || b[k] == GAP || a[k] != b[k] {
            edits += 1;
        }
    }
    (edits, hi - lo, band_req as usize)
}

/// One sample's uniques: encoded sequences, abundances, and k-mer screens.
struct Sample {
    name: String,
    enc: Vec<Vec<u8>>,
    counts: Vec<u64>,
    kmers: Screens,
}

/// Load a derep JSON (`uniques[].sequence` + `count`), gzip-transparent. Only
/// derep JSONs are accepted (the screen operates on per-sample uniques).
fn load_derep(path: &Path, p: &Params, max_uniques: usize, seed: u64) -> io::Result<Sample> {
    let f = File::open(path).with_path(path)?;
    let mut txt = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        MultiGzDecoder::new(f).read_to_string(&mut txt)?;
    } else {
        io::BufReader::new(f).read_to_string(&mut txt)?;
    }
    let v: serde_json::Value =
        serde_json::from_str(&txt).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: not a derep JSON (no `uniques`)", path.display()),
        )
    })?;
    let name = v
        .get("sample")
        .and_then(|s| s.as_str())
        .map(str::to_string)
        .unwrap_or_else(|| {
            path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("?")
                .to_string()
        });
    let mut seqs: Vec<(Vec<u8>, u64)> = uniques
        .iter()
        .filter_map(|e| {
            let s = e.get("sequence")?.as_str()?;
            let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
            Some((encode(s), c))
        })
        .collect();
    // Optional per-sample subsample of uniques to bound the O(n^2) pair count.
    // Random (not abundance-top) keeps the divergence distribution unbiased.
    if max_uniques > 0 && seqs.len() > max_uniques {
        let mut st = seed ^ (seqs.len() as u64).wrapping_mul(0x9E37_79B9);
        // partial Fisher–Yates: move `max_uniques` random picks to the front.
        for i in 0..max_uniques {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
            let j = i + ((st >> 33) as usize) % (seqs.len() - i);
            seqs.swap(i, j);
        }
        seqs.truncate(max_uniques);
    }
    let (enc, counts): (Vec<_>, Vec<_>) = seqs.into_iter().unzip();
    let kmers = Screens::build(&enc, p);
    Ok(Sample {
        name,
        enc,
        counts,
        kmers,
    })
}

/// A `dada` output paired with its derep input, ready for post-inference
/// classification. Input uniques are in denoising input order (abundance-desc),
/// so `map[i]` is the cluster index Raw `i` was corrected to (`None` = failed
/// the abundance test, i.e. did not survive denoising). Centers are indexed by
/// cluster id, carrying the birth metadata that lets us trace *why* an ASV
/// exists (`Prior` = born from a pseudo-pool prior; `birth_pval` = how close
/// that birth was to OMEGA_A).
struct DadaSample {
    name: String,
    // input uniques (denoising input order, aligned with `map`)
    enc: Vec<Vec<u8>>,
    counts: Vec<u64>,
    kmers: Screens,
    map: Vec<Option<usize>>,
    // cluster centers, indexed by cluster id
    c_enc: Vec<Vec<u8>>,
    c_kmers: Screens,
    c_ab: Vec<u64>,
    c_birth: Vec<String>,
    c_birth_pval: Vec<f64>,
}

/// Read a derep JSON's uniques in denoising input order (abundance-descending,
/// matching `load_derep_for_dada`): no subsample — indices must line up with the
/// dada `map`. Returns (encoded seqs, counts).
fn load_derep_aligned(path: &Path) -> io::Result<(Vec<Vec<u8>>, Vec<u64>)> {
    let f = File::open(path).with_path(path)?;
    let mut txt = String::new();
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        MultiGzDecoder::new(f).read_to_string(&mut txt)?;
    } else {
        io::BufReader::new(f).read_to_string(&mut txt)?;
    }
    let v: serde_json::Value =
        serde_json::from_str(&txt).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: not a derep JSON (no `uniques`)", path.display()),
        )
    })?;
    let mut seqs: Vec<(Vec<u8>, u64)> = uniques
        .iter()
        .filter_map(|e| {
            let s = e.get("sequence")?.as_str()?;
            let c = e.get("count").and_then(|c| c.as_u64()).unwrap_or(1);
            Some((encode(s), c))
        })
        .collect();
    // Mirror load_derep_for_dada: only re-sort when the producer didn't declare
    // abundance_desc; the sort is stable so ties keep file order (= input order).
    if v.get("sort_order").and_then(|s| s.as_str()) != Some("abundance_desc") {
        seqs.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    }
    Ok(seqs.into_iter().unzip())
}

/// Read just the `sample` field of a derep JSON (gzip-transparent), for
/// disambiguating candidates. Returns `None` if unreadable or absent.
fn derep_sample_name(path: &Path) -> Option<String> {
    let f = File::open(path).ok()?;
    let reader: Box<dyn Read> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    let v: serde_json::Value = serde_json::from_reader(io::BufReader::new(reader)).ok()?;
    v.get("sample").and_then(|s| s.as_str()).map(str::to_string)
}

/// Resolve the derep JSON for sample `name` inside `derep_dir`.
///
/// `dada-pooled` names its per-sample outputs by the internal sample name, but
/// derep inputs are commonly renamed by the pipeline (e.g. an `S1` sample whose
/// derep is `S1.derep.R1.json.gz`), so an exact `{name}.json[.gz]` match is not
/// guaranteed (issue #72). Resolution order:
///   1. exact `{name}.json` / `{name}.json.gz` (fast, backward compatible);
///   2. otherwise scan the directory for `*.json` / `*.json.gz` whose leading
///      dot-delimited component equals `name`; if exactly one matches use it,
///      and if several do (e.g. R1 and R2 in the same dir) disambiguate by the
///      derep JSON's own `sample` field.
fn find_derep_for_sample(derep_dir: &Path, name: &str) -> io::Result<PathBuf> {
    // 1. Exact match — cheapest, preserves prior behaviour.
    for cand in [
        derep_dir.join(format!("{name}.json")),
        derep_dir.join(format!("{name}.json.gz")),
    ] {
        if cand.exists() {
            return Ok(cand);
        }
    }

    // 2. Glob-like scan for `{name}.*.json[.gz]`.
    let mut matches: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(derep_dir).with_path(derep_dir)? {
        let path = entry?.path();
        let Some(fname) = path.file_name().and_then(|s| s.to_str()) else {
            continue;
        };
        let Some(stem) = fname
            .strip_suffix(".json.gz")
            .or_else(|| fname.strip_suffix(".json"))
        else {
            continue;
        };
        // Whole stem, or its leading dot-component (`S1.derep.R1` → `S1`). Exact
        // component equality avoids `S188` matching an `S1888` file.
        if stem == name || stem.split('.').next() == Some(name) {
            matches.push(path);
        }
    }

    match matches.len() {
        0 => Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "no derep JSON for sample `{name}` in {}",
                derep_dir.display()
            ),
        )),
        1 => Ok(matches.pop().unwrap()),
        _ => {
            // Several filenames share the prefix (e.g. R1/R2 in one dir):
            // fall back to the authoritative internal `sample` field.
            let mut exact: Vec<PathBuf> = matches
                .iter()
                .filter(|p| derep_sample_name(p).as_deref() == Some(name))
                .cloned()
                .collect();
            match exact.len() {
                1 => Ok(exact.pop().unwrap()),
                0 => Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "no derep JSON whose `sample` is `{name}` among {} name-prefixed \
                         candidate(s) in {} (checked each file's internal `sample` field)",
                        matches.len(),
                        derep_dir.display()
                    ),
                )),
                n => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "ambiguous derep JSON for sample `{name}` in {}: {n} files declare \
                         `sample = {name}`",
                        derep_dir.display()
                    ),
                )),
            }
        }
    }
}

/// Load a `dada` output JSON and pair it with its derep input from `derep_dir`.
fn load_dada(dada_path: &Path, derep_dir: &Path, p: &Params) -> io::Result<DadaSample> {
    let f = File::open(dada_path).with_path(dada_path)?;
    let reader: Box<dyn Read> = if dada_path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    let v: serde_json::Value = serde_json::from_reader(io::BufReader::new(reader))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let name = v
        .get("sample")
        .and_then(|s| s.as_str())
        .map(str::to_string)
        .unwrap_or_else(|| {
            dada_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("?")
                .to_string()
        });
    let asvs = v.get("asvs").and_then(|a| a.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: not a dada output JSON (no `asvs`)",
                dada_path.display()
            ),
        )
    })?;
    let (mut c_enc, mut c_ab, mut c_birth, mut c_birth_pval) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for a in asvs {
        let seq = a.get("sequence").and_then(|s| s.as_str()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "asv entry missing `sequence`")
        })?;
        c_enc.push(encode(seq));
        c_ab.push(a.get("abundance").and_then(|x| x.as_u64()).unwrap_or(0));
        c_birth.push(
            a.get("birth_type")
                .and_then(|s| s.as_str())
                .unwrap_or("?")
                .to_string(),
        );
        c_birth_pval.push(
            a.get("birth_pval")
                .and_then(|x| x.as_f64())
                .unwrap_or(f64::NAN),
        );
    }
    let map: Vec<Option<usize>> = v
        .get("map")
        .and_then(|m| m.as_array())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "dada output missing `map`"))?
        .iter()
        .map(|e| e.as_u64().map(|u| u as usize))
        .collect();

    // Locate and load the matching derep input.
    let derep_path = find_derep_for_sample(derep_dir, &name)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", dada_path.display())))?;
    let (enc, counts) = load_derep_aligned(&derep_path)?;
    if enc.len() != map.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{name}: derep has {} uniques but dada map has {} entries — \
                 mismatched inputs (regenerate from the same derep)",
                enc.len(),
                map.len()
            ),
        ));
    }
    let kmers = Screens::build(&enc, p);
    let c_kmers = Screens::build(&c_enc, p);
    Ok(DadaSample {
        name,
        enc,
        counts,
        kmers,
        map,
        c_enc,
        c_kmers,
        c_ab,
        c_birth,
        c_birth_pval,
    })
}

/// Load a `dada-pooled` `_pooled.json[.gz]` record into a single [`DadaSample`].
///
/// Unlike [`load_dada`], the record is self-contained: it carries the merged
/// uniques (with their POOLED abundance) inline alongside the global `map` and
/// ASV list, so there is no derep JSON to locate and no per-sample projection to
/// reconcile. Emitted in denoising input order, so `uniques`/`map` line up with
/// no re-sort — the pool is one population named `__pooled__`.
fn load_dada_pooled(path: &Path, p: &Params) -> io::Result<DadaSample> {
    let f = File::open(path).with_path(path)?;
    let reader: Box<dyn Read> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    let v: serde_json::Value = serde_json::from_reader(io::BufReader::new(reader))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let uniques = v.get("uniques").and_then(|u| u.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: not a dada-pooled record (no `uniques`); pass the `_pooled.json[.gz]` \
                 written by `dada-pooled`",
                path.display()
            ),
        )
    })?;
    let (mut enc, mut counts) = (Vec::with_capacity(uniques.len()), Vec::new());
    for e in uniques {
        let seq = e.get("sequence").and_then(|s| s.as_str()).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "pooled unique missing `sequence`",
            )
        })?;
        enc.push(encode(seq));
        counts.push(e.get("count").and_then(|c| c.as_u64()).unwrap_or(1));
    }
    let asvs = v.get("asvs").and_then(|a| a.as_array()).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: pooled record missing `asvs`", path.display()),
        )
    })?;
    let (mut c_enc, mut c_ab, mut c_birth, mut c_birth_pval) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for a in asvs {
        let seq = a.get("sequence").and_then(|s| s.as_str()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "asv entry missing `sequence`")
        })?;
        c_enc.push(encode(seq));
        c_ab.push(a.get("abundance").and_then(|x| x.as_u64()).unwrap_or(0));
        c_birth.push(
            a.get("birth_type")
                .and_then(|s| s.as_str())
                .unwrap_or("?")
                .to_string(),
        );
        c_birth_pval.push(
            a.get("birth_pval")
                .and_then(|x| x.as_f64())
                .unwrap_or(f64::NAN),
        );
    }
    let map: Vec<Option<usize>> = v
        .get("map")
        .and_then(|m| m.as_array())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "pooled record missing `map`"))?
        .iter()
        .map(|e| e.as_u64().map(|u| u as usize))
        .collect();
    if enc.len() != map.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: {} merged uniques but map has {} entries — corrupt pooled record",
                path.display(),
                enc.len(),
                map.len()
            ),
        ));
    }
    let kmers = Screens::build(&enc, p);
    let c_kmers = Screens::build(&c_enc, p);
    Ok(DadaSample {
        name: "__pooled__".into(),
        enc,
        counts,
        kmers,
        map,
        c_enc,
        c_kmers,
        c_ab,
        c_birth,
        c_birth_pval,
    })
}

/// Post-inference driver for pooled runs: load each `_pooled.json[.gz]` record as
/// one pool population and emit the labelled comparisons (no --derep-dir needed).
fn run_from_dada_pooled(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    let samples: Vec<DadaSample> = inputs
        .iter()
        .map(|path| load_dada_pooled(path, p))
        .collect::<io::Result<_>>()?;

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;
    from_dada_mode(&mut *w, &pool, &samples, p)?;
    w.flush()?;
    Ok(())
}

/// Build the (i, j) pair list for a population of `n` uniques: enumerate all if
/// `n*(n-1)/2 <= max_pairs`, else draw `max_pairs` random pairs (with possible
/// repeats — fine for a calibration scatter).
/// Sample pairs the way `b_compare` forms them: one side abundance-weighted.
///
/// `pairs_for` samples uniformly over all pairs, which is right for a calibration
/// curve — it describes the metric. It is **wrong for predicting a pass rate**,
/// because `b_compare` never compares random pairs: it compares every raw against
/// each cluster CENTRE, and centres are the abundant uniques. Measured on pooled
/// PacBio, the minimizer/k-mer pass ratio is 0.744 on the pairs actually screened
/// and 0.911 on uniform random pairs — so uniform sampling makes the minimizer
/// look **23% less selective than it is**, and a cutoff matched on it overshoots.
///
/// One side is drawn from the abundance distribution (the centre proxy) and the
/// other uniformly (the raw), which is the shape of a raw-vs-centre comparison.
fn weighted_pairs_for(counts: &[u64], max_pairs: usize, seed: u64) -> Vec<(usize, usize)> {
    let n = counts.len();
    if n < 2 {
        return Vec::new();
    }
    let mut cum: Vec<u64> = Vec::with_capacity(n);
    let mut tot = 0u64;
    for &c in counts {
        tot += c.max(1);
        cum.push(tot);
    }
    let mut st = seed ^ 0x9E3779B97F4A7C15;
    let mut next = || {
        st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
        st >> 11
    };
    let mut out = Vec::with_capacity(max_pairs);
    let mut guard = 0usize;
    while out.len() < max_pairs && guard < max_pairs * 8 {
        guard += 1;
        let target = next() % tot.max(1);
        let j = cum.partition_point(|&c| c <= target).min(n - 1);
        let i = (next() as usize) % n;
        if i != j {
            out.push((i, j));
        }
    }
    out
}

fn pairs_for(n: usize, max_pairs: usize, seed: u64) -> Vec<(usize, usize)> {
    let total = n.saturating_mul(n.saturating_sub(1)) / 2;
    if n < 2 {
        return Vec::new();
    }
    if total <= max_pairs {
        let mut v = Vec::with_capacity(total);
        for i in 0..n {
            for j in (i + 1)..n {
                v.push((i, j));
            }
        }
        return v;
    }
    let mut st = seed;
    let mut rnd = |m: usize| {
        st = st.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((st >> 33) as usize) % m
    };
    (0..max_pairs)
        .map(|_| {
            let i = rnd(n);
            let mut j = rnd(n);
            if i == j {
                j = (j + 1) % n;
            }
            (i.min(j), i.max(j))
        })
        .collect()
}

/// Derive the minimizer cutoff that reproduces the k-mer screen's pass rate.
///
/// The full calibration curve exists to relate a screen distance to TRUE
/// divergence, and that is what costs the money: every sampled pair is aligned
/// unbanded (~2.2 M DP cells on 1.5 kb reads, ~34 M alignments on the 95-sample
/// PacBio run). The matched-pass rule never consults true divergence, so this
/// mode computes both screens over the same sampled pairs and stops.
///
/// Both screens are built over the SAME encoded uniques and evaluated on the
/// SAME pairs, which is the whole point: a cutoff is a property of the distance
/// distribution a metric induces, so the two distributions have to be compared
/// on identical input or the quantile is meaningless.
fn run_derive_cutoff(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    let mut kp = p.clone();
    kp.screen_backend = ScreenBackend::Kmer;
    let mut mp = p.clone();
    mp.screen_backend = ScreenBackend::Minimizer;

    // Populations must match the denoising mode. The k-mer screen's pass rate is
    // a property of the population, not the dataset: pooled ITS2 passes 0.70% and
    // the same data per-sample passes 1.93%. Deriving from the wrong one targets
    // the selectivity of a screen nobody ran.
    let mut pops: Vec<(String, Vec<Vec<u8>>, Vec<u64>)> = Vec::new();
    for path in inputs {
        let s = load_derep(path, &kp, p.max_uniques, p.seed)?;
        if p.per_sample {
            pops.push((s.name, s.enc, s.counts));
        } else {
            if pops.is_empty() {
                pops.push(("pool".to_string(), Vec::new(), Vec::new()));
            }
            pops[0].1.extend(s.enc);
            pops[0].2.extend(s.counts);
        }
    }
    pops.retain(|(_, e, _)| e.len() >= 2);
    if pops.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no population has 2 or more uniques to sample pairs from",
        ));
    }

    println!(
        "# derive-cutoff: matched pass rate, no alignment performed\n# pairs: {}",
        if p.derive_uniform_pairs {
            "uniform random (describes the metric; NOT the population b_compare screens)"
        } else {
            "abundance-weighted -- the raw-vs-centre shape b_compare actually screens"
        }
    );
    let mut per_pop: Vec<(String, minimizers::MatchedPass)> = Vec::new();
    let mut kd = Vec::new();
    let mut md = Vec::new();
    for (name, enc, counts) in &pops {
        let kmer = Screens::build(enc, &kp);
        let mini = Screens::build(enc, &mp);
        let pairs = if p.derive_uniform_pairs {
            pairs_for(enc.len(), p.max_pairs, p.seed)
        } else {
            weighted_pairs_for(counts, p.max_pairs, p.seed)
        };
        let mut pkd = Vec::with_capacity(pairs.len());
        let mut pmd = Vec::with_capacity(pairs.len());
        for &(i, j) in &pairs {
            let (li, lj) = (enc[i].len(), enc[j].len());
            pkd.push(kmer.dist(i, &kmer, j, li, lj, p.k));
            pmd.push(mini.dist(i, &mini, j, li, lj, p.k));
        }
        per_pop.push((
            name.clone(),
            minimizers::matched_pass_cutoff(&pkd, &pmd, p.cutoff),
        ));
        kd.extend(pkd);
        md.extend(pmd);
    }

    if p.per_sample && per_pop.len() > 1 {
        // Per-sample spread is the thing to look at: a rule that has to be
        // applied once per run is only usable if the samples agree on it.
        println!("\nper-sample derivation ({} populations):", per_pop.len());
        println!(
            "  {:>28} {:>8} {:>10} {:>8}",
            "sample", "uniques", "kmer pass", "cutoff"
        );
        for ((name, r), (_, enc, _)) in per_pop.iter().zip(pops.iter()) {
            println!(
                "  {:>28} {:>8} {:>9.4}% {:>8.2}",
                &name[name.len().saturating_sub(28)..],
                enc.len(),
                100.0 * r.kmer_pass,
                r.cutoff_rounded
            );
        }
        let mut cs: Vec<f64> = per_pop.iter().map(|(_, r)| r.cutoff_rounded).collect();
        cs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "  per-sample cutoff: median {:.2}, range {:.2}-{:.2}{}",
            cs[cs.len() / 2],
            cs[0],
            cs[cs.len() - 1],
            if cs[cs.len() - 1] - cs[0] > 0.04 {
                "  <- SPREAD > 0.04: one run-wide cutoff is a compromise here"
            } else {
                ""
            }
        );
    }

    let r = minimizers::matched_pass_cutoff(&kd, &md, p.cutoff);
    println!(
        "\nuniques      {} in {} population(s)\npairs        {}",
        pops.iter().map(|(_, e, _)| e.len()).sum::<usize>(),
        pops.len(),
        r.n_pairs
    );
    println!(
        "k-mer screen k={} @ cutoff {:.2}: passes {:.4}% of pairs",
        p.k,
        p.cutoff,
        100.0 * r.kmer_pass
    );
    println!(
        "minimizer    k={} w={}: matching cutoff {:.4} -> **{:.2}** (passes {:.4}%)",
        p.minimizer_k,
        p.minimizer_w,
        r.cutoff,
        r.cutoff_rounded,
        100.0 * r.mini_pass
    );
    // Pair-level AGREEMENT with the k-mer screen at the derived cutoff.
    //
    // Matching the pass RATE says the two screens admit the same NUMBER of pairs.
    // It says nothing about whether they admit the SAME pairs, and that is what
    // decides churn: a pair the k-mer screen aligns and the minimizer shrouds is a
    // raw that may never reach its parent cluster. This is the pair-level recall
    // gate this investigation opened with, and it costs nothing here -- both
    // distances are already computed for every sampled pair.
    //
    // It is also the alignment-free way to compare (k, w) settings. The standing
    // hypothesis for why Illumina churns and HiFi does not is estimator variance
    // from a small sketch -- 74 entries/raw on 231 bp ITS2 against 475 on 1.5 kb
    // HiFi. If that is right, lowering w raises the entry count and DISAGREEMENT
    // here should fall. That is testable in seconds per setting, where the ASV
    // table costs a full pipeline run per arm.
    let c = r.cutoff_rounded;
    let (mut both, mut neither, mut konly, mut monly) = (0usize, 0usize, 0usize, 0usize);
    for (&a, &b) in kd.iter().zip(md.iter()) {
        match (a <= p.cutoff, b <= c) {
            (true, true) => both += 1,
            (false, false) => neither += 1,
            (true, false) => konly += 1,
            (false, true) => monly += 1,
        }
    }
    let n = r.n_pairs.max(1);
    println!("\npair-level agreement with the k-mer screen at cutoff {c:.2}:");
    println!(
        "  both pass          {both:>9} ({:6.3}%)",
        100.0 * both as f64 / n as f64
    );
    println!(
        "  both shroud        {neither:>9} ({:6.3}%)",
        100.0 * neither as f64 / n as f64
    );
    println!(
        "  k-mer only         {konly:>9} ({:6.3}%)  <- SHROUDED by the minimizer; the churn risk",
        100.0 * konly as f64 / n as f64
    );
    println!(
        "  minimizer only     {monly:>9} ({:6.3}%)  <- extra alignments, cost not correctness",
        100.0 * monly as f64 / n as f64
    );
    println!(
        "  DISAGREEMENT       {:>9} ({:6.3}% of pairs)   recall vs k-mer: {:6.3}%",
        konly + monly,
        100.0 * (konly + monly) as f64 / n as f64,
        100.0 * both as f64 / (both + konly).max(1) as f64
    );

    // The pass-rate curve either side, so how sharply the choice matters is
    // visible rather than implied. A flat neighbourhood means the exact value
    // hardly matters; a steep one means it does.
    println!("\nminimizer pass rate near the derived cutoff:");
    for step in -4i32..=4 {
        let c = r.cutoff_rounded + f64::from(step) * 0.01;
        if c <= 0.0 || c >= 1.0 {
            continue;
        }
        let pass = md.iter().filter(|&&d| d <= c).count() as f64 / r.n_pairs as f64;
        println!(
            "  {:.2}{} {:8.4}%",
            c,
            if step == 0 { " <-" } else { "   " },
            100.0 * pass
        );
    }
    println!(
        "\nNOTE: this reproduces the k-mer screen's SELECTIVITY, which is the safe\n\
         target, not the cheapest cutoff that still agrees with it. On PacBio HiFi\n\
         the ASV table is identical from 0.45 to 0.60 and this rule picks 0.50 --\n\
         inside the plateau but at its expensive end, worth ~19% of wall clock.\n\
         Finding the cheap edge of a plateau needs the ASV table, so sweep if you\n\
         can afford to."
    );
    Ok(())
}

pub fn run(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    if inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no input derep JSON(s) given",
        ));
    }
    if p.from_dada_pooled {
        return run_from_dada_pooled(inputs, p);
    }
    if p.from_dada {
        return run_from_dada(inputs, p);
    }
    if p.derive_cutoff {
        return run_derive_cutoff(inputs, p);
    }
    let loaded: Vec<Sample> = inputs
        .iter()
        .map(|path| load_derep(path, p, p.max_uniques, p.seed))
        .collect::<io::Result<_>>()?;

    // Form populations: one per sample (per-sample) or a single merged pool.
    let pops: Vec<Sample> = if p.per_sample {
        loaded
    } else {
        let mut pool = Sample {
            name: "pool".into(),
            enc: Vec::new(),
            counts: Vec::new(),
            kmers: Screens::build(&[], p),
        };
        for s in loaded {
            pool.enc.extend(s.enc);
            pool.counts.extend(s.counts);
            pool.kmers.extend(s.kmers);
        }
        vec![pool]
    };

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;

    if p.nearest_parent {
        nearest_parent_mode(&mut *w, &pool, &pops, p)?;
    } else {
        pairs_mode(&mut *w, &pool, &pops, p)?;
    }
    w.flush()?;
    Ok(())
}

/// All-pairs (random-subsampled) mode: kdist vs divergence over sampled pairs.
fn pairs_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    pops: &[Sample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,kdist,edits,core_len,pct_div,band_req,screened_in,ab_i,ab_j"
    )?;
    let (mut tot, mut scr, mut leak) = (0u64, 0u64, 0u64);
    let mut band_all: Vec<usize> = Vec::new();
    let mut band_scr: Vec<usize> = Vec::new();
    for (pi, s) in pops.iter().enumerate() {
        let n = s.enc.len();
        let seed = p.seed.wrapping_add(pi as u64).wrapping_mul(0x100000001B3);
        let pairs = pairs_for(n, p.max_pairs, seed);
        if p.verbose {
            let total = n.saturating_mul(n.saturating_sub(1)) / 2;
            eprintln!(
                "[kdist] {} : {n} uniques, {total} pairs -> {} computed (k={}, band={}, {} threads)",
                s.name,
                pairs.len(),
                p.k,
                p.band,
                p.threads,
            );
        }
        let rows: Vec<(f64, usize, usize, f64, usize, bool)> = pool.install(|| {
            pairs
                .par_iter()
                .map_init(AlignBuffers::new, |buf, &(i, j)| {
                    let kd = s
                        .kmers
                        .dist(i, &s.kmers, j, s.enc[i].len(), s.enc[j].len(), p.k);
                    align_endsfree_with_buf(&s.enc[i], &s.enc[j], 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (kd, edits, core, pct, band_req, kd < p.cutoff)
                })
                .collect()
        });
        for (idx, &(kd, edits, core, pct, band_req, sin)) in rows.iter().enumerate() {
            let (i, j) = pairs[idx];
            tot += 1;
            band_all.push(band_req);
            if sin {
                scr += 1;
                band_scr.push(band_req);
                if pct > p.leak_pct {
                    leak += 1;
                }
            }
            writeln!(
                w,
                "{},{kd:.4},{edits},{core},{pct:.3},{band_req},{},{},{}",
                s.name, sin as u8, s.counts[i], s.counts[j]
            )?;
        }
    }
    if p.verbose && tot > 0 {
        eprintln!(
            "[kdist] {tot} pairs: screened-in (kdist<{}) {scr} ({:.1}%); of those {leak} are >{}% divergent (leakage)",
            p.cutoff,
            100.0 * scr as f64 / tot as f64,
            p.leak_pct,
        );
        eprintln!("[kdist] {}", band_fit("all pairs", &band_all));
        eprintln!("[kdist] {}", band_fit("screened-in", &band_scr));
    }
    Ok(())
}

/// Post-inference driver: load each (dada output, derep input) pair and emit the
/// labelled pairwise comparisons.
fn run_from_dada(inputs: &[PathBuf], p: &Params) -> io::Result<()> {
    let derep_dir = p.derep_dir.as_deref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "--from-dada requires --derep-dir <DIR> (where the matching derep JSONs live)",
        )
    })?;
    let samples: Vec<DadaSample> = inputs
        .iter()
        .map(|path| load_dada(path, derep_dir, p))
        .collect::<io::Result<_>>()?;

    let mut w: Box<dyn Write> = match &p.output {
        Some(path) => Box::new(BufWriter::new(File::create(path).with_path(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(p.threads)
        .build()
        .map_err(io::Error::other)?;
    from_dada_mode(&mut *w, &pool, &samples, p)?;
    w.flush()?;
    Ok(())
}

/// A row of work: align input unique `i` (or center `i` for center pairs) against
/// a partner center `c`, under a class label.
enum Job {
    /// Member `i` absorbed by its own center `c` (an error copy denoising fixed).
    Member { i: usize, c: usize },
    /// Failed unique `i` (map==None) vs its nearest center `c` by k-mer distance.
    Failed { i: usize, c: usize },
    /// Two surviving ASV centers `a`,`b` (the inter-ASV resolution floor).
    CenterPair { a: usize, b: usize },
}

/// Post-inference classification mode. For every input unique, denoising decided
/// one of three fates; we align it against the relevant center and tag the row so
/// the divergence/k-mer-distance can be read *conditioned on what dada did*:
///   - member  → its own center: real error copies, the within-cluster cloud.
///   - failed  → nearest center: shed by the abundance test but not assigned.
///   - center  → other centers: survivors vs each other (resolution floor).
///
/// `birth_type`/`birth_pval` ride along on the partner center so prior-born ASVs
/// (pseudo-pool) and births near OMEGA_A are visible in the same table.
fn from_dada_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    samples: &[DadaSample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,class,cluster,ab,center_ab,ab_ratio,birth_type,birth_pval,\
         kdist,edits,core_len,pct_div,band_req,screened_in"
    )?;
    for s in samples {
        let ncenters = s.c_enc.len();
        // Build the job list from each unique's fate.
        let mut jobs: Vec<Job> = Vec::with_capacity(s.enc.len());
        let (mut n_center, mut n_member, mut n_failed) = (0u64, 0u64, 0u64);
        for (i, m) in s.map.iter().enumerate() {
            match m {
                None => {
                    // nearest center by k-mer distance (cheap pre-pass, no align)
                    if ncenters > 0 {
                        let c = (0..ncenters)
                            .min_by(|&a, &b| {
                                let ka = s.kmers.dist(
                                    i,
                                    &s.c_kmers,
                                    a,
                                    s.enc[i].len(),
                                    s.c_enc[a].len(),
                                    p.k,
                                );
                                let kb = s.kmers.dist(
                                    i,
                                    &s.c_kmers,
                                    b,
                                    s.enc[i].len(),
                                    s.c_enc[b].len(),
                                    p.k,
                                );
                                ka.partial_cmp(&kb).unwrap()
                            })
                            .unwrap();
                        jobs.push(Job::Failed { i, c });
                    }
                    n_failed += 1;
                }
                Some(c) => {
                    let c = *c;
                    if s.enc[i] == s.c_enc[c] {
                        n_center += 1; // the representative itself — no self-align
                    } else {
                        jobs.push(Job::Member { i, c });
                        n_member += 1;
                    }
                }
            }
        }
        // Inter-center pairs (resolution floor): enumerate / subsample.
        let cpairs = pairs_for(ncenters, p.max_pairs, p.seed ^ 0xC0FFEE);
        for (a, b) in cpairs {
            jobs.push(Job::CenterPair { a, b });
        }
        if p.verbose {
            eprintln!(
                "[kdist] {} : {} uniques ({n_center} centers, {n_member} members, {n_failed} failed), \
                 {ncenters} ASVs, {} jobs (k={}, band={}, {} threads)",
                s.name,
                s.enc.len(),
                jobs.len(),
                p.k,
                p.band,
                p.threads,
            );
        }
        // Align every job in parallel.
        let rows: Vec<(usize, f64, usize, usize, f64, usize)> = pool.install(|| {
            jobs.par_iter()
                .map_init(AlignBuffers::new, |buf, job| {
                    // (screen set, index) pairs rather than borrowed slices:
                    // the sketch backend has no slice to borrow.
                    let (ei, ej, si, i_idx, sj, j_idx, c) = match job {
                        Job::Member { i, c } => {
                            (&s.enc[*i], &s.c_enc[*c], &s.kmers, *i, &s.c_kmers, *c, *c)
                        }
                        Job::Failed { i, c } => {
                            (&s.enc[*i], &s.c_enc[*c], &s.kmers, *i, &s.c_kmers, *c, *c)
                        }
                        Job::CenterPair { a, b } => (
                            &s.c_enc[*a],
                            &s.c_enc[*b],
                            &s.c_kmers,
                            *a,
                            &s.c_kmers,
                            *b,
                            *b,
                        ),
                    };
                    let kd = si.dist(i_idx, sj, j_idx, ei.len(), ej.len(), p.k);
                    align_endsfree_with_buf(ei, ej, 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (c, kd, edits, core, pct, band_req)
                })
                .collect()
        });
        // Emit, pairing each row back with its job for labels/abundances.
        // Track the failed class by abundance: singletons can never seed an ASV
        // under the default (≥2 reads, unless --detect-singletons), so a failed
        // singleton fails for that reason, NOT for being distant — split it out.
        let (mut f_singleton, mut f_singleton_in, mut f_multi, mut f_multi_in) =
            (0u64, 0u64, 0u64, 0u64);
        for (job, &(c, kd, edits, core, pct, band_req)) in jobs.iter().zip(&rows) {
            let (class, ab, center_ab) = match job {
                Job::Member { i, .. } => ("member", s.counts[*i], s.c_ab[c]),
                Job::Failed { i, .. } => ("failed", s.counts[*i], s.c_ab[c]),
                Job::CenterPair { a, .. } => ("center_pair", s.c_ab[*a], s.c_ab[c]),
            };
            if let Job::Failed { .. } = job {
                let within = (kd < p.cutoff) as u64;
                if ab <= 1 {
                    f_singleton += 1;
                    f_singleton_in += within;
                } else {
                    f_multi += 1;
                    f_multi_in += within;
                }
            }
            let ratio = center_ab as f64 / ab.max(1) as f64;
            writeln!(
                w,
                "{},{class},{c},{ab},{center_ab},{ratio:.2},{},{:.3e},{kd:.4},{edits},{core},{pct:.3},{band_req},{}",
                s.name,
                s.c_birth[c],
                s.c_birth_pval[c],
                (kd < p.cutoff) as u8,
            )?;
        }
        if p.verbose {
            let f_total = f_singleton + f_multi;
            if f_total > 0 {
                eprintln!(
                    "[kdist] {} : {f_total} failed | singletons {f_singleton} ({f_singleton_in} within cutoff) \
                     | multi-read {f_multi} ({f_multi_in} within cutoff) — failed singletons are the \
                     --detect-singletons tradeoff, not distance",
                    s.name,
                );
            }
            let priors = s.c_birth.iter().filter(|b| b.as_str() == "Prior").count();
            if priors > 0 {
                eprintln!(
                    "[kdist] {} : {priors}/{ncenters} ASVs born from priors (pseudo); \
                     filter the table on class=center_pair,birth_type=Prior to see their nearest survivor",
                    s.name,
                );
            }
        }
    }
    Ok(())
}

/// Candidate band sizes for the band-fit summary (DADA2 default is 16).
const BAND_SWEEP: [usize; 7] = [2, 4, 8, 16, 32, 64, 128];

/// For each candidate band B, the fraction of alignments whose true path fits
/// within B (band_req <= B) — i.e. that a banded aligner at B would compute
/// correctly. The complement is distorted/effectively-rejected by that band.
fn band_fit(label: &str, band_reqs: &[usize]) -> String {
    let n = band_reqs.len();
    if n == 0 {
        return format!("{label} band-fit: (none)");
    }
    let parts: Vec<String> = BAND_SWEEP
        .iter()
        .map(|&b| {
            let f = band_reqs.iter().filter(|&&r| r <= b).count();
            format!("≤{b}:{:.1}%", 100.0 * f as f64 / n as f64)
        })
        .collect();
    let mx = band_reqs.iter().copied().max().unwrap_or(0);
    format!("{label} band-fit ({n}, max_req {mx}): {}", parts.join(" "))
}

/// Divergence below which a nearest-parent link is treated as a clear
/// error-copy candidate when computing the screen "headroom".
const CLEAR_ERROR_COPY_PCT: f64 = 3.0;

/// Abundance-aware mode: for each unique, find its nearest MORE-abundant
/// neighbour (its candidate error-copy "parent", mirroring DADA2's greedy
/// center-based comparison) by k-mer distance, then align that one pair for the
/// true divergence. The distribution of parent-link kdist is the empirical
/// error-copy distance ceiling; `cutoff − ceiling` is the screen's headroom.
fn nearest_parent_mode(
    w: &mut dyn Write,
    pool: &rayon::ThreadPool,
    pops: &[Sample],
    p: &Params,
) -> io::Result<()> {
    writeln!(
        w,
        "sample,ab,parent_ab,ab_ratio,kdist,edits,core_len,pct_div,band_req,screened_in"
    )?;
    for s in pops {
        let n = s.enc.len();
        // abundance-descending order (stable by index for ties)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| s.counts[b].cmp(&s.counts[a]).then(a.cmp(&b)));
        if p.verbose {
            eprintln!(
                "[kdist] {} : {n} uniques, nearest more-abundant parent each (k={}, band={}, {} threads)",
                s.name, p.k, p.band, p.threads,
            );
        }
        // For each non-top unique (position r), scan the more-abundant prefix
        // order[0..r] for the min-kdist parent, then align that pair.
        let rows: Vec<(usize, usize, f64, usize, usize, f64, usize)> = pool.install(|| {
            (1..n)
                .into_par_iter()
                .map_init(AlignBuffers::new, |buf, r| {
                    let i = order[r];
                    let (mut best_kd, mut parent) = (f64::INFINITY, order[0]);
                    for &c in &order[0..r] {
                        let kd = s
                            .kmers
                            .dist(i, &s.kmers, c, s.enc[i].len(), s.enc[c].len(), p.k);
                        if kd < best_kd {
                            best_kd = kd;
                            parent = c;
                        }
                    }
                    align_endsfree_with_buf(&s.enc[i], &s.enc[parent], 5, -4, -8, p.band, buf);
                    let (edits, core, band_req) = aln_divergence(&buf.al0, &buf.al1);
                    let pct = if core > 0 {
                        100.0 * edits as f64 / core as f64
                    } else {
                        0.0
                    };
                    (i, parent, best_kd, edits, core, pct, band_req)
                })
                .collect()
        });
        // Headroom: among clear error-copy links (<= CLEAR_ERROR_COPY_PCT
        // divergence) the max kdist is the ceiling the cutoff must cover.
        let (mut linked, mut total) = (0u64, 0u64);
        let mut ceiling = 0.0f64;
        let mut kds: Vec<f64> = Vec::with_capacity(rows.len());
        let mut band_ec: Vec<usize> = Vec::new(); // band_req of clear error-copy links
        for &(i, parent, kd, edits, core, pct, band_req) in &rows {
            total += 1;
            if kd < p.cutoff {
                linked += 1;
            }
            if pct <= CLEAR_ERROR_COPY_PCT {
                ceiling = ceiling.max(kd);
                band_ec.push(band_req);
            }
            kds.push(kd);
            let ratio = s.counts[parent] as f64 / s.counts[i].max(1) as f64;
            writeln!(
                w,
                "{},{},{},{ratio:.2},{kd:.4},{edits},{core},{pct:.3},{band_req},{}",
                s.name,
                s.counts[i],
                s.counts[parent],
                (kd < p.cutoff) as u8
            )?;
        }
        if p.verbose && total > 0 {
            kds.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p90 = kds[(kds.len() * 9 / 10).min(kds.len() - 1)];
            eprintln!(
                "[kdist] {} : {total} children | nearest-parent kdist median {:.3} p90 {:.3} | \
                 {linked} ({:.1}%) within cutoff {} | clear-error-copy ceiling {:.3} -> headroom {:.3}",
                s.name,
                kds[kds.len() / 2],
                p90,
                100.0 * linked as f64 / total as f64,
                p.cutoff,
                ceiling,
                p.cutoff - ceiling,
            );
            eprintln!(
                "[kdist] {} : {}",
                s.name,
                band_fit("clear-error-copy", &band_ec)
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Params, ScreenBackend, find_derep_for_sample, load_dada, load_dada_pooled};
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;
    use std::path::PathBuf;

    /// Self-cleaning scratch dir.
    struct TmpDir(PathBuf);
    impl TmpDir {
        fn new(tag: &str) -> Self {
            static SEQ: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
            let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let d = std::env::temp_dir()
                .join(format!("dada2rs_kdist_{tag}_{}_{n}", std::process::id()));
            let _ = std::fs::remove_dir_all(&d);
            std::fs::create_dir_all(&d).unwrap();
            TmpDir(d)
        }
    }
    impl Drop for TmpDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn touch(dir: &std::path::Path, name: &str) {
        std::fs::File::create(dir.join(name)).unwrap();
    }

    /// Write a plain-JSON derep with an explicit internal `sample` field.
    fn write_derep(dir: &std::path::Path, name: &str, sample: &str) {
        let mut f = std::fs::File::create(dir.join(name)).unwrap();
        write!(f, r#"{{"sample":"{sample}","uniques":[]}}"#).unwrap();
    }

    /// A gzipped `dada` output paired with its plain derep must load cleanly:
    /// `load_dada` has to gunzip a `.json.gz` input rather than parse the raw
    /// gzip bytes as JSON (regression for the "expected value, line 1 column 1"
    /// failure on gzipped dada outputs, e.g. those from `dada-pooled --gzip`).
    /// Minimal `Params` for the loader tests: only `k` and the screen fields are
    /// consulted on these paths.
    fn test_params(k: usize) -> Params {
        Params {
            k,
            screen_backend: ScreenBackend::Kmer,
            minimizer_k: crate::minimizers::MINIMIZER_K,
            minimizer_w: crate::minimizers::MINIMIZER_W,
            derive_cutoff: false,
            derive_uniform_pairs: false,
            cutoff: 0.42,
            leak_pct: 10.0,
            band: -1,
            max_pairs: 0,
            max_uniques: 0,
            per_sample: false,
            nearest_parent: false,
            from_dada: false,
            from_dada_pooled: false,
            derep_dir: None,
            threads: 1,
            seed: 1,
            output: None,
            verbose: false,
        }
    }

    #[test]
    fn load_dada_reads_gzipped_output() {
        let d = TmpDir::new("gzdada");
        // Matching derep: one unique so `map.len()` lines up.
        let mut df = std::fs::File::create(d.0.join("F3D0.json")).unwrap();
        write!(
            df,
            r#"{{"sample":"F3D0","sort_order":"abundance_desc",
                 "uniques":[{{"sequence":"ACGTACGT","count":5}}]}}"#
        )
        .unwrap();
        // Gzipped dada output.
        let dada_path = d.0.join("F3D0.dada.json.gz");
        let mut gz = GzEncoder::new(
            std::fs::File::create(&dada_path).unwrap(),
            Compression::default(),
        );
        write!(
            gz,
            r#"{{"sample":"F3D0",
                 "asvs":[{{"sequence":"ACGTACGT","abundance":5,
                           "birth_type":"Abundance","birth_pval":0.0}}],
                 "map":[0]}}"#
        )
        .unwrap();
        gz.finish().unwrap();

        let sample =
            load_dada(&dada_path, &d.0, &test_params(7)).expect("gzipped dada output should load");
        assert_eq!(sample.name, "F3D0");
        assert_eq!(sample.enc.len(), 1);
    }

    /// A `_pooled.json` record loads as one self-contained population: uniques,
    /// global map, and centers come from the single file (no derep_dir), and a
    /// `null` map entry becomes a `None` fate (a globally-failed unique).
    #[test]
    fn load_pooled_record_self_contained() {
        let d = TmpDir::new("pooled");
        let path = d.0.join("_pooled.json");
        let mut f = std::fs::File::create(&path).unwrap();
        write!(
            f,
            r#"{{"tag":"dada-pooled-record","num_uniques":3,"num_asvs":1,
                 "uniques":[{{"sequence":"ACGTACGT","count":10}},
                            {{"sequence":"ACGTACGA","count":3}},
                            {{"sequence":"TTTTTTTT","count":1}}],
                 "map":[0,0,null],
                 "asvs":[{{"sequence":"ACGTACGT","abundance":13,
                           "birth_type":"Initial","birth_pval":0.0}}]}}"#
        )
        .unwrap();

        let s = load_dada_pooled(&path, &test_params(5)).expect("pooled record should load");
        assert_eq!(s.name, "__pooled__");
        assert_eq!(s.enc.len(), 3);
        assert_eq!(s.counts, vec![10, 3, 1]);
        assert_eq!(s.map, vec![Some(0), Some(0), None]); // third unique globally failed
        assert_eq!(s.c_enc.len(), 1); // one center, no derep_dir consulted
        assert_eq!(s.c_ab, vec![13]);
    }

    #[test]
    fn exact_match_preferred() {
        let d = TmpDir::new("exact");
        touch(&d.0, "F3D0.json");
        touch(&d.0, "F3D0.derep.R1.json.gz"); // also present; exact wins
        let got = find_derep_for_sample(&d.0, "F3D0").unwrap();
        assert_eq!(got.file_name().unwrap(), "F3D0.json");
    }

    #[test]
    fn infix_and_gzip_name_resolves() {
        // The issue #72 case: only a renamed, gzipped derep exists.
        let d = TmpDir::new("infix");
        touch(&d.0, "F3D0_S188.derep.R1.json.gz");
        let got = find_derep_for_sample(&d.0, "F3D0_S188").unwrap();
        assert_eq!(got.file_name().unwrap(), "F3D0_S188.derep.R1.json.gz");
    }

    #[test]
    fn prefix_component_is_exact_not_substring() {
        // `S188` must not match an `S1888` file.
        let d = TmpDir::new("prefix");
        touch(&d.0, "F3D0_S1888.derep.R1.json.gz");
        assert!(find_derep_for_sample(&d.0, "F3D0_S188").is_err());
    }

    #[test]
    fn missing_is_not_found() {
        let d = TmpDir::new("missing");
        touch(&d.0, "other.json");
        let e = find_derep_for_sample(&d.0, "F3D0").unwrap_err();
        assert_eq!(e.kind(), std::io::ErrorKind::NotFound);
    }

    #[test]
    fn multiple_candidates_disambiguated_by_sample_field() {
        // R1 and R2 share the `F3D0` prefix in one dir; the internal `sample`
        // field picks the right one.
        let d = TmpDir::new("multi");
        write_derep(&d.0, "F3D0.derep.R1.json", "F3D0");
        write_derep(&d.0, "F3D0.derep.R2.json", "F3D0_R2");
        let got = find_derep_for_sample(&d.0, "F3D0").unwrap();
        assert_eq!(got.file_name().unwrap(), "F3D0.derep.R1.json");
    }
}
