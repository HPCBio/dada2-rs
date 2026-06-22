//! End-to-end CLI tests for multi-input `dada` and `dada-pseudo`.
//!
//! These run the built binary against the committed fixtures in `data/dada2/`
//! and assert two equivalences that define the new subcommands:
//!
//! 1. `dada-pseudo` produces exactly the same per-sample ASVs as the manual
//!    four-step pseudo-pooling recipe (round-1 `dada` → `make-sequence-table`
//!    → `seq-table-to-fasta --prevalence 2` → round-2 `dada --prior`).
//! 2. Multi-input `dada` produces byte-identical per-sample output to running
//!    single-input `dada` once per file.
//!
//! Everything runs with `--threads 1` for determinism.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const BIN: &str = env!("CARGO_BIN_EXE_dada2-rs");

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixture(name: &str) -> PathBuf {
    // Fixtures live under tests/ (tracked); the repo's /data dir is gitignored
    // and so is absent on CI.
    manifest_dir().join("tests/fixtures").join(name)
}

/// Per-test scratch dir under the target tmp area; cleaned and recreated.
fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dada2rs_{}_{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run the binary; panic with stderr on a non-zero exit.
fn run(args: &[&str]) {
    let out = Command::new(BIN)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn {BIN}: {e}"));
    assert!(
        out.status.success(),
        "command failed: dada2-rs {}\n--- stderr ---\n{}",
        args.join(" "),
        String::from_utf8_lossy(&out.stderr),
    );
}

/// Run the binary expecting failure; return stderr.
fn run_expect_err(args: &[&str]) -> String {
    let out = Command::new(BIN).args(args).output().unwrap();
    assert!(
        !out.status.success(),
        "expected failure but command succeeded: dada2-rs {}",
        args.join(" "),
    );
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// An integer field from a dada output JSON's `params` block.
fn param_i64(path: &Path, key: &str) -> i64 {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    v["params"][key]
        .as_i64()
        .unwrap_or_else(|| panic!("no integer params.{key} in {}", path.display()))
}

/// Sorted set of (sequence, abundance) from a `dada`/`dada-pseudo` output JSON.
fn asv_set(path: &Path) -> BTreeSet<(String, i64)> {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    v["asvs"]
        .as_array()
        .unwrap_or_else(|| panic!("no asvs array in {}", path.display()))
        .iter()
        .map(|a| {
            (
                a["sequence"].as_str().unwrap().to_ascii_uppercase(),
                a["abundance"].as_i64().unwrap(),
            )
        })
        .collect()
}

/// Set of (uppercased) sequences in a FASTA file.
fn fasta_seqs(path: &Path) -> BTreeSet<String> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('>') && !l.trim().is_empty())
        .map(|l| l.trim().to_ascii_uppercase())
        .collect()
}

/// A loess error model learned once from the two committed forward fixtures and
/// shared across all tests in this binary (learning is the slow step, so doing
/// it once keeps CI fast). `OnceLock::get_or_init` runs the closure exactly
/// once even though tests execute on parallel threads.
fn shared_err_model() -> PathBuf {
    static ERR: OnceLock<PathBuf> = OnceLock::new();
    ERR.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("dada2rs_shared_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let err = dir.join("err.json");
        run(&[
            "learn-errors",
            fixture("sam1F.fastq.gz").to_str().unwrap(),
            fixture("sam2F.fastq.gz").to_str().unwrap(),
            "--errfun",
            "loess",
            "--threads",
            "1",
            "-o",
            err.to_str().unwrap(),
        ]);
        err
    })
    .clone()
}

#[test]
fn dada_pseudo_matches_manual_recipe() {
    let dir = scratch("pseudo");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    // --- dada-pseudo (one pass) ---
    let pseudo_out = dir.join("pseudo");
    let pseudo_priors = dir.join("priors_pseudo.fasta");
    run(&[
        "dada-pseudo",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--output-dir",
        pseudo_out.to_str().unwrap(),
        "--pseudo-prevalence",
        "2",
        "--priors-out",
        pseudo_priors.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    // --- manual recipe ---
    let r1 = dir.join("r1");
    let r2 = dir.join("r2");
    std::fs::create_dir_all(&r1).unwrap();
    std::fs::create_dir_all(&r2).unwrap();
    let r1_s1 = r1.join("sam1F.json");
    let r1_s2 = r1.join("sam2F.json");
    for (inp, out) in [(&s1, &r1_s1), (&s2, &r1_s2)] {
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ]);
    }
    let seqtab = dir.join("seqtab.json");
    run(&[
        "make-sequence-table",
        r1_s1.to_str().unwrap(),
        r1_s2.to_str().unwrap(),
        "-o",
        seqtab.to_str().unwrap(),
    ]);
    let manual_priors = dir.join("priors_manual.fasta");
    run(&[
        "seq-table-to-fasta",
        seqtab.to_str().unwrap(),
        "--prevalence",
        "2",
        "-o",
        manual_priors.to_str().unwrap(),
    ]);
    let r2_s1 = r2.join("sam1F.json");
    let r2_s2 = r2.join("sam2F.json");
    for (inp, out) in [(&s1, &r2_s1), (&s2, &r2_s2)] {
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--prior",
            manual_priors.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ]);
    }

    // --- priors must be the same set ---
    assert_eq!(
        fasta_seqs(&pseudo_priors),
        fasta_seqs(&manual_priors),
        "selected prior sequences differ between dada-pseudo and the manual recipe",
    );

    // --- per-sample ASVs must match ---
    for sample in ["sam1F.json", "sam2F.json"] {
        let pseudo = asv_set(&pseudo_out.join(sample));
        let manual = asv_set(&r2.join(sample));
        assert_eq!(
            pseudo, manual,
            "ASV (sequence, abundance) sets differ for {sample}",
        );
    }
}

/// dada-pseudo output must be accepted by the downstream consumers that key off
/// the `dada2_rs_command` tag (merge-pairs and make-sequence-table). This is the
/// path that broke in benchmarking: those readers allowlisted only "dada" /
/// "dada-pooled". The reverse error model is shared with the forward one here —
/// we're asserting the tag is *accepted*, not denoising correctness.
#[test]
fn dada_pseudo_output_feeds_downstream_steps() {
    let dir = scratch("pseudo_downstream");
    let err = shared_err_model();
    let f1 = fixture("sam1F.fastq.gz");
    let f2 = fixture("sam2F.fastq.gz");
    let r1 = fixture("sam1R.fastq.gz");
    let r2 = fixture("sam2R.fastq.gz");

    let fwd = dir.join("pseudo_fwd");
    let rev = dir.join("pseudo_rev");
    for (ins, out) in [([&f1, &f2], &fwd), ([&r1, &r2], &rev)] {
        run(&[
            "dada-pseudo",
            ins[0].to_str().unwrap(),
            ins[1].to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--threads",
            "1",
        ]);
    }

    // make-sequence-table directly on dada-pseudo per-sample JSONs.
    let seqtab = dir.join("seqtab.json");
    run(&[
        "make-sequence-table",
        fwd.join("sam1F.json").to_str().unwrap(),
        fwd.join("sam2F.json").to_str().unwrap(),
        "-o",
        seqtab.to_str().unwrap(),
    ]);

    // merge-pairs on dada-pseudo forward + reverse output.
    let merged = dir.join("merged.json");
    run(&[
        "merge-pairs",
        "--fwd-dada",
        fwd.join("sam1F.json").to_str().unwrap(),
        fwd.join("sam2F.json").to_str().unwrap(),
        "--rev-dada",
        rev.join("sam1R.json").to_str().unwrap(),
        rev.join("sam2R.json").to_str().unwrap(),
        "--fwd-fastq",
        f1.to_str().unwrap(),
        f2.to_str().unwrap(),
        "--rev-fastq",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "-o",
        merged.to_str().unwrap(),
    ]);
}

/// dada-pseudo denoises samples with bounded across-sample concurrency
/// (`--sample-jobs`). Per-sample `dada_uniques` is deterministic and round-1
/// prior selection is a set union, so output must be byte-identical regardless
/// of how many samples run concurrently (this also pins it to the serial path).
#[test]
fn dada_pseudo_is_deterministic_across_sample_jobs() {
    let dir = scratch("pseudo_jobs");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    let run_jobs = |jobs: &str, out: &Path| {
        run(&[
            "dada-pseudo",
            s1.to_str().unwrap(),
            s2.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--pseudo-prevalence",
            "2",
            "--threads",
            "4",
            "--sample-jobs",
            jobs,
        ]);
    };
    let (j1, j2) = (dir.join("j1"), dir.join("j2"));
    run_jobs("1", &j1);
    run_jobs("2", &j2);
    for sample in ["sam1F.json", "sam2F.json"] {
        assert_eq!(
            std::fs::read(j1.join(sample)).unwrap(),
            std::fs::read(j2.join(sample)).unwrap(),
            "dada-pseudo output for {sample} differs between --sample-jobs 1 and 2",
        );
    }
}

/// dada-pseudo streams by default (re-reading inputs per round); `--cache-samples`
/// holds all uniques in memory instead. The two must produce byte-identical
/// output — caching changes only *when* uniques are materialized, not the result.
#[test]
fn dada_pseudo_streaming_matches_cached() {
    let dir = scratch("pseudo_lowmem");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    let run_mode = |out: &Path, cache_samples: bool| {
        let mut args = vec![
            "dada-pseudo",
            s1.to_str().unwrap(),
            s2.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--pseudo-prevalence",
            "2",
            "--threads",
            "4",
        ];
        if cache_samples {
            args.push("--cache-samples");
        }
        run(&args);
    };
    // Default = streaming; --cache-samples = the all-in-memory mode.
    let (cache, stream) = (dir.join("cache"), dir.join("stream"));
    run_mode(&cache, true);
    run_mode(&stream, false);
    for sample in ["sam1F.json", "sam2F.json"] {
        assert_eq!(
            std::fs::read(cache.join(sample)).unwrap(),
            std::fs::read(stream.join(sample)).unwrap(),
            "dada-pseudo streaming (default) output for {sample} differs from --cache-samples",
        );
    }
}

/// dada-pooled loads/dereplicates samples concurrently (reassembled by input
/// index) and pools them into one inference; output must be byte-identical
/// regardless of thread count.
#[test]
fn dada_pooled_is_deterministic_across_threads() {
    let dir = scratch("pooled_det");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    let run_at = |t: &str, out: &Path| {
        run(&[
            "dada-pooled",
            s1.to_str().unwrap(),
            s2.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--threads",
            t,
        ]);
    };
    let (t1, t8) = (dir.join("t1"), dir.join("t8"));
    run_at("1", &t1);
    run_at("8", &t8);
    for sample in ["sam1F.json", "sam2F.json"] {
        assert_eq!(
            std::fs::read(t1.join(sample)).unwrap(),
            std::fs::read(t8.join(sample)).unwrap(),
            "dada-pooled output for {sample} differs between --threads 1 and 8",
        );
    }
}

/// merge-pairs parallelizes across samples; `collect` preserves input order, so
/// the output must be byte-identical regardless of thread count. (Same error
/// model is reused for both directions — this checks determinism, not biology.)
#[test]
fn merge_pairs_is_deterministic_across_threads() {
    let dir = scratch("merge_det");
    let err = shared_err_model();
    let f1 = fixture("sam1F.fastq.gz");
    let f2 = fixture("sam2F.fastq.gz");
    let r1 = fixture("sam1R.fastq.gz");
    let r2 = fixture("sam2R.fastq.gz");

    let dada = |inp: &Path, out: &Path| {
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ]);
    };
    let (f1j, f2j) = (dir.join("f1.json"), dir.join("f2.json"));
    let (r1j, r2j) = (dir.join("r1.json"), dir.join("r2.json"));
    dada(&f1, &f1j);
    dada(&f2, &f2j);
    dada(&r1, &r1j);
    dada(&r2, &r2j);

    let merge_at = |t: &str, out: &Path| {
        run(&[
            "merge-pairs",
            "--fwd-dada",
            f1j.to_str().unwrap(),
            f2j.to_str().unwrap(),
            "--rev-dada",
            r1j.to_str().unwrap(),
            r2j.to_str().unwrap(),
            "--fwd-fastq",
            f1.to_str().unwrap(),
            f2.to_str().unwrap(),
            "--rev-fastq",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--threads",
            t,
            "-o",
            out.to_str().unwrap(),
        ]);
    };
    let (m1, m4) = (dir.join("merged_t1.json"), dir.join("merged_t4.json"));
    merge_at("1", &m1);
    merge_at("4", &m4);
    assert_eq!(
        std::fs::read(&m1).unwrap(),
        std::fs::read(&m4).unwrap(),
        "merge-pairs output differs between --threads 1 and --threads 4",
    );
}

/// With an impossibly high --min-overlap every pair fails to merge. By default
/// those pairs are dropped; with --rescue-unmerged they are concatenated and
/// accepted (marked `concatenated: true`).
#[test]
fn merge_pairs_rescue_unmerged_concatenates() {
    let dir = scratch("merge_rescue");
    let err = shared_err_model();
    let f1 = fixture("sam1F.fastq.gz");
    let r1 = fixture("sam1R.fastq.gz");

    let dada = |inp: &Path, out: &Path| {
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ]);
    };
    let (f1j, r1j) = (dir.join("f1.json"), dir.join("r1.json"));
    dada(&f1, &f1j);
    dada(&r1, &r1j);

    let merge = |out: &Path, extra: &[&str]| {
        let mut args = vec![
            "merge-pairs",
            "--fwd-dada",
            f1j.to_str().unwrap(),
            "--rev-dada",
            r1j.to_str().unwrap(),
            "--fwd-fastq",
            f1.to_str().unwrap(),
            "--rev-fastq",
            r1.to_str().unwrap(),
            "--min-overlap",
            "5000",
            "-o",
            out.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        run(&args);
    };

    // Default: nothing merges, nothing rescued.
    let plain = dir.join("plain.json");
    merge(&plain, &[]);
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&plain).unwrap()).unwrap();
    let sample = &v["samples"][0];
    assert_eq!(sample["accepted_pairs"].as_u64().unwrap(), 0);
    assert_eq!(sample["merged"].as_array().unwrap().len(), 0);

    // Rescue: every distinct pair is concatenated and accepted.
    let rescued = dir.join("rescued.json");
    merge(&rescued, &["--rescue-unmerged", "--concat-nnn-len", "10"]);
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&rescued).unwrap()).unwrap();
    let sample = &v["samples"][0];
    assert!(sample["accepted_pairs"].as_u64().unwrap() > 0);
    let merged = sample["merged"].as_array().unwrap();
    assert!(!merged.is_empty());
    for m in merged {
        assert!(m["accept"].as_bool().unwrap());
        assert!(m["concatenated"].as_bool().unwrap());
        // Concatenated sequence carries the N spacer.
        assert!(m["sequence"].as_str().unwrap().contains("NNNNNNNNNN"));
    }
}

#[test]
fn dada_multi_input_matches_per_file_runs() {
    let dir = scratch("multi");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    // Multi-input run (with across-sample concurrency) -> per-sample files in a
    // directory. Asserting this equals the per-file single runs covers both that
    // multi-input matches single-input AND that --sample-jobs concurrency is
    // deterministic/correct.
    let multi = dir.join("multi");
    run(&[
        "dada",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--output-dir",
        multi.to_str().unwrap(),
        "--threads",
        "4",
        "--sample-jobs",
        "2",
    ]);

    // Single-input runs, one per file.
    for (inp, name) in [(&s1, "sam1F.json"), (&s2, "sam2F.json")] {
        let single = dir.join(format!("single_{name}"));
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            single.to_str().unwrap(),
        ]);
        assert_eq!(
            std::fs::read(multi.join(name)).unwrap(),
            std::fs::read(&single).unwrap(),
            "multi-input output for {name} differs from the single-input run",
        );
    }
}

/// `--homo-gap-p`, when unset, must default to `--gap-p` (R's
/// `HOMOPOLYMER_GAP_PENALTY = NULL` semantics); both are recorded in the output
/// `params` block. With neither flag the defaults are -8/-8 (unchanged).
#[test]
fn dada_homo_gap_defaults_to_gap_penalty() {
    let dir = scratch("gap_penalty");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");

    let run_dada = |out: &Path, extra: &[&str]| {
        let mut args = vec![
            "dada",
            s1.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        run(&args);
    };

    // Default: both -8.
    let def = dir.join("def.json");
    run_dada(&def, &[]);
    assert_eq!(param_i64(&def, "gap_p"), -8);
    assert_eq!(param_i64(&def, "homo_gap_p"), -8);

    // --gap-p set, --homo-gap-p unset: homo falls back to gap.
    let g = dir.join("g.json");
    run_dada(&g, &["--gap-p", "-4"]);
    assert_eq!(param_i64(&g, "gap_p"), -4);
    assert_eq!(param_i64(&g, "homo_gap_p"), -4);

    // Both set: independent.
    let gh = dir.join("gh.json");
    run_dada(&gh, &["--gap-p", "-4", "--homo-gap-p", "-1"]);
    assert_eq!(param_i64(&gh, "gap_p"), -4);
    assert_eq!(param_i64(&gh, "homo_gap_p"), -1);

    // Positive penalties are normalized to negative (R dada.R:223-227): a
    // positive --gap-p flips sign and homo falls back to the normalized value;
    // a positive --homo-gap-p flips independently.
    let pos = dir.join("pos.json");
    run_dada(&pos, &["--gap-p", "8"]);
    assert_eq!(param_i64(&pos, "gap_p"), -8);
    assert_eq!(param_i64(&pos, "homo_gap_p"), -8);

    let posh = dir.join("posh.json");
    run_dada(&posh, &["--gap-p", "-4", "--homo-gap-p", "1"]);
    assert_eq!(param_i64(&posh, "gap_p"), -4);
    assert_eq!(param_i64(&posh, "homo_gap_p"), -1);
}

#[test]
fn dada_input_output_guards() {
    let dir = scratch("guards");
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    // No error model needed: these must fail during argument validation,
    // before any denoising. We point --error-model at a path that exists so
    // the guard (not a missing-file error) is what trips.
    let err = shared_err_model();

    // >1 input with -o is rejected.
    let e = run_expect_err(&[
        "dada",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        dir.join("x.json").to_str().unwrap(),
    ]);
    assert!(e.contains("--output"), "unexpected error: {e}");

    // Single input with --output-dir is rejected.
    let e = run_expect_err(&[
        "dada",
        s1.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--output-dir",
        dir.join("d").to_str().unwrap(),
    ]);
    assert!(e.contains("--output-dir"), "unexpected error: {e}");
}

/// Dereplication orders uniques by abundance descending, ties broken lexically
/// by sequence — matching R `derepFastq` (its `qtables2` builds uniques in
/// lexical order, then a stable abundance sort preserves it among ties). This
/// is the order the DADA traversal assumes, so it must be reproduced exactly.
#[test]
fn derep_orders_by_abundance_then_lexical() {
    let dir = scratch("derep_order");
    // Two equal-abundance uniques (AAAA=2, CCCC=2) emitted in NON-lexical
    // first-seen order (CCCC before AAAA); GGGG=3 is the unique max. Expect
    // abundance-desc then lexical: GGGG, AAAA, CCCC — NOT the old first-seen
    // tie-break (which would give GGGG, CCCC, AAAA).
    let read = |id: &str, seq: &str| format!("@{id}\n{seq}\n+\n{}\n", "I".repeat(seq.len()));
    let mut fq = String::new();
    for id in ["c1", "c2"] {
        fq += &read(id, "CCCCCCCCCC");
    }
    for id in ["a1", "a2"] {
        fq += &read(id, "AAAAAAAAAA");
    }
    for id in ["g1", "g2", "g3"] {
        fq += &read(id, "GGGGGGGGGG");
    }
    let fq_path = dir.join("tie.fastq");
    std::fs::write(&fq_path, fq).unwrap();

    let out = dir.join("derep.json");
    run(&[
        "derep",
        fq_path.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);

    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&out).unwrap()).unwrap();
    let order: Vec<(String, i64)> = v["uniques"]
        .as_array()
        .expect("uniques array")
        .iter()
        .map(|u| {
            (
                u["sequence"].as_str().unwrap().to_string(),
                u["count"].as_i64().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        order,
        vec![
            ("GGGGGGGGGG".to_string(), 3),
            ("AAAAAAAAAA".to_string(), 2),
            ("CCCCCCCCCC".to_string(), 2),
        ],
        "derep uniques must be abundance-desc then lexical (R derepFastq order)",
    );
}

/// Denoising a FASTQ directly must produce byte-identical output to denoising
/// its pre-dereplicated JSON. This is the invariant that lets `dada*` consume a
/// derep JSON interchangeably with raw FASTQ — reading the JSON reconstructs
/// exactly the same uniques/counts/quals AND order as in-line dereplication
/// (e.g. streaming `dada-pseudo` round 2 re-reads derep JSON, not FASTQ).
#[test]
fn dada_from_fastq_matches_dada_from_derep_json() {
    let dir = scratch("derep_equiv");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");

    // dada directly from FASTQ
    let from_fastq = dir.join("from_fastq.json");
    run(&[
        "dada",
        s1.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        from_fastq.to_str().unwrap(),
    ]);

    // derep -> JSON, then dada from that JSON
    let derep_json = dir.join("derep.json");
    run(&[
        "derep",
        s1.to_str().unwrap(),
        "-o",
        derep_json.to_str().unwrap(),
    ]);
    let from_json = dir.join("from_json.json");
    run(&[
        "dada",
        derep_json.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        from_json.to_str().unwrap(),
    ]);

    // The `input_file` provenance field intentionally differs (the FASTQ name
    // vs the derep JSON name); strip it before comparing the denoising result.
    let strip_input_file = |path: &std::path::Path| {
        let mut v: serde_json::Value =
            serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        v.as_object_mut().unwrap().remove("input_file");
        v
    };
    assert_eq!(
        strip_input_file(&from_fastq),
        strip_input_file(&from_json),
        "dada output from derep JSON differs from dada output from FASTQ",
    );
}

/// `--failed-uniques` (issue #60) must emit a header plus exactly one row per
/// `map == null` unique (a unique that failed to denoise), with that unique's
/// sequence and in-sample abundance. The single-sample dada `map` is the clean
/// reference signal, so the TSV row count must equal the JSON null count and
/// every TSV sequence must be a `map == null` input unique.
#[test]
fn dada_failed_uniques_matches_map_nulls() {
    let dir = scratch("failed_uniques");
    let err = shared_err_model();
    let s1 = fixture("sam1F.fastq.gz");

    let out_json = dir.join("d.json");
    let fu_tsv = dir.join("failed.tsv");
    run(&[
        "dada",
        s1.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "-o",
        out_json.to_str().unwrap(),
        "--failed-uniques",
        fu_tsv.to_str().unwrap(),
    ]);

    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&out_json).unwrap()).unwrap();
    let map = v["map"].as_array().unwrap();
    let null_count = map.iter().filter(|m| m.is_null()).count();
    assert!(null_count > 0, "fixture should have some failed uniques");

    let tsv = std::fs::read_to_string(&fu_tsv).unwrap();
    let mut lines = tsv.lines();
    assert_eq!(
        lines.next().unwrap(),
        "sequence\tsample\treads",
        "TSV must start with the header row",
    );
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len(), null_count, "one TSV row per map==null unique",);
    for row in rows {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 3, "row must have sequence/sample/reads");
        assert_eq!(cols[1], "sam1F", "sample column");
        assert!(cols[2].parse::<u32>().unwrap() >= 1, "reads >= 1");
    }
}
