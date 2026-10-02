//! Output must be byte-identical across thread counts and `--sample-jobs`.

mod common;

use std::path::Path;

use common::{fixture, run, scratch, shared_err_model};

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
