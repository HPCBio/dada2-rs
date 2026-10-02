//! Equivalences between the denoising entry points.
//!
//! `dada-pseudo` must reproduce the manual four-step pseudo-pooling recipe
//! (round-1 `dada` → `make-sequence-table` → `seq-table-to-fasta --prevalence 2`
//! → round-2 `dada --prior`), and multi-input `dada` must match one single-input
//! run per file. Everything runs with `--threads 1` for determinism.

mod common;

use std::path::Path;

use common::{asv_set, err_model, fasta_seqs, fixture, run, scratch};

#[test]
fn dada_pseudo_matches_manual_recipe() {
    let dir = scratch("pseudo");
    let err = err_model();
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

/// `--reestimate-err-between-rounds` must be strictly opt-in: with the flag off,
/// output is unchanged, and with it on the run records that its round 2 did not
/// use the supplied error model.
///
/// The flag exists to emulate R DADA2's `pool="pseudo"`, which re-fits `err`
/// from round 1's transitions inside the self-consistency loop (issue #100).
/// Whether that is the better behaviour is unresolved, so what is pinned here is
/// the *contract*: default behaviour untouched, and a re-estimated run is
/// self-describing in its output params. The fixtures are two shallow samples,
/// so this asserts the plumbing, not the magnitude of any difference.
#[test]
fn dada_pseudo_reestimate_err_is_opt_in_and_recorded() {
    let dir = scratch("pseudo_reestimate");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");

    let run_mode = |out: &std::path::Path, reestimate: bool| {
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
            "1",
        ];
        if reestimate {
            args.push("--reestimate-err-between-rounds");
        }
        run(&args);
    };

    let off = dir.join("off");
    let on = dir.join("on");
    run_mode(&off, false);
    run_mode(&on, true);

    for sample in ["sam1F.json", "sam2F.json"] {
        let off_json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(off.join(sample)).unwrap()).unwrap();
        let on_json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(on.join(sample)).unwrap()).unwrap();
        // The output envelope is flat, tagged by a `dada2_rs_command` key.
        let params = |v: &serde_json::Value| v["params"].clone();

        // Flag off: the provenance field is absent entirely, so pre-existing
        // outputs stay byte-identical to before the flag was added.
        assert!(
            params(&off_json)
                .get("reestimated_err_between_rounds")
                .is_none(),
            "{sample}: default run must not emit reestimated_err_between_rounds",
        );
        // Flag on: the run says so, because --error-model alone no longer
        // describes what round 2 actually used.
        assert_eq!(
            params(&on_json)["reestimated_err_between_rounds"],
            serde_json::json!(true),
            "{sample}: re-estimated run must record it in params",
        );
    }
}

/// dada-pseudo streams by default (re-reading inputs per round); `--cache-samples`
/// holds all uniques in memory instead. The two must produce byte-identical
/// output — caching changes only *when* uniques are materialized, not the result.
#[test]
fn dada_pseudo_streaming_matches_cached() {
    let dir = scratch("pseudo_lowmem");
    let err = err_model();
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

#[test]
fn dada_multi_input_matches_per_file_runs() {
    let dir = scratch("multi");
    let err = err_model();
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

/// Per-ASV sums of derep counts over a dada output's `map`, plus the number of
/// uniques the map leaves unassigned.
fn mapped_abundance(dada_json: &Path, derep_json: &Path) -> (Vec<i64>, usize) {
    let d: serde_json::Value = serde_json::from_slice(&std::fs::read(dada_json).unwrap()).unwrap();
    let u: serde_json::Value = serde_json::from_slice(&std::fs::read(derep_json).unwrap()).unwrap();
    let map = d["map"].as_array().expect("dada output has a map");
    let uniques = u["uniques"].as_array().unwrap();
    assert_eq!(
        map.len(),
        uniques.len(),
        "map must cover every derep unique"
    );
    let mut sums = vec![0i64; d["asvs"].as_array().unwrap().len()];
    let mut unmapped = 0;
    for (m, un) in map.iter().zip(uniques) {
        match m.as_u64() {
            Some(a) => sums[a as usize] += un["count"].as_i64().unwrap(),
            None => unmapped += 1,
        }
    }
    (sums, unmapped)
}

/// An ASV's reported abundance counts only the reads `map` assigns to it, as
/// R's `clustering$abundance` does (`b_make_clustering_df` sums correct raws
/// only). Reporting every member's reads overstated it by the reads that fail
/// `OMEGA_C` (#204). Checked for `dada` and `dada-pseudo`; the test requires an
/// unassigned unique, since without one both rules agree.
#[test]
fn asv_abundance_counts_only_mapped_reads() {
    let dir = scratch("mapped_abund");
    let err = err_model();
    let (s1, s2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let mut unmapped_total = 0;

    let per = dir.join("per");
    std::fs::create_dir_all(&per).unwrap();
    let pseudo = dir.join("pseudo");
    run(&[
        "dada-pseudo",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--output-dir",
        pseudo.to_str().unwrap(),
        "--threads",
        "1",
    ]);
    for (fq, name) in [(&s1, "sam1F"), (&s2, "sam2F")] {
        let derep = dir.join(format!("{name}.derep.json"));
        run(&["derep", fq.to_str().unwrap(), "-o", derep.to_str().unwrap()]);
        let dada = per.join(format!("{name}.json"));
        run(&[
            "dada",
            fq.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            dada.to_str().unwrap(),
        ]);
        for out in [dada, pseudo.join(format!("{name}.json"))] {
            let (sums, unmapped) = mapped_abundance(&out, &derep);
            let reported: Vec<i64> = asv_set_ordered(&out);
            assert_eq!(
                reported,
                sums,
                "{}: abundance != reads mapped to it",
                out.display()
            );
            unmapped_total += unmapped;
        }
    }
    assert!(
        unmapped_total > 0,
        "no unassigned uniques in the fixtures, so this test cannot distinguish the two rules"
    );
}

/// ASV abundances in output order.
fn asv_set_ordered(path: &Path) -> Vec<i64> {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
    v["asvs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| a["abundance"].as_i64().unwrap())
        .collect()
}
