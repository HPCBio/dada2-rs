//! Equivalences between the denoising entry points.
//!
//! `dada-pseudo` must reproduce the manual four-step pseudo-pooling recipe
//! (round-1 `dada` → `make-sequence-table` → `seq-table-to-fasta --prevalence 2`
//! → round-2 `dada --prior`), and multi-input `dada` must match one single-input
//! run per file. Everything runs with `--threads 1` for determinism.

mod common;

use std::path::Path;

use common::{BIN, asv_set, err_model, fasta_seqs, fixture, normalized, run, scratch};

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

/// `--pool-tiebreak` (#260): `first-seen`, the default, makes pool order follow
/// input order, as R's `combineDereps2` does; `lexical` makes it independent of
/// input order. Each output records the rule and its sample's input position,
/// and `--pooled-record` lists the inputs in order, so a run can be reproduced.
#[test]
fn pool_tiebreak_records_order_and_rule() {
    let dir = scratch("pool_tiebreak");
    let err = err_model();
    let (s1, s2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let pooled = |tag: &str, inputs: [&Path; 2], extra: &[&str]| -> serde_json::Value {
        let out = dir.join(tag);
        let rec = dir.join(format!("{tag}.pooled.json"));
        let mut args = vec![
            "dada-pooled",
            "--error-model",
            err.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--pooled-record",
            rec.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        args.extend(inputs.iter().map(|p| p.to_str().unwrap()));
        run(&args);
        serde_json::from_slice(&std::fs::read(rec).unwrap()).unwrap()
    };
    let order = |rec: &serde_json::Value| -> Vec<String> {
        rec["uniques"]
            .as_array()
            .unwrap()
            .iter()
            .map(|u| u["sequence"].as_str().unwrap().to_string())
            .collect()
    };

    let fwd = pooled("fs_fwd", [&s1, &s2], &[]);
    let rev = pooled("fs_rev", [&s2, &s1], &[]);
    assert_eq!(fwd["pool_tiebreak"], "first-seen");
    assert_ne!(
        order(&fwd),
        order(&rev),
        "first-seen must follow input order"
    );
    let inputs: Vec<&str> = fwd["inputs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|i| i["sample"].as_str().unwrap())
        .collect();
    assert_eq!(inputs, ["sam1F", "sam2F"]);

    let lex_fwd = pooled("lex_fwd", [&s1, &s2], &["--pool-tiebreak", "lexical"]);
    let lex_rev = pooled("lex_rev", [&s2, &s1], &["--pool-tiebreak", "lexical"]);
    assert_eq!(lex_fwd["pool_tiebreak"], "lexical");
    assert_eq!(
        order(&lex_fwd),
        order(&lex_rev),
        "lexical must not depend on input order"
    );

    // Per-sample outputs carry the rule and their position in the input order.
    for (sample, index) in [("sam2F", 0), ("sam1F", 1)] {
        let v: serde_json::Value = serde_json::from_slice(
            &std::fs::read(dir.join("fs_rev").join(format!("{sample}.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(v["params"]["pool_tiebreak"], "first-seen");
        assert_eq!(v["pool_input_index"], index);
    }
}

/// `--end-model` is opt-in (issue #289): naming the default changes nothing,
/// and a model the error model was not learned with is named in a warning.
#[test]
fn end_model_default_unchanged_and_mismatch_warned() {
    let dir = scratch("end_model");
    let err = err_model();
    let inp = fixture("sam1F.fastq.gz");
    let dada = |extra: &[&str], out: &Path| {
        let mut args = vec![
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        let o = std::process::Command::new(BIN)
            .args(&args)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8_lossy(&o.stderr).into_owned()
    };

    let default = dir.join("default.json");
    let named = dir.join("named.json");
    let global = dir.join("global.json");
    let quiet = dada(&[], &default);
    dada(&["--end-model", "ends-free"], &named);
    assert_eq!(normalized(&default), normalized(&named));
    assert!(!quiet.contains("end_model"), "{quiet}");

    let warned = dada(&["--end-model", "global"], &global);
    assert!(
        warned.contains("end_model = Global (err model: EndsFree)"),
        "{warned}"
    );
}

/// The experimental indel-aware λ (#291) is off by default: an indel-only
/// variant is absorbed however abundant, as in R. With `--indel-rates` a 3-nt
/// deletion splits, while a homopolymer-length deletion stays absorbed under
/// the default ceiling and splits only with `--hp-ceiling off`.
#[test]
fn indel_lambda_is_opt_in_and_spares_homopolymer_changes() {
    let dir = scratch("indel_lambda");
    let err = err_model();
    let err = err.to_str().unwrap();
    // A centre with no long runs except one A4, and two indel-only variants.
    let unit = "ACGTTGCAGGCTAACCGTAGCATCGGATCCTAGGCTTAACGTCAGTCAGGATCCATGCAAGT";
    let centre = format!("{unit}{}AAAA{}", &unit[..40], &unit[10..60]);
    let del3 = format!("{}{}", &centre[..30], &centre[33..]);
    let a4 = centre.find("AAAA").unwrap();
    let hp_del = format!("{}{}", &centre[..a4], &centre[a4 + 1..]);
    assert!(
        !centre[..33].contains("AAAA"),
        "deleted bases must not be a run"
    );
    let mut fq = String::new();
    for (i, (s, n)) in [(&centre, 300), (&del3, 40), (&hp_del, 40)]
        .iter()
        .enumerate()
    {
        for k in 0..*n {
            fq += &format!("@r{i}_{k}\n{s}\n+\n{}\n", "I".repeat(s.len()));
        }
    }
    let fq_path = dir.join("in.fastq");
    std::fs::write(&fq_path, fq).unwrap();
    let fq_path = fq_path.to_str().unwrap();
    // Low indel rates, high enough deletion rate in runs of 4 to split with `off`.
    let mut rates = String::from("sample\ttable\tkey\topportunities\tevents\n");
    rates += "x\tins_q\t40\t10000000\t10\nx\tdel_hp\t1\t10000000\t10\nx\tdel_hp\t4\t100000\t1\n";
    rates += "x\tdel_fq\t40\t10000000\t10\nx\tins_len\t1\t0\t100\nx\tdel_len\t1\t0\t95\nx\tdel_len\t3\t0\t5\n";
    let rates_path = dir.join("rates.tsv");
    std::fs::write(&rates_path, rates).unwrap();
    let rates_path = rates_path.to_str().unwrap();

    let asvs = |name: &str, extra: &[&str]| -> (Vec<String>, serde_json::Value) {
        let out = dir.join(format!("{name}.json"));
        let base = [
            "dada",
            fq_path,
            "--error-model",
            err,
            "-o",
            out.to_str().unwrap(),
            "--threads",
            "1",
        ];
        run(&[&base[..], extra].concat());
        let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&out).unwrap()).unwrap();
        let seqs = v["asvs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| a["sequence"].as_str().unwrap().to_string())
            .collect();
        (seqs, v["params"].clone())
    };
    let (default, params) = asvs("default", &[]);
    assert_eq!(
        default,
        vec![centre.clone()],
        "R parity: indel-only variants are absorbed"
    );
    assert!(
        params.get("indel_lambda").is_none(),
        "default output must not change"
    );

    let (on, params) = asvs("on", &["--indel-rates", rates_path]);
    assert!(on.contains(&del3), "the 3-nt deletion splits: {on:?}");
    assert!(
        !on.contains(&hp_del),
        "the homopolymer change stays under the ceiling"
    );
    assert_eq!(params["indel_lambda"]["hp_ceiling"], 0.05);
    assert_eq!(params["indel_lambda"]["pair_window"], 5);

    let (off, _) = asvs("off", &["--indel-rates", rates_path, "--hp-ceiling", "off"]);
    assert!(
        off.contains(&hp_del),
        "with `off` it is scored at the run-of-4 rate: {off:?}"
    );
}

/// `learn-errors --indel-rates-out` (#291) learns indel rates after the
/// substitution model: the error-model JSON stays byte-identical (the
/// substitution model is still R's), the loop stops on a repeat, and `dada
/// --indel-rates` accepts the learned rates.
#[test]
fn learn_errors_learns_indel_rates_without_touching_the_error_model() {
    let dir = scratch("learn_indels");
    let (f1, f2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let (f1, f2) = (f1.to_str().unwrap(), f2.to_str().unwrap());
    let p = |n: &str| dir.join(n).to_str().unwrap().to_string();
    run(&["learn-errors", f1, f2, "-o", &p("plain.json")]);
    run(&[
        "learn-errors",
        f1,
        f2,
        "-o",
        &p("with.json"),
        "--indel-rates-out",
        &p("rates.tsv"),
        "--indel-history",
        &p("history.tsv"),
    ]);
    assert_eq!(
        normalized(&dir.join("plain.json")),
        normalized(&dir.join("with.json")),
        "indel learning changed the error model",
    );
    let history = std::fs::read_to_string(p("history.tsv")).unwrap();
    assert!(
        history.contains("# indel_converged: true"),
        "{}",
        &history[..400.min(history.len())]
    );
    let iterations: std::collections::BTreeSet<&str> = history
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("sample\t"))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert!(
        iterations.len() >= 2 && iterations.contains("iter0"),
        "{iterations:?}"
    );
    let rates = std::fs::read_to_string(p("rates.tsv")).unwrap();
    assert!(rates.lines().any(|l| l.starts_with("learned\tdel_fq\t")));
    run(&[
        "dada",
        f1,
        "--error-model",
        &p("with.json"),
        "--indel-rates",
        &p("rates.tsv"),
        "-o",
        &p("d.json"),
    ]);
}
