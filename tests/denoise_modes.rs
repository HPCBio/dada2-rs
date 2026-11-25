//! Equivalences between the denoising entry points.
//!
//! `dada-pseudo` must reproduce the manual four-step pseudo-pooling recipe
//! (round-1 `dada` → `make-sequence-table` → `seq-table-to-fasta --prevalence 2`
//! → round-2 `dada --prior`), and multi-input `dada` must match one single-input
//! run per file. Everything runs with `--threads 1` for determinism.

mod common;

use std::path::Path;

use common::{
    BIN, assert_same_outputs, asv_set, err_model, fasta_seqs, fixture, normalized, run,
    run_expect_err, scratch,
};

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
    rates += "x\tins_q_novel\t40\t10000000\t10\nx\tins_q_run\t40\t10000000\t10\n";
    rates += "x\tdel_fq_L1\t40\t10000000\t10\nx\tdel_fq_L4\t40\t100000\t1\n";
    rates += "x\tins_len\t1\t0\t100\nx\tdel_len\t1\t0\t95\nx\tdel_len\t3\t0\t5\n";
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
    assert!(
        rates.contains("# indel_fit: pooled\ntable\tkey\trate\n"),
        "{rates}"
    );
    assert!(rates.lines().any(|l| l.starts_with("del_fq_L1\t")));
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

/// `learn-errors --indel-fit-cmd` hands each iteration's counts to a command
/// and reads its fitted rates back, as `--errfun-cmd` does for substitutions;
/// a failing command's stderr reaches the user.
#[test]
fn learn_errors_fits_indel_rates_with_an_external_command() {
    let dir = scratch("learn_indels_external");
    let f1 = fixture("sam1F.fastq.gz");
    let f1 = f1.to_str().unwrap();
    let p = |n: &str| dir.join(n).to_str().unwrap().to_string();
    // A stand-in fit: flat rates, and a copy of the counts it was given.
    let mut rates = String::from("table\tkey\trate\n");
    for t in [
        "ins_q_novel",
        "ins_q_run",
        "del_fq_L1",
        "del_fq_L2",
        "del_fq_L3",
        "del_fq_L4",
        "del_fq_L5",
    ] {
        rates += &format!("{t}\t0\t1e-4\n");
    }
    rates += "ins_r\t0\t0.1\ndel_r\t0\t0.1\n";
    std::fs::write(p("flat.tsv"), &rates).unwrap();
    std::fs::write(
        p("fit.sh"),
        format!("cp \"$1\" {}\ncp {} \"$2\"\n", p("seen.tsv"), p("flat.tsv")),
    )
    .unwrap();
    let fit = format!("sh {}", p("fit.sh"));
    run(&[
        "learn-errors",
        f1,
        "-o",
        &p("err.json"),
        "--indel-rates-out",
        &p("rates.tsv"),
        "--indel-max-iter",
        "2",
        "--indel-fit-cmd",
        &fit,
    ]);
    let seen = std::fs::read_to_string(p("seen.tsv")).unwrap();
    assert!(
        seen.lines().any(|l| l.starts_with("learned\tdel_fq_L1\t")),
        "{seen}"
    );
    let out = std::fs::read_to_string(p("rates.tsv")).unwrap();
    assert!(out.contains(&format!("# indel_fit: {fit}\n")), "{out}");
    assert!(out.lines().any(|l| l == "ins_r\t0\t1e-1"), "{out}");

    std::fs::write(p("fail.sh"), "echo 'mgcv is missing' >&2\nexit 3\n").unwrap();
    let o = std::process::Command::new(BIN)
        .args([
            "learn-errors",
            f1,
            "-o",
            &p("err2.json"),
            "--indel-rates-out",
            &p("r2.tsv"),
        ])
        .args(["--indel-fit-cmd", &format!("sh {}", p("fail.sh"))])
        .output()
        .unwrap();
    assert!(!o.status.success());
    let e = String::from_utf8_lossy(&o.stderr);
    assert!(e.contains("mgcv is missing"), "{e}");
}

/// ASV lineage (#306): every ASV carries its `cluster` id, `birth_from` (null
/// only for the initial cluster) and its nearest other ASV, and `dada-pooled`
/// numbers clusters pool-wide, so a sequence keeps its id across samples.
#[test]
fn dada_outputs_carry_asv_lineage() {
    let dir = scratch("asv_lineage");
    let (f1, f2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let (f1, f2) = (f1.to_str().unwrap(), f2.to_str().unwrap());
    let err = err_model();
    let err = err.to_str().unwrap();
    let p = |n: &str| dir.join(n).to_str().unwrap().to_string();
    run(&["dada", f1, "--error-model", err, "-o", &p("one.json")]);
    run(&[
        "dada-pooled",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("pooled"),
    ]);
    run(&[
        "dada-pseudo",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("pseudo"),
    ]);

    let asvs = |path: &str| -> Vec<serde_json::Value> {
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        v["asvs"].as_array().unwrap().clone()
    };
    let check = |path: &str, pool_wide: bool| {
        let a = asvs(path);
        assert!(a.len() > 1, "{path}: need several ASVs");
        let ids: std::collections::BTreeSet<u64> =
            a.iter().map(|x| x["cluster"].as_u64().unwrap()).collect();
        for x in &a {
            let id = x["cluster"].as_u64().unwrap();
            let initial = x["birth_type"] == "Initial";
            assert_eq!(x["birth_from"].is_null(), initial, "{path}: {x}");
            let near = x["nearest"].as_u64().expect("nearest present");
            assert_ne!(near, id, "{path}: nearest is itself");
            if !pool_wide {
                assert!(
                    ids.contains(&near),
                    "{path}: nearest {near} not in this file"
                );
            }
            assert!(
                x["nearest_hamming"].as_u64().unwrap() + x["nearest_indels"].as_u64().unwrap() > 0
            );
        }
    };
    check(&p("one.json"), false);
    for d in ["pooled", "pseudo"] {
        check(&p(&format!("{d}/sam1F.json")), d == "pooled");
        check(&p(&format!("{d}/sam2F.json")), d == "pooled");
    }
    let id_of = |path: &str| -> std::collections::BTreeMap<String, u64> {
        asvs(path)
            .iter()
            .map(|x| {
                (
                    x["sequence"].as_str().unwrap().to_string(),
                    x["cluster"].as_u64().unwrap(),
                )
            })
            .collect()
    };
    let (s1, s2) = (
        id_of(&p("pooled/sam1F.json")),
        id_of(&p("pooled/sam2F.json")),
    );
    let shared: Vec<_> = s1.keys().filter(|s| s2.contains_key(*s)).collect();
    assert!(!shared.is_empty());
    for s in shared {
        assert_eq!(s1[s], s2[s], "pooled cluster ids differ across samples");
    }
}

/// The dada `params` name the alignment model off R's defaults, and only then,
/// so R-parity output keeps its shape (#313).
#[test]
fn dada_params_record_the_alignment_model() {
    let dir = scratch("align_params");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let params = |name: &str, extra: &[&str]| -> serde_json::Value {
        let out = dir.join(format!("{name}.json"));
        let base = [
            "dada",
            s1.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--threads",
            "1",
        ];
        run(&[&base[..], extra].concat());
        let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&out).unwrap()).unwrap();
        v["params"].clone()
    };
    let default = params("default", &[]);
    for key in ["end_model", "align_backend", "wfa_max_edits", "gap_open_p"] {
        assert!(default.get(key).is_none(), "default output gained {key}");
    }
    let affine = params("affine", &["--end-model", "global", "--gap-open-p", "1"]);
    assert_eq!(affine["end_model"], "global");
    assert_eq!(affine["gap_open_p"], -1);
    assert!(affine.get("align_backend").is_none());

    #[cfg(feature = "wfa")]
    {
        let wfa = params("wfa", &["--align-backend", "wfa2", "--wfa-max-edits", "0"]);
        assert_eq!(wfa["align_backend"], "wfa2");
        assert_eq!(wfa["wfa_max_edits"], 0);
    }
}

/// `--partition` is opt-in (issue #323): naming the default `divisive` changes
/// nothing on `dada` or `dada-pooled`, and records nothing. `--birth-rule`
/// and `--reassign` need `single-pass`, `--reassign` does not combine with
/// the UNOISE rule, both birth rules are recorded and thread-independent,
/// and `dada-pseudo` does not take the flag.
#[test]
fn partition_default_unchanged_and_guarded() {
    let dir = scratch("partition");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    let (s1, s2, err) = (
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        err.to_str().unwrap(),
    );

    let single_default = dir.join("single_default.json");
    let single_named = dir.join("single_named.json");
    for (out, extra) in [(&single_default, None), (&single_named, Some("divisive"))] {
        let mut args = vec![
            "dada",
            s1,
            "--error-model",
            err,
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ];
        if let Some(p) = extra {
            args.extend_from_slice(&["--partition", p]);
        }
        run(&args);
    }
    assert_eq!(normalized(&single_default), normalized(&single_named));
    assert!(!normalized(&single_named).contains("\"partition\""));

    let pooled_default = dir.join("pooled_default");
    let pooled_named = dir.join("pooled_named");
    for (out, extra) in [(&pooled_default, None), (&pooled_named, Some("divisive"))] {
        let mut args = vec![
            "dada-pooled",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            "1",
            "--output-dir",
            out.to_str().unwrap(),
        ];
        if let Some(p) = extra {
            args.extend_from_slice(&["--partition", p]);
        }
        run(&args);
    }
    assert_same_outputs(
        &pooled_default,
        &pooled_named,
        "dada-pooled --partition divisive",
    );

    let out = dir.join("rejected.json");
    let out = out.to_str().unwrap();
    let base = ["dada", s1, "--error-model", err, "-o", out];
    let stderr = run_expect_err(&[&base[..], &["--birth-rule", "unoise"]].concat());
    assert!(
        stderr.contains("--birth-rule applies only with --partition single-pass"),
        "{stderr}"
    );
    let stderr = run_expect_err(&[&base[..], &["--reassign"]].concat());
    assert!(
        stderr.contains("--reassign applies only with --partition single-pass"),
        "{stderr}"
    );
    let stderr = run_expect_err(
        &[
            &base[..],
            &[
                "--partition",
                "single-pass",
                "--birth-rule",
                "unoise",
                "--reassign",
            ],
        ]
        .concat(),
    );
    assert!(
        stderr.contains("does not apply with --birth-rule unoise"),
        "{stderr}"
    );

    // single-pass runs, is recorded, and does not depend on the thread count.
    let sp1 = dir.join("single_pass_t1");
    let sp4 = dir.join("single_pass_t4");
    for (out, threads) in [(&sp1, "1"), (&sp4, "4")] {
        run(&[
            "dada-pooled",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            threads,
            "--output-dir",
            out.to_str().unwrap(),
            "--partition",
            "single-pass",
        ]);
    }
    assert_same_outputs(&sp1, &sp4, "single-pass, 1 vs 4 threads");
    let any = sp1.join(&common::json_outputs(&sp1)[0]);
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&any).unwrap()).unwrap();
    assert_eq!(v["params"]["partition"], "single-pass");
    assert_eq!(v["params"]["birth_rule"], "pvalue");

    // The UNOISE rule runs, is recorded, and does not depend on the thread count.
    let un1 = dir.join("unoise_t1");
    let un4 = dir.join("unoise_t4");
    for (out, threads) in [(&un1, "1"), (&un4, "4")] {
        run(&[
            "dada-pooled",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            threads,
            "--output-dir",
            out.to_str().unwrap(),
            "--partition",
            "single-pass",
            "--birth-rule",
            "unoise",
            "--min-abund",
            "8",
        ]);
    }
    assert_same_outputs(&un1, &un4, "unoise, 1 vs 4 threads");
    let any = un1.join(&common::json_outputs(&un1)[0]);
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&any).unwrap()).unwrap();
    assert_eq!(v["params"]["birth_rule"], "unoise");
    let stderr = run_expect_err(&[
        "dada-pseudo",
        s1,
        "--error-model",
        err,
        "-o",
        out,
        "--partition",
        "divisive",
    ]);
    assert!(stderr.contains("--partition"), "{stderr}");
}

/// `--metrics-json` attributes the UNOISE compare phase as it does the other
/// partitions (#328): the per-call timings are summed, so almost none of the
/// phase is left unattributed, and the aligned count matches the run shape.
#[test]
fn unoise_compare_is_attributed() {
    let dir = scratch("unoise_metrics");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    let metrics = dir.join("unoise.metrics.json");
    run(&[
        "dada-pooled",
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        "--error-model",
        err.to_str().unwrap(),
        "--threads",
        "4",
        "--output-dir",
        dir.join("out").to_str().unwrap(),
        "--partition",
        "single-pass",
        "--birth-rule",
        "unoise",
        "--metrics-json",
        metrics.to_str().unwrap(),
        "--metrics-attribution",
    ]);
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&metrics).unwrap()).unwrap();
    let r = &m["runs"][0];
    let c = &r["compare"];
    let (nalign, nshroud) = (
        r["run"]["nalign"].as_u64().unwrap(),
        r["run"]["nshroud"].as_u64().unwrap(),
    );
    assert!(nalign > nshroud, "no pairs aligned");
    assert_eq!(c["screened"].as_u64().unwrap(), nalign);
    assert_eq!(c["aligned"].as_u64().unwrap(), nalign - nshroud);
    let total = c["total"].as_f64().unwrap();
    let a = &c["attribution"];
    let unattributed = a["unattributed"].as_f64().unwrap();
    assert!(a["map"].as_f64().unwrap() > 0.0, "map unmeasured: {c}");
    assert!(
        unattributed <= 0.1 * total,
        "{unattributed} s of {total} s unattributed"
    );
    assert!(
        a["stored"].as_u64().unwrap() > 0,
        "no best-centre updates: {c}"
    );
    let split = &c["split"];
    assert!(
        split["screen"].as_f64().unwrap() > 0.0,
        "screen unmeasured: {c}"
    );
    assert!(
        split["dp_kernel"].as_f64().unwrap() > 0.0,
        "DP unmeasured: {c}"
    );
}

/// `--birth-rule unoise` runs without `--error-model` (#349). The rule decides
/// nothing from the model, so the ASVs, abundances and `map` match a run with
/// one, and the values only a model can give are `null`. Every other rule, and
/// the flags that read the model, still require it.
#[test]
fn unoise_runs_without_error_model() {
    let dir = scratch("unoise_no_model");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let (s1, err) = (s1.to_str().unwrap(), err.to_str().unwrap());
    let p = |n: &str| dir.join(n).to_str().unwrap().to_string();
    let unoise = [
        "dada",
        s1,
        "--partition",
        "single-pass",
        "--birth-rule",
        "unoise",
    ];
    let outputs = |tag: &str| -> Vec<String> {
        [
            "--aux-outputs",
            "--indel-events",
            "{t}.indels.tsv",
            "--cluster-trace",
            "{t}.trace.json",
            "-o",
            "{t}.json",
        ]
        .iter()
        .map(|a| match a.strip_prefix("{t}") {
            Some(rest) => p(&format!("{tag}{rest}")),
            None => a.to_string(),
        })
        .collect()
    };
    let (mo, no) = (outputs("m"), outputs("n"));
    let mo: Vec<&str> = mo.iter().map(String::as_str).collect();
    let no: Vec<&str> = no.iter().map(String::as_str).collect();
    run(&[&unoise[..], &["--error-model", err], &mo].concat());
    run(&[&unoise[..], &no].concat());
    let read = |n: &str| -> serde_json::Value {
        serde_json::from_str(&std::fs::read_to_string(p(n)).unwrap()).unwrap()
    };
    let (m, n) = (read("m.json"), read("n.json"));
    assert_eq!(m["asvs"], n["asvs"], "ASVs differ without a model");
    assert_eq!(m["map"], n["map"], "map differs without a model");
    let mut ms = m["aux"].clone();
    let mut ns = n["aux"].clone();
    for (with, without) in ms["cluster_stats"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .zip(ns["cluster_stats"].as_array_mut().unwrap())
    {
        assert!(with["pval"].is_f64(), "cluster pval missing with a model");
        assert!(
            without["pval"].is_null(),
            "cluster pval without a model: {without}"
        );
        with.as_object_mut().unwrap().remove("pval");
        without.as_object_mut().unwrap().remove("pval");
    }
    assert_eq!(ms, ns, "aux outputs differ beyond the cluster p-values");
    let (mt, nt) = (read("m.trace.json"), read("n.trace.json"));
    assert!(mt["err_in"].is_array() && nt.get("err_in").is_none());
    let members = |t: &serde_json::Value| {
        t["clusters"]
            .as_array()
            .unwrap()
            .iter()
            .flat_map(|c| c["members"].as_array().unwrap().clone())
            .collect::<Vec<_>>()
    };
    let nm = members(&nt);
    assert_eq!(members(&mt).len(), nm.len());
    for f in ["lambda", "e_reads", "pval"] {
        assert!(
            members(&mt).iter().all(|x| x[f].is_f64()),
            "{f} missing with a model"
        );
        assert!(nm.iter().all(|x| x[f].is_null()), "{f} without a model");
    }
    // The indel table's `p` and `lambda` columns (7 and 8) are `NA`.
    let columns = |n: &str| -> Vec<(String, String)> {
        std::fs::read_to_string(p(n))
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with("sample\t"))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[6].to_string(), f[7].to_string())
            })
            .collect()
    };
    let (mi, ni) = (columns("m.indels.tsv"), columns("n.indels.tsv"));
    assert!(!ni.is_empty() && mi.len() == ni.len(), "indel rows differ");
    assert!(mi.iter().all(|(p, l)| p != "NA" && l != "NA"));
    assert!(
        ni.iter().all(|(p, l)| p == "NA" && l == "NA"),
        "{:?}",
        &ni[..1]
    );

    let stderr = run_expect_err(&["dada", s1, "-o", &p("x.json")]);
    assert!(
        stderr.contains("--error-model is required unless --birth-rule unoise"),
        "{stderr}"
    );
    for flag in [
        &["--inherit-err-params"][..],
        &["--use-err-in"],
        &["--indel-rates", err],
    ] {
        let stderr = run_expect_err(&[&unoise[..], flag, &["-o", &p("x.json")]].concat());
        assert!(
            stderr.contains(&format!("{} requires --error-model", flag[0])),
            "{stderr}"
        );
    }
}

/// Under `--birth-rule unoise` an omitted `--min-abund` is 8, vsearch
/// `--cluster_unoise`'s `--minsize` default, even with `--inherit-err-params`
/// (the error model carries R's 1). An explicit value still wins, and the
/// p-value rule keeps R's default of 1.
#[test]
fn unoise_min_abund_defaults_to_eight() {
    let dir = scratch("unoise_min_abund");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    let (s1, s2, err) = (
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        err.to_str().unwrap(),
    );
    let min_abund = |name: &str, extra: &[&str]| {
        let out = dir.join(name);
        let mut args = vec![
            "dada-pooled",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            "1",
            "--output-dir",
            out.to_str().unwrap(),
            "--partition",
            "single-pass",
        ];
        args.extend_from_slice(extra);
        run(&args);
        let any = out.join(&common::json_outputs(&out)[0]);
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&any).unwrap()).unwrap();
        (out, v["params"]["min_abund"].as_u64().unwrap())
    };

    let (omitted, n) = min_abund("omitted", &["--birth-rule", "unoise"]);
    assert_eq!(n, 8);
    let (explicit, n) = min_abund("explicit8", &["--birth-rule", "unoise", "--min-abund", "8"]);
    assert_eq!(n, 8);
    assert_same_outputs(
        &omitted,
        &explicit,
        "unoise, omitted vs explicit --min-abund 8",
    );
    let (_, n) = min_abund(
        "inherit",
        &["--birth-rule", "unoise", "--inherit-err-params"],
    );
    assert_eq!(n, 8);
    let (_, n) = min_abund("explicit1", &["--birth-rule", "unoise", "--min-abund", "1"]);
    assert_eq!(n, 1);
    let (_, n) = min_abund("pvalue", &[]);
    assert_eq!(n, 1);
}

/// `--batch-align` is on by default (#337) and changes no output: every mode
/// matches `--batch-align false`, and bare `--batch-align` still parses. Runs
/// at 4 threads, because batching applies only to the parallel compare.
#[test]
fn batch_align_default_matches_per_pair() {
    let dir = scratch("batch_default");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s2 = fixture("sam2F.fastq.gz");
    let (s1, s2, err) = (
        s1.to_str().unwrap(),
        s2.to_str().unwrap(),
        err.to_str().unwrap(),
    );
    let p = |n: &str| dir.join(n).to_str().unwrap().to_string();

    // The default run must actually batch, or the comparisons below prove nothing.
    run(&[
        "dada-pooled",
        s1,
        s2,
        "--error-model",
        err,
        "--threads",
        "4",
        "--output-dir",
        &p("pooled_default"),
        "--metrics-json",
        &p("pooled.metrics.json"),
    ]);
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(p("pooled.metrics.json")).unwrap()).unwrap();
    let deferred = m["runs"][0]["compare"]["batch"]["deferred"]
        .as_u64()
        .unwrap();
    assert!(deferred > 0, "default run deferred no pairs");
    run(&[
        "dada-pooled",
        s1,
        s2,
        "--error-model",
        err,
        "--threads",
        "4",
        "--output-dir",
        &p("pooled_off"),
        "--batch-align",
        "false",
    ]);
    assert_same_outputs(
        &dir.join("pooled_default"),
        &dir.join("pooled_off"),
        "dada-pooled",
    );

    // The UNOISE rule's comparison batches too (#345), and must defer pairs
    // here for the comparison to test anything.
    let unoise = |name: &str, extra: &[&str]| {
        let base = [
            "dada-pooled",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            "4",
            "--partition",
            "single-pass",
            "--birth-rule",
            "unoise",
            "--min-abund",
            "2",
            "--output-dir",
        ];
        run(&[&base[..], &[p(name).as_str()], extra].concat());
    };
    let unoise_metrics = p("unoise.metrics.json");
    unoise("unoise_default", &["--metrics-json", &unoise_metrics]);
    unoise("unoise_off", &["--batch-align", "false"]);
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&unoise_metrics).unwrap()).unwrap();
    let deferred = m["runs"][0]["compare"]["batch"]["deferred"]
        .as_u64()
        .unwrap();
    assert!(deferred > 0, "default unoise run deferred no pairs");
    assert_same_outputs(
        &dir.join("unoise_default"),
        &dir.join("unoise_off"),
        "dada-pooled unoise",
    );

    for (name, flags) in [
        ("default", &[][..]),
        ("off", &["--batch-align", "false"][..]),
        ("bare", &["--batch-align"][..]),
    ] {
        let base = ["dada", s1, "--error-model", err, "--threads", "4", "-o"];
        run(&[
            &base[..],
            &[p(&format!("dada_{name}.json")).as_str()],
            flags,
        ]
        .concat());
        let base = [
            "dada-pseudo",
            s1,
            s2,
            "--error-model",
            err,
            "--threads",
            "4",
            "--output-dir",
        ];
        run(&[&base[..], &[p(&format!("pseudo_{name}")).as_str()], flags].concat());
        let base = ["learn-errors", s1, s2, "--threads", "4", "-o"];
        run(&[&base[..], &[p(&format!("le_{name}.json")).as_str()], flags].concat());
    }
    for other in ["off", "bare"] {
        assert_eq!(
            normalized(&dir.join("dada_default.json")),
            normalized(&dir.join(format!("dada_{other}.json"))),
            "dada, default vs {other}"
        );
        assert_same_outputs(
            &dir.join("pseudo_default"),
            &dir.join(format!("pseudo_{other}")),
            &format!("dada-pseudo, default vs {other}"),
        );
        assert_eq!(
            normalized(&dir.join("le_default.json")),
            normalized(&dir.join(format!("le_{other}.json"))),
            "learn-errors, default vs {other}"
        );
    }
}
