//! Hand-offs between pipeline steps: derep ordering, FASTQ vs derep-JSON input,
//! the failed-uniques TSV, merge-pairs, and downstream acceptance of each
//! command's output tag.

mod common;

use std::path::Path;

use common::{err_model, fixture, run, scratch};

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
    let err = err_model();
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
    let err = err_model();
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

/// An indel TSV without its `#` run-identity header.
fn tsv_body(path: impl AsRef<Path>) -> String {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect()
}

/// `--indel-events` / `--indel-opportunities` (issue #291) must not change the
/// dada JSON, and the `sub_q` table must reproduce `trans` from `--aux-outputs`:
/// both are read-weighted counts over the same final-subs pass.
#[test]
fn dada_indel_tables_leave_output_unchanged_and_match_trans() {
    let dir = scratch("indel_events");
    let err = err_model();
    let s1 = fixture("sam1F.fastq.gz");
    let s1 = s1.to_str().unwrap();
    let err = err.to_str().unwrap();

    let plain = dir.join("plain.json");
    run(&[
        "dada",
        s1,
        "--error-model",
        err,
        "--aux-outputs",
        "-o",
        plain.to_str().unwrap(),
    ]);
    let with = dir.join("with.json");
    let events = dir.join("events.tsv");
    let opp = dir.join("opp.tsv");
    run(&[
        "dada",
        s1,
        "--error-model",
        err,
        "--aux-outputs",
        "-o",
        with.to_str().unwrap(),
        "--indel-events",
        events.to_str().unwrap(),
        "--indel-opportunities",
        opp.to_str().unwrap(),
    ]);
    assert_eq!(
        std::fs::read(&plain).unwrap(),
        std::fs::read(&with).unwrap(),
        "indel tables changed the dada JSON",
    );

    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&with).unwrap()).unwrap();
    let trans: Vec<u64> = v["aux"]["transitions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap())
        .collect();
    let ncol = v["aux"]["transitions_ncol"].as_u64().unwrap() as usize;
    let total: u64 = trans.iter().sum();
    let subs: u64 = (0..16)
        .filter(|t| t / 4 != t % 4)
        .map(|t| trans[t * ncol..(t + 1) * ncol].iter().sum::<u64>())
        .sum();

    let tsv = tsv_body(&opp);
    let mut lines = tsv.lines();
    assert_eq!(
        lines.next().unwrap(),
        "sample\ttable\tkey\topportunities\tevents"
    );
    let (mut opp_sum, mut ev_sum) = (0u64, 0u64);
    for row in lines {
        let c: Vec<&str> = row.split('\t').collect();
        if c[1] == "sub_q" {
            opp_sum += c[3].parse::<u64>().unwrap();
            ev_sum += c[4].parse::<u64>().unwrap();
        }
    }
    assert_eq!((opp_sum, ev_sum), (total, subs), "sub_q must equal trans");

    let ev = std::fs::read_to_string(&events).unwrap();
    assert!(ev.starts_with("# dada2_rs_command: dada\n# dada2_rs_version: "));
    assert!(ev.contains("\n# error_matrix: err_out\n"));
    assert!(tsv_body(&events).starts_with("sample\tcluster\tcenter_len\tunique\t"));
}

/// `dada-pooled` indel tables: per-sample outputs unchanged, events attributed
/// to the samples a merged unique occurs in, rate tables pool-level. On one
/// input the pool is that sample, so the tables must equal `dada`'s.
#[test]
fn dada_pooled_indel_tables() {
    let dir = scratch("indel_pooled");
    let err = err_model();
    let err = err.to_str().unwrap();
    let (f1, f2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let (f1, f2) = (f1.to_str().unwrap(), f2.to_str().unwrap());
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();

    run(&[
        "dada-pooled",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("plain"),
    ]);
    run(&[
        "dada-pooled",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("with"),
        "--indel-events",
        &p("events.tsv"),
        "--indel-opportunities",
        &p("opp.tsv"),
    ]);
    for s in ["sam1F", "sam2F"] {
        assert_eq!(
            std::fs::read(dir.join("plain").join(format!("{s}.json"))).unwrap(),
            std::fs::read(dir.join("with").join(format!("{s}.json"))).unwrap(),
            "indel tables changed {s}.json",
        );
    }
    assert!(
        std::fs::read_to_string(p("events.tsv"))
            .unwrap()
            .starts_with("# dada2_rs_command: dada-pooled\n")
    );
    let events = tsv_body(p("events.tsv"));
    let samples: std::collections::BTreeSet<&str> = events
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(samples, ["sam1F", "sam2F"].into_iter().collect());
    let opp = tsv_body(p("opp.tsv"));
    assert!(opp.lines().skip(1).all(|l| l.starts_with("pooled\t")));

    // One input: pooled tables == dada's, up to the sample label.
    run(&[
        "dada-pooled",
        f1,
        "--error-model",
        err,
        "--output-dir",
        &p("one"),
        "--indel-events",
        &p("one_ev.tsv"),
        "--indel-opportunities",
        &p("one_opp.tsv"),
    ]);
    run(&[
        "dada",
        f1,
        "--error-model",
        err,
        "-o",
        &p("d.json"),
        "--indel-events",
        &p("d_ev.tsv"),
        "--indel-opportunities",
        &p("d_opp.tsv"),
    ]);
    assert_eq!(tsv_body(p("one_ev.tsv")), tsv_body(p("d_ev.tsv")),);
    assert_eq!(
        tsv_body(p("one_opp.tsv")).replace("pooled\t", "sam1F\t"),
        tsv_body(p("d_opp.tsv")),
    );
}

/// Per-sample modes collect each sample's indel tables into one file sorted by
/// sample: multi-input `dada` must equal the single-input runs concatenated,
/// and `dada-pseudo` (round 2) must leave its per-sample JSON unchanged.
#[test]
fn per_sample_modes_indel_tables() {
    let dir = scratch("indel_per_sample");
    let err = err_model();
    let err = err.to_str().unwrap();
    let (f1, f2) = (fixture("sam1F.fastq.gz"), fixture("sam2F.fastq.gz"));
    let (f1, f2) = (f1.to_str().unwrap(), f2.to_str().unwrap());
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let read = |name: &str| std::fs::read_to_string(p(name)).unwrap();
    let read_tsv = |name: &str| tsv_body(p(name));

    let flags = |ev: &str, op: &str| {
        [
            "--indel-events".to_string(),
            p(ev),
            "--indel-opportunities".to_string(),
            p(op),
        ]
    };
    let mut expect_ev = String::new();
    let mut expect_op = String::new();
    for (i, f) in [f1, f2].into_iter().enumerate() {
        let [a, b, c, d] = flags(&format!("ev{i}.tsv"), &format!("op{i}.tsv"));
        run(&[
            "dada",
            f,
            "--error-model",
            err,
            "-o",
            &p(&format!("d{i}.json")),
            &a,
            &b,
            &c,
            &d,
        ]);
        let (ev, op) = (
            read_tsv(&format!("ev{i}.tsv")),
            read_tsv(&format!("op{i}.tsv")),
        );
        let skip = |t: &str| {
            if i == 0 {
                t.to_string()
            } else {
                t.split_once('\n').unwrap().1.to_string()
            }
        };
        expect_ev += &skip(&ev);
        expect_op += &skip(&op);
    }
    let [a, b, c, d] = flags("multi_ev.tsv", "multi_op.tsv");
    run(&[
        "dada",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("multi"),
        "--threads",
        "2",
        &a,
        &b,
        &c,
        &d,
    ]);
    assert_eq!(read_tsv("multi_ev.tsv"), expect_ev);
    assert_eq!(read_tsv("multi_op.tsv"), expect_op);

    run(&[
        "dada-pseudo",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("ps_plain"),
    ]);
    let [a, b, c, d] = flags("ps_ev.tsv", "ps_op.tsv");
    run(&[
        "dada-pseudo",
        f1,
        f2,
        "--error-model",
        err,
        "--output-dir",
        &p("ps_with"),
        &a,
        &b,
        &c,
        &d,
    ]);
    for s in ["sam1F", "sam2F"] {
        assert_eq!(
            read(&format!("ps_plain/{s}.json")),
            read(&format!("ps_with/{s}.json")),
            "indel tables changed pseudo {s}.json",
        );
    }
    let header = read("ps_ev.tsv");
    assert!(header.starts_with("# dada2_rs_command: dada-pseudo\n"));
    assert!(header.contains("\n# round: 2\n# reestimated_err_between_rounds: false\n"));
    for t in ["ps_ev.tsv", "ps_op.tsv"] {
        let samples: std::collections::BTreeSet<String> = read_tsv(t)
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect();
        assert_eq!(
            samples,
            ["sam1F", "sam2F"].map(String::from).into_iter().collect(),
            "{t}"
        );
    }

    // R re-estimates the error model between rounds; round 1 feeds the refit,
    // so the tables must not disturb it either.
    let re = "--reestimate-err-between-rounds";
    let base = [
        "dada-pseudo",
        f1,
        f2,
        "--error-model",
        err,
        re,
        "--output-dir",
    ];
    run(&[&base[..], &[&p("re_plain")]].concat());
    run(&[&base[..], &[&p("re_with"), &a, &b, &c, &d]].concat());
    for s in ["sam1F", "sam2F"] {
        assert_eq!(
            read(&format!("re_plain/{s}.json")),
            read(&format!("re_with/{s}.json")),
            "{re} changed {s}.json",
        );
    }
    assert!(read("ps_op.tsv").contains("\n# reestimated_err_between_rounds: true\n"));
}

/// dada-pseudo output must be accepted by the downstream consumers that key off
/// the `dada2_rs_command` tag (merge-pairs and make-sequence-table). This is the
/// path that broke in benchmarking: those readers allowlisted only "dada" /
/// "dada-pooled". The reverse error model is shared with the forward one here —
/// we're asserting the tag is *accepted*, not denoising correctness.
#[test]
fn dada_pseudo_output_feeds_downstream_steps() {
    let dir = scratch("pseudo_downstream");
    let err = err_model();
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

/// With an impossibly high --min-overlap every pair fails to merge. By default
/// those pairs are dropped; with --rescue-unmerged they are concatenated and
/// accepted (marked `concatenated: true`).
#[test]
fn merge_pairs_rescue_unmerged_concatenates() {
    let dir = scratch("merge_rescue");
    let err = err_model();
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

/// A single-input cluster trace carries the output JSON's sample name, so
/// `--sample-name` reaches it and a `.fastq.gz` input does not leave `.fastq`
/// behind (#67).
#[test]
fn dada_cluster_trace_uses_output_sample_name() {
    let dir = scratch("trace_label");
    let err = err_model();
    let input = fixture("sam1F.fastq.gz");
    for (tag, name) in [("stem", None), ("named", Some("S1"))] {
        let out = dir.join(format!("{tag}.json"));
        let trace = dir.join(format!("{tag}.trace.json"));
        let mut args = vec![
            "dada",
            "--error-model",
            err.to_str().unwrap(),
            "--cluster-trace",
            trace.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ];
        if let Some(n) = name {
            args.extend(["--sample-name", n]);
        }
        run(&args);
        let read = |p: &Path| -> serde_json::Value {
            serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap()
        };
        let expected = name.unwrap_or("sam1F");
        assert_eq!(read(&out)["sample"], expected);
        assert_eq!(read(&trace)["sample"], expected);
    }
}

/// `--sample-names` is comma-separated on every subcommand that takes it, as
/// on dada-pooled. make-sequence-table and merge-pairs used to take a
/// space-separated list, so a comma list became one name and failed the
/// length check.
#[test]
fn sample_names_are_comma_separated() {
    let dir = scratch("sample_names_commas");
    let err = err_model();
    let (f1, r1) = (fixture("sam1F.fastq.gz"), fixture("sam1R.fastq.gz"));
    let (f2, r2) = (fixture("sam2F.fastq.gz"), fixture("sam2R.fastq.gz"));
    let dada = |inp: &Path, out: &Path| {
        run(&[
            "dada",
            inp.to_str().unwrap(),
            "--error-model",
            err.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]);
    };
    let j = |n: &str| dir.join(n);
    for (inp, out) in [
        (&f1, "f1.json"),
        (&r1, "r1.json"),
        (&f2, "f2.json"),
        (&r2, "r2.json"),
    ] {
        dada(inp, &j(out));
    }
    let samples = |p: &Path, key: &str| -> Vec<String> {
        let v: serde_json::Value = serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap();
        v["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().or(s[key].as_str()).unwrap().to_string())
            .collect()
    };

    let seqtab = j("seqtab.json");
    run(&[
        "make-sequence-table",
        "--sample-names",
        "A,B",
        "-o",
        seqtab.to_str().unwrap(),
        j("f1.json").to_str().unwrap(),
        j("f2.json").to_str().unwrap(),
    ]);
    assert_eq!(samples(&seqtab, "sample"), ["A", "B"]);

    let merged = j("merged.json");
    run(&[
        "merge-pairs",
        "--fwd-dada",
        j("f1.json").to_str().unwrap(),
        j("f2.json").to_str().unwrap(),
        "--rev-dada",
        j("r1.json").to_str().unwrap(),
        j("r2.json").to_str().unwrap(),
        "--fwd-fastq",
        f1.to_str().unwrap(),
        f2.to_str().unwrap(),
        "--rev-fastq",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "--sample-names",
        "A,B",
        "-o",
        merged.to_str().unwrap(),
    ]);
    assert_eq!(samples(&merged, "sample"), ["A", "B"]);
}

/// `dada*` run from a derep JSON records the FASTQ behind it as `source_fastq`,
/// so `merge-pairs`' provenance warning stays silent on aligned lists and
/// still fires when the FASTQ lists are swapped (#111).
#[test]
fn merge_pairs_provenance_from_derep_json() {
    let dir = scratch("merge_provenance");
    let err = err_model();
    let names = ["sam1F", "sam2F", "sam1R", "sam2R"];
    for name in names {
        let fq = fixture(&format!("{name}.fastq.gz"));
        run(&[
            "derep",
            fq.to_str().unwrap(),
            "-o",
            dir.join(format!("{name}.derep.json")).to_str().unwrap(),
        ]);
    }
    let derep = |name: &str| dir.join(format!("{name}.derep.json"));
    let fastq = |name: &str| fixture(&format!("{name}.fastq.gz"));
    const WARNING: &str = "check that the file lists line up";

    for mode in ["dada", "dada-pooled", "dada-pseudo"] {
        let (fwd, rev) = (dir.join(format!("{mode}_F")), dir.join(format!("{mode}_R")));
        for (pair, out) in [(["sam1F", "sam2F"], &fwd), (["sam1R", "sam2R"], &rev)] {
            run(&[
                mode,
                derep(pair[0]).to_str().unwrap(),
                derep(pair[1]).to_str().unwrap(),
                "--error-model",
                err.to_str().unwrap(),
                "--output-dir",
                out.to_str().unwrap(),
                "--threads",
                "1",
            ]);
        }
        let v: serde_json::Value =
            serde_json::from_slice(&std::fs::read(fwd.join("sam1F.json")).unwrap()).unwrap();
        assert_eq!(v["input_file"], "sam1F.derep.json", "{mode}");
        assert_eq!(v["source_fastq"], "sam1F.fastq.gz", "{mode}");

        let merge = |fwd_fq: [&str; 2]| {
            let out = std::process::Command::new(common::BIN)
                .args([
                    "merge-pairs",
                    "--fwd-dada",
                    fwd.join("sam1F.json").to_str().unwrap(),
                    fwd.join("sam2F.json").to_str().unwrap(),
                    "--rev-dada",
                    rev.join("sam1R.json").to_str().unwrap(),
                    rev.join("sam2R.json").to_str().unwrap(),
                    "--fwd-fastq",
                    fastq(fwd_fq[0]).to_str().unwrap(),
                    fastq(fwd_fq[1]).to_str().unwrap(),
                    "--rev-fastq",
                    fastq("sam1R").to_str().unwrap(),
                    fastq("sam2R").to_str().unwrap(),
                    "-o",
                    dir.join(format!("{mode}_merged.json")).to_str().unwrap(),
                ])
                .output()
                .unwrap();
            String::from_utf8_lossy(&out.stderr).into_owned()
        };
        let aligned = merge(["sam1F", "sam2F"]);
        assert!(!aligned.contains(WARNING), "{mode}: {aligned}");
        let swapped = merge(["sam2F", "sam1F"]);
        assert!(swapped.contains(WARNING), "{mode}: {swapped}");
    }
}
