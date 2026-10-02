//! Hand-offs between pipeline steps: derep ordering, FASTQ vs derep-JSON input,
//! the failed-uniques TSV, merge-pairs, and downstream acceptance of each
//! command's output tag.

mod common;

use std::path::Path;

use common::{fixture, run, scratch, shared_err_model};

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
