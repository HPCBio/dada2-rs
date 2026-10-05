//! Denoising: `dada`, `dada-pooled`, `dada-pseudo`, and the helpers they share.

use std::{
    fs::File,
    io,
    path::{Path, PathBuf},
    sync::Mutex,
};

use flate2::read::MultiGzDecoder;

use dada2_rs::{
    cli, cluster_trace, containers, dada, derep, error_models, failed_uniques, learn_errors,
    metrics, minimizers, misc, nwalign, sequence_table,
};

use containers::BirthType;
use derep::dereplicate;
use learn_errors::{ErrFun, LearnedErrParams};
use metrics::{MeasureLevel, MetricsDocument};
use misc::{Tagged, read_fasta_records, read_tagged_json};
use nwalign::{AlignBackend, AlignParams, ScreenBackend};
use sequence_table::{HashAlgo, SequenceTable};
use serde::Serialize;

use super::common::{
    WFA_MAX_EDITS_DEFAULT, check_input_paths, fastq_stem, file_basename, note_homopolymer_gapping,
    select_sequences,
};
use error_models::LoessConfig;
use misc::WithPath;

pub(crate) fn run_dada(args: cli::DadaArgs) -> io::Result<()> {
    let cli::DadaArgs {
        input,
        error_model,
        use_err_in,
        sample_name,
        prior,
        inherit_err_params,
        phred_offset,
        threads,
        sample_jobs,
        denoise,
        experimental,
        aux_outputs,
        cluster_trace,
        trace_no_members,
        trace_min_abund,
        output,
        output_dir,
        failed_uniques: failed_uniques_path,
        compact,
        gzip,
        metrics_json,
        metrics_attribution,
        verbose,
    } = args;
    // Wall clock for `--metrics-json`, started before any I/O so the
    // document measures the subcommand and not just the denoiser.
    let t_start = std::time::Instant::now();
    let measure_level = resolve_measure_level(verbose, metrics_json.as_ref(), metrics_attribution);
    check_input_paths("input", &input)?;
    // ---- Load uniques from FASTQ or a derep/sample JSON ----
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    // ---- Multi-input: independent per-sample (serial, NOT pooled) ----
    if input.len() > 1 {
        // Reject single-sample-only flags.
        let mut bad: Vec<&str> = Vec::new();
        if output.is_some() {
            bad.push("--output/-o");
        }
        if sample_name.is_some() {
            bad.push("--sample-name");
        }
        if aux_outputs {
            bad.push("--aux-outputs");
        }
        if cluster_trace.is_some() {
            bad.push("--cluster-trace");
        }
        if trace_no_members {
            bad.push("--trace-no-members");
        }
        // `trace_min_abund` has a default of 1; only flag if non-default.
        if trace_min_abund != 1 {
            bad.push("--trace-min-abund");
        }
        if !bad.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{} {} single-sample-only and cannot be combined with multiple inputs; use --output-dir for per-sample output",
                    bad.join(", "),
                    if bad.len() == 1 { "is" } else { "are" },
                ),
            ));
        }
        let output_dir = output_dir.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "more than one input given: --output-dir is required (one {sample}.json per sample is written there)",
            )
        })?;
        std::fs::create_dir_all(&output_dir)?;

        // Load the optional prior set once.
        let prior_set = prior.as_deref().map(read_prior_set).transpose()?;

        // Resolve parameters once (shared across samples).
        let resolved = resolve_dada_params(
            &error_model,
            use_err_in,
            inherit_err_params,
            threads,
            false, // aux_outputs (rejected above for multi-input)
            false, // pool
            verbose,
            measure_level,
            denoise,
            experimental,
        )?;

        let jobs = concurrent_samples(sample_jobs, threads, input.len());
        if verbose {
            eprintln!(
                "[dada] denoising {} sample(s), {jobs} concurrent",
                input.len()
            );
        }
        let collect_failed = failed_uniques_path.is_some();
        let sink = SampleSink::new("dada", &output_dir, gzip, verbose, None);
        for_each_sample_concurrent(input.len(), jobs, threads, |i, sub_pool| {
            let path = &input[i];
            let (mut raw_inputs, json_sample) =
                load_sample_raws(path, phred_offset, sub_pool, verbose)?;
            if let Some(ref set) = prior_set {
                let n_marked = mark_priors(&mut raw_inputs, set);
                if verbose {
                    eprintln!(
                        "[dada] {} of {} unique(s) marked as prior in {}",
                        n_marked,
                        raw_inputs.len(),
                        path.display(),
                    );
                }
            }
            let sample = sample_label(None, json_sample, path);
            let out = denoise_and_serialize(
                "dada",
                &sample,
                &file_basename(path),
                &raw_inputs,
                &resolved.params,
                &resolved.run,
                sub_pool,
                compact,
                collect_failed,
                verbose,
                (jobs > 1).then(|| sample.clone()),
                None,
            )?;
            sink.write(&sample, out)
        })?;
        sink.finish(
            failed_uniques_path.as_deref(),
            metrics_json.as_deref(),
            t_start,
            measure_level,
        )?;
        return Ok(());
    }

    // ---- Single input: -o/stdout, plus --sample-name, aux outputs and trace ----
    if output_dir.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--output-dir is only valid with multiple inputs; for a single input use --output/-o or stdout",
        ));
    }
    let input = &input[0];

    let (mut raw_inputs, json_sample) = load_sample_raws(input, phred_offset, &pool, verbose)?;

    // ---- Mark prior sequences ----
    if let Some(ref prior_path) = prior {
        let prior_seqs = read_prior_set(prior_path)?;
        let n_marked = mark_priors(&mut raw_inputs, &prior_seqs);
        if verbose {
            eprintln!(
                "[dada] {} of {} unique(s) marked as prior from {}",
                n_marked,
                raw_inputs.len(),
                prior_path.display(),
            );
        }
    }

    // ---- Load error model + resolve parameters (shared helper) ----
    let resolved = resolve_dada_params(
        &error_model,
        use_err_in,
        inherit_err_params,
        threads,
        aux_outputs,
        false, // pool
        verbose,
        measure_level,
        denoise,
        experimental,
    )?;
    let sample = sample_label(sample_name, json_sample, input);
    let trace = cluster_trace.as_deref().map(|path| TraceRequest {
        path,
        params: cluster_trace::TraceParams {
            no_members: trace_no_members,
            min_abund: trace_min_abund,
        },
    });
    let (json, failed, run_metrics) = denoise_and_serialize(
        "dada",
        &sample,
        &file_basename(input),
        &raw_inputs,
        &resolved.params,
        &resolved.run,
        &pool,
        compact,
        failed_uniques_path.is_some(),
        verbose,
        None,
        trace.as_ref(),
    )?;

    if let Some(ref fu_path) = failed_uniques_path {
        write_failed_uniques("dada", fu_path, failed, verbose)?;
    }
    if let (Some(mpath), Some(m)) = (metrics_json.as_ref(), run_metrics) {
        write_run_metrics(
            "dada",
            mpath,
            vec![(sample, None, m)],
            t_start,
            measure_level,
            verbose,
        )?;
    }

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

pub(crate) fn run_dada_pooled(args: cli::DadaPooledArgs) -> io::Result<()> {
    let cli::DadaPooledArgs {
        input,
        error_model,
        use_err_in,
        prior,
        inherit_err_params,
        sample_names,
        output_dir,
        phred_offset,
        threads,
        denoise,
        pool_tiebreak,
        experimental,
        failed_uniques: failed_uniques_path,
        pooled_record,
        cluster_trace,
        trace_no_members,
        trace_min_abund,
        compact,
        gzip,
        metrics_json,
        metrics_attribution,
        verbose,
    } = args;
    let t_start = std::time::Instant::now();
    let measure_level = resolve_measure_level(verbose, metrics_json.as_ref(), metrics_attribution);
    check_input_paths("input", &input)?;

    let n_samples = input.len();
    check_sample_names(sample_names.as_deref(), n_samples)?;

    std::fs::create_dir_all(&output_dir)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    // ---- Streaming per-sample dereplication + merge (#41) ----
    // Load one sample, fold it into the pool, drop it, so per-sample dereps
    // are never all resident. Both steps are serial, unlike `run_dada`; they
    // are timed apart for the phase report.
    let t_derep = std::time::Instant::now();
    let mut t_merge_acc = std::time::Duration::ZERO;
    let mut derep_cost = DerepLoadCost::default();

    let mut json_samples: Vec<Option<String>> = vec![None; n_samples];
    let mut merged = derep::DerepPool::new();

    // (reads, uniques) per input, for `pipeline.inputs`. The unnamed
    // `[derep]` line is suppressed here: this is its named home.
    let mut input_counts: Vec<(u64, usize)> = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let (derep, name) =
            load_derep_for_dada(&input[i], phred_offset, &pool, false, &mut derep_cost)?;
        json_samples[i] = name;
        input_counts.push((
            derep.uniques.iter().map(|(_, c)| c).sum(),
            derep.uniques.len(),
        ));

        let t_m = std::time::Instant::now();
        merged.add(&derep);
        t_merge_acc += t_m.elapsed();
        // `derep` (this sample's full quals + seqs) drops here.
    }

    // Split the interleaved wall into load vs fold for the phase report:
    // total loop time minus folding ≈ the serial derep/load front.
    let t_merge = t_merge_acc;
    let t_derep = t_derep.elapsed().saturating_sub(t_merge);
    // One getrusage each, so read regardless of verbosity.
    let rss_after_derep_merge = misc::peak_rss_kb();
    if verbose {
        if let Some(split) = derep_cost.report(t_derep) {
            eprintln!("{split}");
        }
        eprintln!(
            "[dada-pooled] peak RSS after derep+merge: {} MB",
            rss_after_derep_merge / 1024
        );
    }
    let derep::PooledDerep {
        seqs: merged_seqs,
        qual_sums: merged_qual_sum,
        abundance: merged_total,
        local_to_merged,
        sample_counts: sample_unique_counts,
    } = merged.finish(pool_tiebreak);

    // Resolve sample names: CLI override > JSON-embedded > filename stem.
    let sample_names: Vec<String> = match sample_names {
        Some(names) => names,
        None => input
            .iter()
            .zip(json_samples)
            .map(|(p, js)| sample_label(None, js, p))
            .collect(),
    };

    let n_merged = merged_seqs.len();
    if verbose {
        eprintln!(
            "[dada-pooled] {} sample(s) → {} merged unique(s), {} total reads; ties in {} order",
            n_samples,
            n_merged,
            merged_total.iter().sum::<u32>(),
            pool_tiebreak.label(),
        );
    }

    // ---- Build merged RawInput list ----
    // Consume the merge vecs (move, don't clone): each merged sequence
    // and its u32 qual-sum row are moved straight into a RawInput, so the
    // merge intermediates are freed as raw_inputs is built rather than
    // sitting resident through the dada call (issue #39). quals are
    // already u32 sums — no conversion.
    let mut raw_inputs: Vec<dada::RawInput> = merged_seqs
        .into_iter()
        .zip(merged_qual_sum)
        .zip(merged_total)
        .map(|((seq_bytes, quals), abundance)| {
            let sequence: String = String::from_utf8(seq_bytes)
                .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
            dada::RawInput {
                quals: Some(quals),
                seq: sequence,
                abundance,
                prior: false,
            }
        })
        .collect();

    if raw_inputs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "All input FASTQ files contain no reads",
        ));
    }

    // ---- Mark prior sequences ----
    if let Some(ref prior_path) = prior {
        let prior_seqs = read_prior_set(prior_path)?;
        let n_marked = mark_priors(&mut raw_inputs, &prior_seqs);
        if verbose {
            eprintln!(
                "[dada-pooled] {} of {} merged unique(s) marked as prior from {}",
                n_marked,
                raw_inputs.len(),
                prior_path.display(),
            );
        }
    }

    // ---- Load error model + resolve parameters (shared helper) ----
    let resolved = resolve_dada_params(
        &error_model,
        use_err_in,
        inherit_err_params,
        threads,
        false, // aux_outputs
        true,  // pool
        verbose,
        measure_level,
        denoise,
        experimental,
    )?;
    let dada_params = resolved.params;
    let mut run_params = resolved.run;
    run_params.n_prior = raw_inputs.iter().filter(|r| r.prior).count();
    run_params.pool_tiebreak = Some(pool_tiebreak);

    // ---- Run DADA once on the merged table ----
    let rss_after_merge = misc::peak_rss_kb();
    if verbose {
        eprintln!(
            "[dada-pooled] peak RSS after merge: {} MB",
            rss_after_merge / 1024
        );
    }
    let t_dada = std::time::Instant::now();
    let mut result = pool
        .install(|| dada::dada_uniques(&raw_inputs, &dada_params))
        .map_err(io::Error::other)?;
    let t_dada = t_dada.elapsed();
    let pooled_metrics = result.metrics.take();
    let rss_after_dada = misc::peak_rss_kb();
    if verbose {
        eprintln!(
            "[dada-pooled] peak RSS after dada: {} MB",
            rss_after_dada / 1024
        );
    }

    if verbose {
        eprintln!(
            "[dada-pooled] {} ASV(s) from {} merged unique(s); {} aligns, {} shrouded",
            result.clusters.len(),
            raw_inputs.len(),
            result.nalign,
            result.nshroud,
        );
    }

    // ---- Cluster trace (diagnostics; write-only, does not affect ASVs) ----
    if let Some(ref path) = cluster_trace {
        let trace = TraceRequest {
            path,
            params: cluster_trace::TraceParams {
                no_members: trace_no_members,
                min_abund: trace_min_abund,
            },
        };
        write_cluster_trace(
            "dada-pooled",
            &trace,
            "pooled",
            &raw_inputs,
            &result,
            &dada_params,
            compact,
            verbose,
        )?;
    }

    // ---- Per-sample output ----
    let t_output = std::time::Instant::now();

    // Failed-to-denoise uniques (issue #60). Pooled denoising runs once on
    // the merged unique table, so "failed" is a global property
    // (`result.map[mu] == None`); rows are collected per sample the failed
    // merged unique appears in, carrying that sample's read count.
    let collect_failed = failed_uniques_path.is_some();
    let mut failed_rows: Vec<failed_uniques::Row> = Vec::new();

    for (s, sample_name) in sample_names.iter().enumerate() {
        // Sum per-cluster reads for this sample by walking its local uniques.
        let mut cluster_reads: Vec<u32> = vec![0u32; result.clusters.len()];
        for (lu, &mu) in local_to_merged[s].iter().enumerate() {
            if let Some(c) = result.map[mu] {
                cluster_reads[c] += sample_unique_counts[s][lu];
            } else if collect_failed {
                failed_rows.push(failed_uniques::Row {
                    sequence: raw_inputs[mu].seq.clone(),
                    sample: sample_name.clone(),
                    reads: sample_unique_counts[s][lu],
                });
            }
        }

        // Filter to clusters present in this sample; renumber globally → locally.
        let mut global_to_local: Vec<Option<usize>> = vec![None; result.clusters.len()];
        let mut asvs: Vec<AsvEntry> = Vec::new();
        for (c, cluster) in result.clusters.iter().enumerate() {
            if cluster_reads[c] == 0 {
                continue;
            }
            global_to_local[c] = Some(asvs.len());
            asvs.push(asv_entry_from_cluster(cluster, cluster_reads[c]));
        }

        let total_reads: u32 = cluster_reads.iter().sum();

        // Per-sample local-unique → local-cluster map (mirrors single-sample dada).
        let map: Vec<Option<usize>> = (0..sample_unique_counts[s].len())
            .map(|lu| {
                let mu = local_to_merged[s][lu];
                result.map[mu].and_then(|c| global_to_local[c])
            })
            .collect();

        let n_asvs = asvs.len();
        let out = DadaOutput {
            sample: sample_name.clone(),
            input_file: file_basename(&input[s]),
            num_asvs: n_asvs,
            total_reads,
            asvs,
            stats: DadaStats {
                nalign: result.nalign,
                nshroud: result.nshroud,
            },
            params: run_params,
            map,
            aux: None,
            pool_input_index: Some(s),
        };

        let json = to_json(&Tagged::new("dada-pooled", out), compact)?;

        let path = output_dir.join(if gzip {
            format!("{sample_name}.json.gz")
        } else {
            format!("{sample_name}.json")
        });
        misc::write_maybe_gz(&path, json.as_bytes())?;
        if verbose {
            eprintln!(
                "[dada-pooled] wrote {} ({} ASV(s), {} reads)",
                path.display(),
                n_asvs,
                total_reads
            );
        }
    }

    // ---- Pooled record (kdist-calibrate --from-dada-pooled) ----
    // The per-sample JSONs above are projections of ONE global partition:
    // inference ran once on the merged unique table, so each unique's
    // fate (`result.map`) and every ASV are pool-level facts. A pool-level
    // screen assessment must not re-aggregate the per-sample splits (that
    // double-counts shared sequences and re-derives the failed singleton
    // split from local counts). When asked (`--pooled-record`), emit one
    // self-contained record carrying the merged uniques with their POOLED
    // abundance, the global `map`, and the global ASVs, so kdist can screen
    // the pool without touching --derep-dir. Off by default and written to
    // an explicit path (NOT --output-dir) so it never lands in the
    // per-sample `*.json.gz` glob downstream steps rely on. gzip follows the
    // path's own `.gz` extension (write_maybe_gz), not `--gzip`.
    if let Some(ref rec_path) = pooled_record {
        #[derive(Serialize)]
        struct PooledUnique {
            sequence: String,
            count: u32,
        }
        #[derive(Serialize)]
        struct PooledInput {
            sample: String,
            input_file: String,
        }
        #[derive(Serialize)]
        struct PooledRecord {
            num_uniques: usize,
            num_asvs: usize,
            uniques: Vec<PooledUnique>,
            map: Vec<Option<usize>>,
            asvs: Vec<AsvEntry>,
            /// Inputs in the order pooled, which decides tied uniques (#260).
            inputs: Vec<PooledInput>,
            pool_tiebreak: derep::PoolTiebreak,
        }
        let pooled_uniques: Vec<PooledUnique> = raw_inputs
            .iter()
            .map(|r| PooledUnique {
                sequence: r.seq.clone(),
                count: r.abundance,
            })
            .collect();
        let pooled_asvs: Vec<AsvEntry> = result
            .clusters
            .iter()
            .map(|c| asv_entry_from_cluster(c, c.abundance))
            .collect();
        let record = PooledRecord {
            num_uniques: pooled_uniques.len(),
            num_asvs: pooled_asvs.len(),
            uniques: pooled_uniques,
            map: result.map.clone(),
            asvs: pooled_asvs,
            inputs: sample_names
                .iter()
                .zip(&input)
                .map(|(sample, path)| PooledInput {
                    sample: sample.clone(),
                    input_file: file_basename(path),
                })
                .collect(),
            pool_tiebreak,
        };
        let json = to_json(&Tagged::new("dada-pooled-record", record), compact)?;
        misc::write_maybe_gz(rec_path, json.as_bytes())?;
        if verbose {
            eprintln!(
                "[dada-pooled] wrote pooled record {} ({} merged unique(s), {} ASV(s))",
                rec_path.display(),
                n_merged,
                result.clusters.len(),
            );
        }
    }

    if let Some(ref fu_path) = failed_uniques_path {
        write_failed_uniques("dada-pooled", fu_path, failed_rows, verbose)?;
    }
    let t_output = t_output.elapsed();
    if let (Some(mpath), Some(m)) = (metrics_json.as_ref(), pooled_metrics) {
        let mut doc = MetricsDocument::new(t_start.elapsed(), measure_level);
        doc.pipeline.derep = Some(t_derep.as_secs_f64());
        doc.pipeline.merge = Some(t_merge.as_secs_f64());
        doc.pipeline.dada = Some(t_dada.as_secs_f64());
        doc.pipeline.output = Some(t_output.as_secs_f64());
        doc.pipeline.derep_detail = derep_cost.detail();
        doc.pipeline.inputs = Some(
            sample_names
                .iter()
                .zip(&input_counts)
                .zip(&derep_cost.per_sample)
                .map(|((name, &(reads, uniques)), t)| metrics::DerepInput {
                    sample: name.clone(),
                    reads,
                    uniques,
                    load_seconds: t.as_secs_f64(),
                })
                .collect(),
        );
        // `peak_rss_kb` returns 0 when getrusage fails; absent, not zero.
        let rss = [rss_after_derep_merge, rss_after_merge, rss_after_dada];
        doc.pipeline.peak_rss_mb = rss.iter().all(|&kb| kb > 0).then(|| metrics::PeakRss {
            after_derep_merge: rss[0] as f64 / 1024.0,
            after_merge: rss[1] as f64 / 1024.0,
            after_dada: rss[2] as f64 / 1024.0,
        });
        // One invocation: pooled denoises the merged table once, so the
        // per-sample outputs are slices of a single run, not runs.
        doc.push("__pooled__", None, m);
        write_metrics_json(mpath, &doc)?;
        if verbose {
            eprintln!("[dada-pooled] wrote run metrics to {}", mpath.display());
        }
    }
    if verbose {
        let total = (t_derep + t_merge + t_dada + t_output)
            .as_secs_f64()
            .max(1e-9);
        let pct = |d: std::time::Duration| 100.0 * d.as_secs_f64() / total;
        eprintln!(
            "[dada-pooled] phase wall times: derep={:.1}s ({:.0}%, serial)  merge={:.1}s ({:.0}%, serial)  run_dada={:.1}s ({:.0}%, parallel)  output={:.1}s ({:.0}%, serial)",
            t_derep.as_secs_f64(),
            pct(t_derep),
            t_merge.as_secs_f64(),
            pct(t_merge),
            t_dada.as_secs_f64(),
            pct(t_dada),
            t_output.as_secs_f64(),
            pct(t_output),
        );
    }
    Ok(())
}

pub(crate) fn run_dada_pseudo(args: cli::DadaPseudoArgs) -> io::Result<()> {
    let cli::DadaPseudoArgs {
        input,
        error_model,
        use_err_in,
        inherit_err_params,
        sample_names,
        output_dir,
        pseudo_prevalence,
        pseudo_min_abundance,
        priors_out,
        reestimate_err_between_rounds,
        phred_offset,
        threads,
        sample_jobs,
        cache_samples,
        denoise,
        experimental,
        failed_uniques: failed_uniques_path,
        compact,
        gzip,
        metrics_json,
        metrics_attribution,
        verbose,
    } = args;
    let t_start = std::time::Instant::now();
    let measure_level = resolve_measure_level(verbose, metrics_json.as_ref(), metrics_attribution);
    check_input_paths("input", &input)?;
    use std::collections::{HashMap, HashSet};

    // Streaming is the default (faster AND lighter on large runs — the
    // retained all-samples cache is pure overhead); --cache-samples opts
    // back into holding every sample's uniques across both rounds.
    let low_memory = !cache_samples;

    let n_samples = input.len();
    let jobs = concurrent_samples(sample_jobs, threads, n_samples);
    check_sample_names(sample_names.as_deref(), n_samples)?;

    std::fs::create_dir_all(&output_dir)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    // Resolve parameters once (shared across both rounds and all samples).
    let mut resolved = resolve_dada_params(
        &error_model,
        use_err_in,
        inherit_err_params,
        threads,
        // aux_outputs: round 1 only needs the extra final-subs pass when
        // we are going to re-fit from its transition counts. Flipped back
        // off before round 2 below, so the flag costs nothing when unset.
        reestimate_err_between_rounds,
        false, // pool (pseudo is per-sample, not pooled)
        verbose,
        measure_level,
        denoise,
        experimental,
    )?;

    // ---- Validate the re-estimation request up front ----
    // Both conditions are cheap to check and expensive to hit late: the
    // transition matrix is only collected per quality column when quals
    // are in use, and re-fitting needs the errfun the model was built
    // with. Failing here beats failing after round 1.
    if reestimate_err_between_rounds {
        if !resolved.params.use_quals {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--reestimate-err-between-rounds requires quality scores \
                 (transition counts collapse to a single column with \
                 --use-quals false)",
            ));
        }
        if resolved.err_params.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "--reestimate-err-between-rounds needs the error model's \
                 `params` block to know which errfun to re-fit with; this \
                 model has none (re-run learn-errors to produce one)",
            ));
        }
    }
    // Accumulated 16 × nq transition counts from round 1, R's
    // `accumulateTrans`. Counts are integers, so folding them in whatever
    // order samples happen to finish is exact — the result does not
    // depend on --sample-jobs.
    let trans_acc: Mutex<Vec<u32>> = Mutex::new(vec![0u32; 16 * resolved.nq]);

    // ---- Round 1: denoise each sample independently (no priors) ----
    // Produces sample_names + per-sample round-1 ASVs. The default
    // (streaming) drops each sample's uniques after denoising and
    // re-reads them in round 2, capping peak memory at `jobs` samples in
    // flight instead of all of them; --cache-samples opts into keeping
    // every sample's uniques resident across both rounds.
    if verbose {
        eprintln!(
            "[dada-pseudo] round 1: {n_samples} sample(s), no priors ({jobs} concurrent{})",
            if low_memory { ", streaming" } else { "" },
        );
    }
    type R1 = (
        Vec<String>,
        Vec<Vec<(String, u32)>>,
        Option<Vec<Vec<dada::RawInput>>>,
    );
    let (sample_names, round1_asvs, mut sample_raws_opt): R1 = if !low_memory {
        // Cached: pre-load all uniques, then denoise from the cache.
        let mut sample_raws: Vec<Vec<dada::RawInput>> = Vec::with_capacity(n_samples);
        let mut json_samples: Vec<Option<String>> = Vec::with_capacity(n_samples);
        for path in &input {
            let (raws, js) = load_sample_raws(path, phred_offset, &pool, verbose)?;
            sample_raws.push(raws);
            json_samples.push(js);
        }
        let names: Vec<String> = match &sample_names {
            Some(n) => n.clone(),
            None => input
                .iter()
                .zip(&json_samples)
                .map(|(p, js)| sample_label(None, js.clone(), p))
                .collect(),
        };
        let collected: Mutex<IndexedAsvs> = Mutex::new(Vec::with_capacity(n_samples));
        for_each_sample_concurrent(n_samples, jobs, threads, |s, sub_pool| {
            let result = sub_pool
                .install(|| dada::dada_uniques(&sample_raws[s], &resolved.params))
                .map_err(io::Error::other)?;
            accumulate_round1_trans(&trans_acc, &result)?;
            let asvs = result_to_asvs(&result);
            if verbose {
                eprintln!(
                    "[dada-pseudo]   round 1 {}: {} ASV(s)",
                    names[s],
                    asvs.len()
                );
            }
            collected.lock().unwrap().push((s, asvs));
            Ok(())
        })?;
        let mut r1: Vec<Vec<(String, u32)>> = vec![Vec::new(); n_samples];
        for (s, asvs) in collected.into_inner().unwrap() {
            r1[s] = asvs;
        }
        (names, r1, Some(sample_raws))
    } else {
        // Streaming: load + denoise + drop per sample; capture name + ASVs.
        let collected: Mutex<IndexedNamedAsvs> = Mutex::new(Vec::with_capacity(n_samples));
        for_each_sample_concurrent(n_samples, jobs, threads, |s, sub_pool| {
            let (raws, js) = load_sample_raws(&input[s], phred_offset, sub_pool, verbose)?;
            let result = sub_pool
                .install(|| dada::dada_uniques(&raws, &resolved.params))
                .map_err(io::Error::other)?;
            accumulate_round1_trans(&trans_acc, &result)?;
            let asvs = result_to_asvs(&result);
            let name = sample_label(sample_names.as_ref().map(|n| n[s].clone()), js, &input[s]);
            if verbose {
                eprintln!("[dada-pseudo]   round 1 {}: {} ASV(s)", name, asvs.len());
            }
            collected.lock().unwrap().push((s, name, asvs));
            Ok(())
        })?;
        let mut names = vec![String::new(); n_samples];
        let mut r1: Vec<Vec<(String, u32)>> = vec![Vec::new(); n_samples];
        for (s, name, asvs) in collected.into_inner().unwrap() {
            names[s] = name;
            r1[s] = asvs;
        }
        (names, r1, None)
    };

    // ---- Build a sequence table from round-1 ASVs (samples × sequences) ----
    let mut seq_index: HashMap<String, usize> = HashMap::new();
    let mut sequences: Vec<String> = Vec::new();
    for asvs in &round1_asvs {
        for (seq, _) in asvs {
            if !seq_index.contains_key(seq) {
                seq_index.insert(seq.clone(), sequences.len());
                sequences.push(seq.clone());
            }
        }
    }
    let nseq = sequences.len();
    let mut counts: Vec<Vec<u64>> = vec![vec![0u64; nseq]; n_samples];
    for (s, asvs) in round1_asvs.iter().enumerate() {
        for (seq, abund) in asvs {
            let j = seq_index[seq];
            counts[s][j] += *abund as u64;
        }
    }
    let table = SequenceTable {
        samples: sample_names.clone(),
        sequence_ids: sequences.iter().map(|s| HashAlgo::Sha1.digest(s)).collect(),
        sequences,
        counts,
    };

    // ---- Select round-2 priors via R's PSEUDO rule ----
    let selected = select_sequences(&table, Some(pseudo_prevalence), pseudo_min_abundance);
    let prior_set: HashSet<String> = selected
        .iter()
        .map(|&j| table.sequences[j].to_ascii_uppercase())
        .collect();

    if verbose {
        eprintln!(
            "[dada-pseudo] selected {} prior sequence(s) from {} round-1 ASV(s) (prevalence>={}{})",
            prior_set.len(),
            table.sequences.len(),
            pseudo_prevalence,
            match pseudo_min_abundance {
                Some(m) => format!(" or total-abundance>={m}"),
                None => String::new(),
            },
        );
    }

    // ---- Optional priors FASTA dump ----
    if let Some(ref priors_path) = priors_out {
        if let Some(parent) = priors_path.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent)?;
        }
        let mut fasta = String::new();
        for (i, &j) in selected.iter().enumerate() {
            fasta.push_str(&format!(">prior{}\n{}\n", i + 1, table.sequences[j]));
        }
        std::fs::write(priors_path, fasta)?;
        if verbose {
            eprintln!("[dada-pseudo] wrote priors to {}", priors_path.display());
        }
    }

    if verbose && (n_samples < 2 || prior_set.is_empty()) {
        eprintln!(
            "[dada-pseudo] note: {} — round 2 is equivalent to round 1 (no priors applied)",
            if n_samples < 2 {
                "fewer than 2 samples".to_string()
            } else {
                "no sequences met the prior-selection threshold".to_string()
            },
        );
    }

    // ---- Optional: re-estimate the error model from round 1 (#100) ----
    // R DADA2's pool="pseudo" runs both rounds inside the self-consistency
    // loop, which re-fits `err` from the just-finished round's transitions
    // (dada.R:371-378) with no selfConsist guard — so its round 2 uses a
    // model derived from round 1, not the supplied one. Off by default:
    // the published definition of pseudo-pooling is priors-only, and
    // whether R's re-fit is intended is still open (issue #100).
    if reestimate_err_between_rounds {
        let trans = std::mem::take(&mut *trans_acc.lock().unwrap());
        let err_params = resolved
            .err_params
            .as_ref()
            .expect("validated above: err_params present");
        let errfun = errfun_from_learned(err_params)?;
        let new_err = errfun
            .apply(&trans, resolved.nq)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if new_err.len() != resolved.params.err_mat.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "re-fitted error model is {} entries, expected {}",
                    new_err.len(),
                    resolved.params.err_mat.len()
                ),
            ));
        }
        if verbose {
            let max_delta = resolved
                .params
                .err_mat
                .iter()
                .zip(&new_err)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            eprintln!(
                "[dada-pseudo] re-estimated error model from round 1 \
                 (errfun {}, max |delta| = {max_delta:.2e}); round 2 will \
                 NOT use the supplied model",
                err_params.errfun,
            );
        }
        resolved.params.err_mat = new_err;
        resolved.run.reestimated_err_between_rounds = true;
    }
    // Round 1's extra alignment pass is done with; do not pay for it again.
    resolved.params.aux_outputs = false;

    // ---- Round 2: re-denoise each sample with the priors flagged ----
    if verbose {
        eprintln!(
            "[dada-pseudo] round 2: denoising with priors ({jobs} concurrent{})",
            if low_memory { ", streaming" } else { "" },
        );
    }
    let collect_failed = failed_uniques_path.is_some();
    // Round 2 only. Round 1 runs through a different path that does not
    // serialize per-sample output, so its metrics are not collected
    // here; the document says `round: 2` rather than implying it covers
    // both turns.
    let sink = SampleSink::new("dada-pseudo", &output_dir, gzip, verbose, Some(2));
    // Cached: mark priors serially up front (cheap), then denoise concurrently
    // with shared read access.
    if let Some(sample_raws) = sample_raws_opt.as_mut() {
        for (raws, sample_name) in sample_raws.iter_mut().zip(&sample_names) {
            mark_round2_priors(raws, &prior_set, sample_name, verbose);
        }
    }
    let cached = sample_raws_opt.as_ref();
    for_each_sample_concurrent(n_samples, jobs, threads, |s, sub_pool| {
        let sample_name = &sample_names[s];
        // Streaming: re-load the sample and mark its priors here, so peak
        // memory stays bounded by `jobs` samples.
        let loaded;
        let raws = match cached {
            Some(sample_raws) => &sample_raws[s],
            None => {
                let (mut raws, _js) = load_sample_raws(&input[s], phred_offset, sub_pool, verbose)?;
                mark_round2_priors(&mut raws, &prior_set, sample_name, verbose);
                loaded = raws;
                &loaded
            }
        };
        let out = denoise_and_serialize(
            "dada-pseudo",
            sample_name,
            &file_basename(&input[s]),
            raws,
            &resolved.params,
            &resolved.run,
            sub_pool,
            compact,
            collect_failed,
            verbose,
            (jobs > 1).then(|| sample_name.to_string()),
            None,
        )?;
        sink.write(sample_name, out)
    })?;
    sink.finish(
        failed_uniques_path.as_deref(),
        metrics_json.as_deref(),
        t_start,
        measure_level,
    )?;
    Ok(())
}

/// `screen_backend` is omitted from output when it is the default, so a k-mer run's
/// JSON is unchanged from before the backend existed.
fn is_default_screen(b: &ScreenBackend) -> bool {
    *b == ScreenBackend::Kmer
}

#[derive(Serialize, Copy, Clone)]
struct DadaRunParams {
    omega_a: f64,
    omega_c: f64,
    omega_p: f64,
    min_fold: f64,
    min_hamming: u32,
    min_abund: u32,
    detect_singletons: bool,
    band: i32,
    homo_gap_p: i32,
    gap_p: i32,
    match_score: i32,
    mismatch: i32,
    max_clust: usize,
    greedy: bool,
    use_quals: bool,
    kdist_cutoff: f64,
    kmer_size: usize,
    use_kmers: bool,
    /// Which pre-alignment screen ran. Recorded so a run's output says what
    /// produced it: the two backends give closely-agreeing but not identical
    /// tables, so "which screen" is part of a result's provenance, not a detail.
    ///
    /// Omitted for the default `kmer` backend, which keeps every existing output
    /// byte-identical (AGENTS.md: do not alter the flat JSON output shape). Its
    /// absence therefore means `kmer`.
    #[serde(skip_serializing_if = "is_default_screen")]
    screen_backend: ScreenBackend,
    /// Sketch parameters, emitted only under the minimizer backend so every
    /// existing k-mer output keeps its exact shape (AGENTS.md: do not alter the
    /// flat JSON output shape).
    #[serde(skip_serializing_if = "Option::is_none")]
    minimizer_k: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    minimizer_w: Option<usize>,
    /// Denoising mode: false = independent per-sample (`dada`), true = full
    /// pooling across samples (`dada-pooled`, R DADA2 pool=TRUE).
    pool: bool,
    /// Number of unique input sequences flagged as priors (0 when no `--prior`).
    n_prior: usize,
    /// True when `dada-pseudo --reestimate-err-between-rounds` re-fitted the
    /// error model from round 1, so round 2 did NOT use the supplied model.
    /// Skipped when false, keeping every other mode's output byte-identical to
    /// before the flag existed — but a run whose error model changed mid-flight
    /// must say so, since the `--error-model` path alone no longer describes it.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    reestimated_err_between_rounds: bool,
    /// How `dada-pooled` ordered tied uniques (#260); pooled runs only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pool_tiebreak: Option<derep::PoolTiebreak>,
}

/// `true` when `path` looks like a JSON file (`.json` or `.json.gz`).
fn is_json_path(path: &Path) -> bool {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    ext == "json"
        || path
            .file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
            == Some("json")
}

/// Error-model JSON shape shared by `dada`, `dada-pooled`, and `dada-pseudo`.
#[derive(serde::Deserialize)]
struct ErrorModelJson {
    nq: usize,
    err_in: Vec<Vec<f64>>,
    err_out: Vec<Vec<f64>>,
    #[serde(default)]
    params: Option<LearnedErrParams>,
}

/// Bundle of resolved DADA parameters plus the scalars needed to serialize a
/// `DadaRunParams` provenance block. Produced by [`resolve_dada_params`] so the
/// `dada`, `dada-pooled`, and `dada-pseudo` handlers all share one code path.
struct ResolvedDada {
    params: dada::DadaParams,
    nq: usize,
    run: DadaRunParams,
    /// The error model's own `params` block, when it carried one. Kept so a
    /// caller can rebuild the errfun the model was fitted with — needed by
    /// `dada-pseudo --reestimate-err-between-rounds`, which must re-fit round 2's
    /// model the same way the original was fitted.
    err_params: Option<LearnedErrParams>,
}

/// Resolve how much instrumentation to collect from the flags that ask for it.
///
/// `--verbose` keeps implying full attribution, so today's stderr output stays
/// byte-identical and the archived `docs/findings/data/*.txt` remain comparable
/// with new runs (issue #162). `--metrics-json` on its own collects only the
/// free tier, which is why it is safe to leave on for a production run.
fn resolve_measure_level(
    verbose: bool,
    metrics_json: Option<&PathBuf>,
    metrics_attribution: bool,
) -> MeasureLevel {
    if verbose || metrics_attribution {
        MeasureLevel::Attribution
    } else if metrics_json.is_some() {
        MeasureLevel::Phases
    } else {
        MeasureLevel::Off
    }
}

/// Load the error model and resolve every DADA parameter via the three-tier
/// precedence (CLI explicit > inherited from err-model `params` > built-in
/// default), warning when an explicit value differs from the model's.
///
/// `aux_outputs` and `pool` are handler-specific and passed in. `n_prior` is
/// filled in later by the caller (priors are marked after this point), so it is
/// left at 0 here.
#[allow(clippy::too_many_arguments)]
fn resolve_dada_params(
    error_model: &Path,
    use_err_in: bool,
    inherit_err_params: bool,
    threads: usize,
    aux_outputs: bool,
    pool: bool,
    verbose: bool,
    measure: MeasureLevel,
    denoise: cli::DadaDenoiseArgs,
    experimental: cli::ExperimentalArgs,
) -> io::Result<ResolvedDada> {
    let cli::DadaDenoiseArgs {
        omega_a,
        omega_c,
        omega_p,
        min_fold,
        min_hamming,
        min_abund,
        detect_singletons,
        max_clust,
        greedy,
        use_quals,
        band,
        gap_p,
        homo_gap_p,
        match_score,
        mismatch,
        align_backend,
        kdist_cutoff,
        kmer_size,
        no_kmer_screen,
    } = denoise;
    let cli::ExperimentalArgs {
        screen_backend,
        minimizer_k,
        minimizer_w,
        screen_audit,
        wfa_max_edits,
    } = experimental;
    let em: ErrorModelJson = read_tagged_json(error_model, &["learn-errors", "errors-from-sample"])
        .with_path(error_model)?;
    let nq = em.nq;
    let rows = if use_err_in { &em.err_in } else { &em.err_out };
    if rows.len() != 16 || rows.iter().any(|r| r.len() != nq) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Error model matrix must be 16 × {nq}, got {} rows",
                rows.len()
            ),
        ));
    }
    let err_mat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();

    if inherit_err_params && em.params.is_none() {
        eprintln!(
            "[dada] warning: --inherit-err-params requested but error model {} has no `params` block (likely produced by a pre-provenance version); falling back to built-in defaults",
            error_model.display(),
        );
    }

    let p = em.params.as_ref();
    macro_rules! resolve {
        ($cli:expr, $em_field:ident, $default:expr) => {{
            match ($cli, inherit_err_params, p) {
                (Some(v), _, _) => v,
                (None, true, Some(em_params)) => em_params.$em_field,
                _ => $default,
            }
        }};
    }
    let omega_a = resolve!(omega_a, omega_a, 1e-40);
    // `omega_c` is intentionally not inherited from the err model:
    // learn-errors uses 0 (R DADA2 convention), dada uses 1e-40.
    let omega_c = omega_c.unwrap_or(1e-40);
    let omega_p = resolve!(omega_p, omega_p, 1e-4);
    let min_fold = resolve!(min_fold, min_fold, 1.0);
    let min_hamming = resolve!(min_hamming, min_hamming, 1);
    let min_abund = resolve!(min_abund, min_abund, 1);
    let detect_singletons = resolve!(detect_singletons, detect_singletons, false);
    let band = resolve!(band, band, 16);
    let gap_p = resolve!(gap_p, gap_p, -8);
    // R normalizes a positive gap penalty to negative (dada.R:223) before the
    // homopolymer comparison; do the same so `--gap-p 8` == `-8`.
    let gap_p = if gap_p > 0 { -gap_p } else { gap_p };
    let match_score = resolve!(match_score, match_score, 5);
    let mismatch = resolve!(mismatch, mismatch, -4);
    // R's HOMOPOLYMER_GAP_PENALTY = NULL tracks GAP_PENALTY: default tier is
    // the resolved gap_p, not a literal -8. A positive value is normalized to
    // negative as well (dada.R:227), matching the gap_p handling above.
    let homo_gap_p = resolve!(homo_gap_p, homo_gap_p, gap_p);
    let homo_gap_p = if homo_gap_p > 0 {
        -homo_gap_p
    } else {
        homo_gap_p
    };
    note_homopolymer_gapping(verbose, gap_p, homo_gap_p);
    let max_clust = resolve!(max_clust, max_clust, 0);
    let greedy = resolve!(greedy, greedy, true);
    let use_quals = resolve!(use_quals, use_quals, true);
    // The cutoff default is per-BACKEND: 0.42 is calibrated for the frequency
    // vector and over-screens the sketch ~3x (see minimizers::MINIMIZER_KDIST_CUTOFF).
    // An explicit --kdist-cutoff, or an inherited one, still wins.
    let backend_default_cutoff = match screen_backend.unwrap_or_default() {
        ScreenBackend::Minimizer => minimizers::MINIMIZER_KDIST_CUTOFF,
        ScreenBackend::Kmer => 0.42,
    };
    let kdist_cutoff = resolve!(kdist_cutoff, kdist_cutoff, backend_default_cutoff);
    let kmer_size = resolve!(kmer_size, kmer_size, 5);
    let use_kmers = match (no_kmer_screen, inherit_err_params, p) {
        (Some(no), _, _) => !no,
        (None, true, Some(em_params)) => em_params.use_kmers,
        _ => true,
    };
    let backend = resolve!(align_backend, backend, AlignBackend::Nw);
    // WFA edit-budget cap default (issue #51): generous enough never to truncate
    // a real error-copy alignment (those are ~99.9% identical), tight enough to
    // bound divergent non-error-copy pairs that slip past the k-mer screen.
    let wfa_max_edits = resolve!(wfa_max_edits, wfa_max_edits, WFA_MAX_EDITS_DEFAULT);

    // ---- Consistency warnings (only when NOT inheriting) ----
    if !inherit_err_params && let Some(em_params) = p {
        let mut mismatches: Vec<String> = Vec::new();
        macro_rules! check {
            ($name:literal, $cli_val:expr, $em_val:expr) => {
                if $cli_val != $em_val {
                    mismatches.push(format!(
                        "  {} = {:?} (err model: {:?})",
                        $name, $cli_val, $em_val
                    ));
                }
            };
        }
        check!("omega_a", omega_a, em_params.omega_a);
        check!("omega_p", omega_p, em_params.omega_p);
        check!("min_fold", min_fold, em_params.min_fold);
        check!("min_hamming", min_hamming, em_params.min_hamming);
        check!("min_abund", min_abund, em_params.min_abund);
        check!(
            "detect_singletons",
            detect_singletons,
            em_params.detect_singletons
        );
        check!("band", band, em_params.band);
        check!("gap_p", gap_p, em_params.gap_p);
        check!("homo_gap_p", homo_gap_p, em_params.homo_gap_p);
        check!("kdist_cutoff", kdist_cutoff, em_params.kdist_cutoff);
        check!("kmer_size", kmer_size, em_params.kmer_size);
        check!("use_kmers", use_kmers, em_params.use_kmers);
        check!("align_backend", backend, em_params.backend);
        check!("wfa_max_edits", wfa_max_edits, em_params.wfa_max_edits);
        // The screen gates which pairs reach build_trans_mat, so it shapes the
        // fitted model as surely as kdist_cutoff does -- up to 90.9% relative
        // difference in err_out on the MiSeq SOP. Applying a model across a
        // screen change is a provenance error worth naming.
        check!(
            "screen_backend",
            screen_backend.unwrap_or_default(),
            em_params.screen_backend
        );
        if em_params.screen_backend == ScreenBackend::Minimizer
            || screen_backend.unwrap_or_default() == ScreenBackend::Minimizer
        {
            check!(
                "minimizer_k",
                minimizer_k.unwrap_or(minimizers::MINIMIZER_K),
                em_params.minimizer_k
            );
            check!(
                "minimizer_w",
                minimizer_w.unwrap_or(minimizers::MINIMIZER_W),
                em_params.minimizer_w
            );
        }
        if !mismatches.is_empty() {
            eprintln!(
                "[dada] warning: {} dada parameter(s) differ from error model {}; pass --inherit-err-params to adopt the err model's values:",
                mismatches.len(),
                error_model.display(),
            );
            for line in &mismatches {
                eprintln!("{line}");
            }
        }
    }

    let align_params = AlignParams {
        backend,
        wfa_max_edits,
        match_score,
        mismatch,
        gap_p,
        homo_gap_p,
        use_kmers,
        kdist_cutoff,
        kmer_size,
        // Deliberately NOT resolved through the error model's `params` block,
        // unlike every neighbour here. The screen backend is experimental and
        // absent from models written by any released version, so inheriting it
        // would silently resolve to `Kmer` and override an explicit
        // `--screen-backend minimizer`. Consequence: a model learned under one
        // screen and applied under the other is not flagged the way a
        // `kdist_cutoff` mismatch is. Revisit if the backend is promoted.
        screen_backend: screen_backend.unwrap_or_default(),
        minimizer_k: minimizer_k.unwrap_or(minimizers::MINIMIZER_K),
        minimizer_w: minimizer_w.unwrap_or(minimizers::MINIMIZER_W),
        screen_audit,
        band,
        vectorized: true,
        gapless: true,
    };

    if verbose {
        eprintln!("[dada] {}", nwalign::backend_repr(&align_params));
        let (alloc, alloc_warn) = misc::cpu_allocation_repr(threads);
        eprintln!("[dada] cpu allocation: {alloc}");
        if let Some(w) = alloc_warn {
            eprintln!("[dada] {w}");
        }
        // Resolved gate values, not an echo of the environment (#145).
        for line in dada2_rs::gates::report() {
            eprintln!("[dada] {line}");
        }
    }

    // `progress_tag` stays unset here: the concurrent paths set it per sample
    // (see denoise_and_serialize).
    let params = dada::DadaParams {
        aux_outputs,
        ..dada::DadaParams::new(
            align_params,
            err_mat,
            nq,
            dada::DenoiseOpts {
                omega_a,
                omega_c,
                omega_p,
                detect_singletons,
                max_clust,
                min_fold,
                min_hamming,
                min_abund,
                use_quals,
                greedy,
            },
            threads,
            verbose,
            measure,
        )
    };

    let run = DadaRunParams {
        screen_backend: align_params.screen_backend,
        minimizer_k: match align_params.screen_backend {
            ScreenBackend::Minimizer => Some(align_params.minimizer_k),
            ScreenBackend::Kmer => None,
        },
        minimizer_w: match align_params.screen_backend {
            ScreenBackend::Minimizer => Some(align_params.minimizer_w),
            ScreenBackend::Kmer => None,
        },
        omega_a,
        omega_c,
        omega_p,
        min_fold,
        min_hamming,
        min_abund,
        detect_singletons,
        band,
        homo_gap_p,
        gap_p,
        match_score,
        mismatch,
        max_clust,
        greedy,
        use_quals,
        kdist_cutoff,
        kmer_size,
        use_kmers,
        pool,
        n_prior: 0,
        // Set only by dada-pseudo, after round 1, if the re-fit actually ran.
        reestimated_err_between_rounds: false,
        pool_tiebreak: None,
    };

    Ok(ResolvedDada {
        params,
        nq,
        run,
        err_params: em.params,
    })
}

/// Fold one round-1 sample's transition counts into the shared accumulator
/// (R's `accumulateTrans`). A no-op unless `aux_outputs` was set, i.e. unless
/// `dada-pseudo --reestimate-err-between-rounds` asked for the re-fit.
///
/// Integer addition, so the fold is exact and independent of the order samples
/// complete in — the accumulated matrix does not vary with `--sample-jobs`.
fn accumulate_round1_trans(
    acc: &std::sync::Mutex<Vec<u32>>,
    result: &dada::DadaResult,
) -> io::Result<()> {
    let Some(aux) = result.aux.as_ref() else {
        return Ok(());
    };
    let mut acc = acc.lock().unwrap();
    if aux.transitions.len() != acc.len() {
        // Unreachable in practice: compute_aux sizes the matrix from
        // `params.err_ncol` (the same nq the accumulator was built from) whenever
        // quals are in use, and the caller rejects --use-quals false up front.
        // Checked anyway so a future change to that sizing fails loudly instead
        // of silently truncating the fit's input.
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "round-1 transition matrix is {} entries, expected {}",
                aux.transitions.len(),
                acc.len()
            ),
        ));
    }
    for (slot, add) in acc.iter_mut().zip(&aux.transitions) {
        *slot += add;
    }
    Ok(())
}

/// Rebuild the [`ErrFun`] an error model was fitted with, from the `params`
/// block the model carries.
///
/// Used by `dada-pseudo --reestimate-err-between-rounds`, which must re-fit
/// round 2's model the same way the supplied one was fitted. R applies a fixed
/// `errorEstimationFunction` (default `loessErrfun`) independent of how `err`
/// was produced; taking the recorded errfun is the closest analogue that also
/// keeps the loess surface and clamps consistent with the original fit.
fn errfun_from_learned(p: &LearnedErrParams) -> io::Result<ErrFun> {
    let config = p.loess.as_ref().map(LoessConfig::from).unwrap_or_default();
    let missing = |what: &str| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "error model records errfun \"{}\" but no {what}; cannot re-fit \
                 between pseudo rounds",
                p.errfun
            ),
        )
    };
    Ok(match p.errfun.as_str() {
        "loess" => ErrFun::Loess { config },
        "noqual" => ErrFun::Noqual {
            pseudocount: p.errfun_pseudocount.unwrap_or(1.0),
            config,
        },
        "binned-qual" => ErrFun::BinnedQual {
            bins: p.errfun_bins.clone().ok_or_else(|| missing("bins"))?,
            config,
        },
        "pacbio" => ErrFun::PacBio { config },
        "external" => ErrFun::External {
            command: p.errfun_cmd.clone().ok_or_else(|| missing("command"))?,
        },
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("error model records unknown errfun \"{other}\""),
            ));
        }
    })
}

/// Per-sample round-1 ASVs tagged with their sample index (collected unordered
/// from the concurrent round-1 denoise, then reassembled by index).
type IndexedAsvs = Vec<(usize, Vec<(String, u32)>)>;

/// Like [`IndexedAsvs`] but also carrying the resolved sample name — used by the
/// streaming dada-pseudo round 1, which loads names on the fly (no pre-load).
type IndexedNamedAsvs = Vec<(usize, String, Vec<(String, u32)>)>;

/// Run `f` over samples `0..n` with bounded across-sample concurrency: spawn
/// `jobs` workers, each owning a rayon sub-pool of ~`threads / jobs` threads,
/// pulling sample indices from a shared counter. `f` gets the worker's sub-pool
/// so it can `sub_pool.install(|| dada_uniques(...))`, pinning the per-sample
/// comparison map to that sub-pool. A single sample's map is often too small to
/// feed many threads, so this keeps every core fed. `jobs == 1` reproduces the
/// serial, full-pool behavior exactly. Returns the first error, if any.
fn for_each_sample_concurrent(
    n: usize,
    jobs: usize,
    threads: usize,
    f: impl Fn(usize, &rayon::ThreadPool) -> io::Result<()> + Sync,
) -> io::Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering};

    let jobs = jobs.clamp(1, n.max(1));
    // Split `threads` across `jobs` sub-pools, spreading the remainder so the
    // pools differ by at most one thread (e.g. 20 threads / 3 jobs -> 7,7,6).
    let base = (threads / jobs).max(1);
    let rem = threads % jobs;
    let pools: Vec<rayon::ThreadPool> = (0..jobs)
        .map(|j| {
            let t = if j < rem { base + 1 } else { base };
            rayon::ThreadPoolBuilder::new().num_threads(t).build()
        })
        .collect::<Result<_, _>>()
        .map_err(io::Error::other)?;

    let next = AtomicUsize::new(0);
    let err: Mutex<Option<io::Error>> = Mutex::new(None);
    std::thread::scope(|scope| {
        for pool in &pools {
            let (next, err, f) = (&next, &err, &f);
            scope.spawn(move || {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= n || err.lock().unwrap().is_some() {
                        break;
                    }
                    if let Err(e) = f(i, pool) {
                        let mut slot = err.lock().unwrap();
                        slot.get_or_insert(e);
                        break;
                    }
                }
            });
        }
    });
    match err.into_inner().unwrap() {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

/// Load one sample's dereplicated uniques as `RawInput`s (FASTQ or derep/sample
/// JSON), returning them alongside any JSON-embedded sample name. Errors if the
/// sample has no uniques. Used by the multi-sample `dada`/`dada-pseudo` paths.
fn load_sample_raws(
    path: &Path,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
) -> io::Result<(Vec<dada::RawInput>, Option<String>)> {
    let (derep, json_sample) = load_derep_for_dada(
        path,
        phred_offset,
        pool,
        verbose,
        &mut DerepLoadCost::default(),
    )?;
    let raws: Vec<dada::RawInput> = derep
        .uniques
        .into_iter()
        .zip(derep.quals)
        .map(|((seq, count), quals)| {
            let sequence = String::from_utf8(seq)
                .unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned());
            dada::RawInput {
                seq: sequence,
                abundance: count as u32,
                prior: false,
                quals: Some(quals),
            }
        })
        .collect();
    if raws.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{}: no uniques found", path.display()),
        ));
    }
    Ok((raws, json_sample))
}

/// Map a DADA result's clusters to (decoded sequence, reads) ASV pairs.
fn result_to_asvs(result: &dada::DadaResult) -> Vec<(String, u32)> {
    result
        .clusters
        .iter()
        .map(|c| {
            let sequence: String = c
                .sequence
                .iter()
                .map(|&b| misc::nt_decode(b) as char)
                .collect();
            (sequence, c.abundance)
        })
        .collect()
}

/// Mark each RawInput whose (uppercased) sequence is in `prior_set` as a prior.
/// Returns the number marked.
fn mark_priors(
    raws: &mut [dada::RawInput],
    prior_set: &std::collections::HashSet<String>,
) -> usize {
    let mut n = 0usize;
    for inp in raws.iter_mut() {
        if prior_set.contains(&inp.seq.to_ascii_uppercase()) {
            inp.prior = true;
            n += 1;
        }
    }
    n
}

/// One sample's output JSON from `dada`, `dada-pooled` or `dada-pseudo`.
#[derive(Serialize)]
struct DadaOutput {
    sample: String,
    /// Original input file name (no directory) for provenance.
    input_file: String,
    num_asvs: usize,
    total_reads: u32,
    asvs: Vec<AsvEntry>,
    stats: DadaStats,
    params: DadaRunParams,
    map: Vec<Option<usize>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    aux: Option<AuxJson>,
    /// This sample's position in `dada-pooled`'s input order, which decides
    /// tied uniques under `--pool-tiebreak first-seen` (#260); pooled only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pool_input_index: Option<usize>,
}

#[derive(Serialize)]
struct ClusterStatJson {
    sequence: String,
    abundance: u32,
    n0: u32,
    n1: u32,
    nunq: u32,
    pval: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    birth_from: Option<usize>,
    birth_pval: f64,
    birth_fold: f64,
    birth_ham: u32,
    birth_e: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    birth_qave: Option<f64>,
}

#[derive(Serialize)]
struct BirthSubJson {
    cluster: usize,
    pos: u16,
    nt0: char,
    nt1: char,
    #[serde(skip_serializing_if = "Option::is_none")]
    qual: Option<u8>,
}

/// `dada --aux-outputs`: R-parity per-cluster diagnostics.
#[derive(Serialize)]
struct AuxJson {
    cluster_stats: Vec<ClusterStatJson>,
    cluster_quality: Vec<Vec<f64>>,
    cluster_quality_maxlen: usize,
    birth_subs: Vec<BirthSubJson>,
    transitions: Vec<u32>,
    transitions_ncol: usize,
}

impl From<&dada::DadaAux> for AuxJson {
    fn from(a: &dada::DadaAux) -> Self {
        let cluster_stats_j = a
            .cluster_stats
            .iter()
            .map(|c| {
                let sequence: String = c
                    .sequence
                    .iter()
                    .map(|&b| misc::nt_decode(b) as char)
                    .collect();
                ClusterStatJson {
                    sequence,
                    abundance: c.abundance,
                    n0: c.n0,
                    n1: c.n1,
                    nunq: c.nunq,
                    pval: c.pval,
                    birth_from: c.birth_from,
                    birth_pval: c.birth_pval,
                    birth_fold: c.birth_fold,
                    birth_ham: c.birth_ham,
                    birth_e: c.birth_e,
                    birth_qave: c.birth_qave,
                }
            })
            .collect();
        let birth_subs_j = a
            .birth_subs
            .iter()
            .map(|r| BirthSubJson {
                cluster: r.cluster,
                pos: r.pos,
                nt0: r.nt0 as char,
                nt1: r.nt1 as char,
                qual: r.qual,
            })
            .collect();
        AuxJson {
            cluster_stats: cluster_stats_j,
            cluster_quality: a.cluster_quality.clone(),
            cluster_quality_maxlen: a.cluster_quality_maxlen,
            birth_subs: birth_subs_j,
            transitions: a.transitions.clone(),
            transitions_ncol: a.transitions_ncol,
        }
    }
}

/// The `--prior` FASTA as a set of uppercased sequences.
fn read_prior_set(path: &Path) -> io::Result<std::collections::HashSet<String>> {
    Ok(read_fasta_records(path)
        .with_path(path)?
        .into_iter()
        .map(|(_, seq)| String::from_utf8_lossy(&seq).to_ascii_uppercase())
        .collect())
}

/// A sample's output name: the one given on the command line, else the one
/// embedded in a derep/sample JSON input, else the FASTQ file stem.
fn sample_label(cli: Option<String>, json: Option<String>, path: &Path) -> String {
    cli.or(json).unwrap_or_else(|| fastq_stem(path))
}

/// Reject a `--sample-names` list that does not name every input exactly once.
fn check_sample_names(names: Option<&[String]>, n_samples: usize) -> io::Result<()> {
    match names {
        Some(names) if names.len() != n_samples => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--sample-names has {} entries but {} input file(s) given",
                names.len(),
                n_samples
            ),
        )),
        _ => Ok(()),
    }
}

/// How many samples to denoise at once, each on a `threads / jobs` sub-pool.
/// Samples are independent and single-pass, so this keeps every core fed and
/// bounds memory to `jobs` samples in flight. The default, round(threads/4) ≈
/// 4 threads per sample, is where the sample-jobs sweep's wall-time curve
/// plateaus: more samples in flight fill cores better than the single-sample
/// efficiency curve suggests. At <= 4 threads it is 1, the serial path.
fn concurrent_samples(sample_jobs: Option<usize>, threads: usize, n_samples: usize) -> usize {
    sample_jobs
        .unwrap_or_else(|| ((threads as f64 / 4.0).round() as usize).max(1))
        .clamp(1, n_samples.max(1))
}

/// One sample's entry in a `--metrics-json` document: name, round, metrics.
type MetricRun = (String, Option<u8>, metrics::RunMetrics);

/// What `denoise_and_serialize` returns for one sample.
type SampleOutput = (
    String,
    Vec<failed_uniques::Row>,
    Option<metrics::RunMetrics>,
);

/// The outputs of a run that denoises samples independently (`dada` with
/// several inputs, `dada-pseudo` round 2): one `{sample}.json[.gz]` each,
/// written as samples finish, then the failed-uniques TSV and the metrics
/// document for the whole run.
struct SampleSink<'a> {
    tag: &'static str,
    output_dir: &'a Path,
    gzip: bool,
    verbose: bool,
    /// The metrics `round` every sample is recorded under.
    round: Option<u8>,
    failed: Mutex<Vec<failed_uniques::Row>>,
    // Samples finish out of order under `--sample-jobs`, so this is sorted
    // before it is written.
    metrics: Mutex<Vec<MetricRun>>,
}

impl<'a> SampleSink<'a> {
    fn new(
        tag: &'static str,
        output_dir: &'a Path,
        gzip: bool,
        verbose: bool,
        round: Option<u8>,
    ) -> Self {
        SampleSink {
            tag,
            output_dir,
            gzip,
            verbose,
            round,
            failed: Mutex::new(Vec::new()),
            metrics: Mutex::new(Vec::new()),
        }
    }

    /// Keep one sample's failed rows and metrics, and write its JSON.
    fn write(&self, sample: &str, (json, failed, run_metrics): SampleOutput) -> io::Result<()> {
        self.failed.lock().unwrap().extend(failed);
        if let Some(m) = run_metrics {
            self.metrics
                .lock()
                .unwrap()
                .push((sample.to_string(), self.round, m));
        }
        let out_path = self.output_dir.join(if self.gzip {
            format!("{sample}.json.gz")
        } else {
            format!("{sample}.json")
        });
        misc::write_maybe_gz(&out_path, json.as_bytes())?;
        if self.verbose {
            eprintln!("[{}] wrote {}", self.tag, out_path.display());
        }
        Ok(())
    }

    /// Write the run-level failed-uniques TSV and metrics document, if asked.
    fn finish(
        self,
        failed_uniques: Option<&Path>,
        metrics_json: Option<&Path>,
        t_start: std::time::Instant,
        measure_level: MeasureLevel,
    ) -> io::Result<()> {
        if let Some(path) = failed_uniques {
            let rows = self.failed.into_inner().unwrap();
            write_failed_uniques(self.tag, path, rows, self.verbose)?;
        }
        if let Some(path) = metrics_json {
            let mut runs = self.metrics.into_inner().unwrap();
            // Concurrent samples finish out of order; sort so two runs
            // of the same inputs produce byte-identical documents.
            runs.sort_by(|a, b| a.0.cmp(&b.0));
            write_run_metrics(self.tag, path, runs, t_start, measure_level, self.verbose)?;
        }
        Ok(())
    }
}

/// Write the `--failed-uniques` TSV (issue #60).
fn write_failed_uniques(
    tag: &str,
    path: &Path,
    rows: Vec<failed_uniques::Row>,
    verbose: bool,
) -> io::Result<()> {
    let n = failed_uniques::write_tsv(path, rows)?;
    if verbose {
        eprintln!(
            "[{tag}] wrote {n} failed-unique row(s) to {}",
            path.display()
        );
    }
    Ok(())
}

/// Write a `--metrics-json` document whose only pipeline phase is `dada`.
fn write_run_metrics(
    tag: &str,
    path: &Path,
    runs: Vec<MetricRun>,
    t_start: std::time::Instant,
    measure_level: MeasureLevel,
    verbose: bool,
) -> io::Result<()> {
    let mut doc = MetricsDocument::new(t_start.elapsed(), measure_level);
    doc.pipeline.dada = Some(t_start.elapsed().as_secs_f64());
    for (sample, round, m) in runs {
        doc.push(sample, round, m);
    }
    write_metrics_json(path, &doc)?;
    if verbose {
        eprintln!("[{tag}] wrote run metrics to {}", path.display());
    }
    Ok(())
}

/// Flag one sample's round-2 priors for `dada-pseudo`.
fn mark_round2_priors(
    raws: &mut [dada::RawInput],
    prior_set: &std::collections::HashSet<String>,
    sample_name: &str,
    verbose: bool,
) {
    let n_marked = mark_priors(raws, prior_set);
    if verbose {
        eprintln!(
            "[dada-pseudo]   round 2 {sample_name}: {n_marked} of {} unique(s) flagged as prior",
            raws.len(),
        );
    }
}

/// A `--cluster-trace` request: where to write it and what to include.
struct TraceRequest<'a> {
    path: &'a Path,
    params: cluster_trace::TraceParams,
}

/// Write the cluster trace of a finished denoising run.
#[allow(clippy::too_many_arguments)]
fn write_cluster_trace(
    tag: &str,
    trace: &TraceRequest,
    sample: &str,
    raw_inputs: &[dada::RawInput],
    result: &dada::DadaResult,
    params: &dada::DadaParams,
    compact: bool,
    verbose: bool,
) -> io::Result<()> {
    if let Some(parent) = trace.path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    cluster_trace::write_trace(
        trace.path,
        sample,
        None, // no iteration: this is the final dada run
        raw_inputs,
        result,
        Some(&params.err_mat),
        params.err_ncol,
        trace.params,
        compact,
    )?;
    if verbose {
        eprintln!("[{tag}] cluster trace written to {}", trace.path.display());
    }
    Ok(())
}

/// One ASV record in a dada-family output JSON. Shared by the single-input
/// `dada`, `dada-pooled`, and the multi-sample `dada` / `dada-pseudo` paths.
#[derive(Serialize)]
struct AsvEntry {
    sequence: String,
    abundance: u32,
    birth_type: String,
    birth_pval: f64,
    birth_fold: f64,
    birth_e: f64,
}

#[derive(Serialize)]
struct DadaStats {
    nalign: u64,
    nshroud: u64,
}

/// Stable string label for a cluster's birth type (matches R DADA2's `$type`).
fn birth_type_str(bt: &BirthType) -> &'static str {
    match bt {
        BirthType::Initial => "Initial",
        BirthType::Abundance => "Abundance",
        BirthType::Prior => "Prior",
        BirthType::Singleton => "Singleton",
    }
}

/// Build an [`AsvEntry`] from a cluster summary, decoding its center sequence.
/// `abundance` is passed explicitly because the pooled per-sample path uses a
/// recomputed per-sample read count rather than `cluster.reads`.
fn asv_entry_from_cluster(cluster: &dada::ClusterSummary, abundance: u32) -> AsvEntry {
    AsvEntry {
        sequence: cluster
            .sequence
            .iter()
            .map(|&b| misc::nt_decode(b) as char)
            .collect(),
        abundance,
        birth_type: birth_type_str(&cluster.birth_type).to_string(),
        birth_pval: cluster.birth_pval,
        birth_fold: cluster.birth_fold,
        birth_e: cluster.birth_e,
    }
}

/// Serialize a value to JSON, compact or pretty per `compact`.
/// Write the `--metrics-json` document. Always pretty-printed and always
/// gzip-free: this is a small file read by humans and by `dev/` scripts, not a
/// bulk artifact, and `jq` over a plain file is the point of it (issue #162).
fn write_metrics_json(path: &Path, doc: &metrics::MetricsDocument) -> io::Result<()> {
    let mut body = serde_json::to_string_pretty(doc).map_err(io::Error::other)?;
    body.push('\n');
    std::fs::write(path, body)
}

fn to_json<T: Serialize>(value: &T, compact: bool) -> io::Result<String> {
    if compact {
        serde_json::to_string(value)
    } else {
        serde_json::to_string_pretty(value)
    }
    .map_err(io::Error::other)
}

#[allow(clippy::too_many_arguments)]
/// Denoise one sample and serialize its dada JSON, writing its cluster trace
/// first when `trace` asks for one. When `collect_failed` is set,
/// also returns the uniques that failed to denoise (`map == null`) as
/// [`failed_uniques::Row`]s tagged with `sample`; otherwise the row vec is empty.
fn denoise_and_serialize(
    tag: &'static str,
    sample: &str,
    input_file: &str,
    raw_inputs: &[dada::RawInput],
    params: &dada::DadaParams,
    run_params: &DadaRunParams,
    pool: &rayon::ThreadPool,
    compact: bool,
    collect_failed: bool,
    verbose: bool,
    // Sample label for the bud-round progress records. `Some` only when samples
    // are denoised concurrently, where the records would otherwise be
    // unattributable; `None` keeps the text byte-identical to R's (#172).
    progress_tag: Option<String>,
    trace: Option<&TraceRequest>,
) -> io::Result<SampleOutput> {
    // Only clone when a tag is actually needed: the copy carries `err_mat`,
    // and a serial run has nothing to disambiguate anyway.
    let tagged;
    let params = match progress_tag {
        Some(tag) => {
            tagged = dada::DadaParams {
                progress_tag: Some(tag),
                ..params.clone()
            };
            &tagged
        }
        None => params,
    };
    let mut result = pool
        .install(|| dada::dada_uniques(raw_inputs, params))
        .map_err(io::Error::other)?;
    let run_metrics = result.metrics.take();

    if verbose {
        eprintln!(
            "[{tag}] {} ASV(s) from {} unique input(s); {} aligns, {} shrouded",
            result.clusters.len(),
            raw_inputs.len(),
            result.nalign,
            result.nshroud,
        );
    }
    if let Some(t) = trace {
        write_cluster_trace(
            tag, t, sample, raw_inputs, &result, params, compact, verbose,
        )?;
    }

    let total_reads: u32 = result.clusters.iter().map(|c| c.abundance).sum();
    let asvs: Vec<AsvEntry> = result
        .clusters
        .iter()
        .map(|c| asv_entry_from_cluster(c, c.abundance))
        .collect();

    let mut run_params = *run_params;
    run_params.n_prior = raw_inputs.iter().filter(|r| r.prior).count();

    let failed = if collect_failed {
        result
            .map
            .iter()
            .enumerate()
            .filter(|(_, m)| m.is_none())
            .map(|(i, _)| failed_uniques::Row {
                sequence: raw_inputs[i].seq.clone(),
                sample: sample.to_string(),
                reads: raw_inputs[i].abundance,
            })
            .collect()
    } else {
        Vec::new()
    };

    let out = DadaOutput {
        sample: sample.to_string(),
        input_file: input_file.to_string(),
        num_asvs: asvs.len(),
        total_reads,
        asvs,
        stats: DadaStats {
            nalign: result.nalign,
            nshroud: result.nshroud,
        },
        params: run_params,
        aux: result.aux.as_ref().map(AuxJson::from),
        pool_input_index: None,
        map: result.map,
    };

    Ok((
        to_json(&Tagged::new(tag, out), compact)?,
        failed,
        run_metrics,
    ))
}

/// Where a pooled run's serial load front spends its time (issue #127).
///
/// `derep` is serial by design (issue #41 streams one sample at a time to hold
/// the pooled memory peak down), so its share of wall grows with every thread
/// added to `run_dada` — 7% of pooled wall at 24 threads, 12.5% at 128. This
/// splits it finely enough to tell a filesystem problem from a format one.
#[derive(Default)]
struct DerepLoadCost {
    /// I/O + gunzip (JSON inputs only).
    read: std::time::Duration,
    /// serde deserialization (JSON inputs only).
    parse: std::time::Duration,
    /// Sorting + conversion to `Derep` (JSON), or the whole streaming
    /// dereplication (FASTQ, where the stages cannot be separated).
    build: std::time::Duration,
    /// Bytes loaded, for a throughput figure: uncompressed for JSON, on-disk
    /// (so compressed, for `.gz`) for FASTQ, which is streamed through the decoder.
    bytes: u64,
    /// Per-sample totals, for the straggler spread.
    per_sample: Vec<std::time::Duration>,
    /// True if any input was FASTQ, whose `build` is not comparable to JSON's.
    any_fastq: bool,
    /// Inputs that were JSON, to tell an all-FASTQ run from a mixed one.
    n_json: usize,
}

/// Build a [`derep::Derep`] for `dada` / `dada-pooled` from either a FASTQ file
/// (uncompressed or gzipped) or a derep/sample JSON file.
///
/// JSON inputs are defensively sorted by abundance descending — DADA2 assumes
/// the most-abundant unique is at index 0.  The `map` (read → unique) field is
/// only populated from the FASTQ path; JSON inputs leave it empty since neither
/// `dada` nor `dada-pooled` consult it.
///
/// Returns the dereplicated table plus the JSON's embedded `sample` field
/// when present; FASTQ inputs always return `None` for the name.
fn load_derep_for_dada(
    path: &Path,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
    cost: &mut DerepLoadCost,
) -> io::Result<(derep::Derep, Option<String>)> {
    let t_sample = std::time::Instant::now();
    let r = load_derep_for_dada_inner(path, phred_offset, pool, verbose, cost);
    cost.per_sample.push(t_sample.elapsed());
    r
}

fn load_derep_for_dada_inner(
    path: &Path,
    phred_offset: u8,
    pool: &rayon::ThreadPool,
    verbose: bool,
    cost: &mut DerepLoadCost,
) -> io::Result<(derep::Derep, Option<String>)> {
    if is_json_path(path) {
        #[derive(serde::Deserialize)]
        struct UniqueEntryJson {
            sequence: String,
            count: u64,
            /// Per-position integer Phred SUM; mean recovered as sum/count on demand.
            qual_sum: Vec<u32>,
        }
        #[derive(serde::Deserialize)]
        struct SampleJson {
            /// Carried on the struct so the tag is validated from the *same*
            /// parse that produces the data. Checking it separately means a
            /// second full scan of the document, which on this path (hundreds
            /// of MB per pooled run) costs 37% of parse time — see #133.
            #[serde(default)]
            dada2_rs_command: Option<String>,
            #[serde(default)]
            sample: Option<String>,
            #[serde(default)]
            sort_order: Option<String>,
            uniques: Vec<UniqueEntryJson>,
        }
        let mut jc = misc::JsonReadCost::default();
        let parsed: SampleJson =
            misc::read_json_timed(path, &["derep", "sample"], &mut jc).with_path(path)?;
        misc::check_json_tag(
            path,
            parsed.dada2_rs_command.as_deref(),
            &["derep", "sample"],
        )
        .with_path(path)?;
        cost.read += jc.read;
        cost.parse += jc.parse;
        cost.bytes += jc.bytes;
        cost.n_json += 1;
        let t_build = std::time::Instant::now();
        let sample_name = parsed.sample;
        let mut entries = parsed.uniques;
        // Skip the defensive sort when the producer has declared the order.
        // Older JSONs without `sort_order` get sorted, matching prior behaviour.
        if parsed.sort_order.as_deref() != Some("abundance_desc") {
            entries.sort_by_key(|a| std::cmp::Reverse(a.count));
        }
        let mut uniques: Vec<(Vec<u8>, u64)> = Vec::with_capacity(entries.len());
        let mut quals: Vec<Vec<u32>> = Vec::with_capacity(entries.len());
        for u in entries {
            if !u.qual_sum.is_empty() && u.qual_sum.len() != u.sequence.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}: qual_sum length {} != sequence length {}",
                        path.display(),
                        u.qual_sum.len(),
                        u.sequence.len(),
                    ),
                ));
            }
            uniques.push((u.sequence.into_bytes(), u.count));
            quals.push(u.qual_sum);
        }
        cost.build += t_build.elapsed();
        Ok((
            derep::Derep {
                uniques,
                quals,
                map: Vec::new(),
            },
            sample_name,
        ))
    } else if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        // FASTQ is a single streaming pass — read, decompress and dereplicate
        // are interleaved and cannot be attributed separately, so the whole
        // cost lands in `build`.
        cost.any_fastq = true;
        cost.bytes += std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let t_build = std::time::Instant::now();
        let derep = dereplicate(
            MultiGzDecoder::new(File::open(path).with_path(path)?),
            phred_offset,
            pool,
            verbose,
        )?;
        cost.build += t_build.elapsed();
        Ok((derep, None))
    } else {
        cost.any_fastq = true;
        cost.bytes += std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let t_build = std::time::Instant::now();
        let derep = dereplicate(
            File::open(path).with_path(path)?,
            phred_offset,
            pool,
            verbose,
        )?;
        cost.build += t_build.elapsed();
        Ok((derep, None))
    }
}

impl DerepLoadCost {
    /// `--metrics-json` form of [`Self::report`]. `None` when nothing was loaded.
    fn detail(&self) -> Option<metrics::DerepDetail> {
        let n = self.per_sample.len();
        if n == 0 {
            return None;
        }
        let mut sorted = self.per_sample.clone();
        sorted.sort_unstable();
        let mb = self.bytes as f64 / 1_048_576.0;
        let rate = |d: std::time::Duration| (d.as_secs_f64() > 0.0).then(|| mb / d.as_secs_f64());
        let (input_kind, mb_per_s) = match (self.any_fastq, self.n_json > 0) {
            (false, _) => ("json", rate(self.read)),
            (true, false) => ("fastq", rate(self.build)),
            (true, true) => ("mixed", None),
        };
        let json = (self.n_json > 0).then_some(());
        Some(metrics::DerepDetail {
            input_kind: input_kind.to_string(),
            samples: n,
            bytes: self.bytes,
            read: json.map(|_| self.read.as_secs_f64()),
            parse: json.map(|_| self.parse.as_secs_f64()),
            build: self.build.as_secs_f64(),
            mb_per_s,
            per_sample: metrics::SpreadSecs {
                min: sorted[0].as_secs_f64(),
                median: sorted[n / 2].as_secs_f64(),
                max: sorted[n - 1].as_secs_f64(),
            },
        })
    }

    /// `--verbose` report, mirroring `compare split`. Returns `None` when
    /// nothing was loaded.
    ///
    /// The throughput figure is the point of this: it separates a filesystem
    /// problem (network storage, cold cache) from a format one (gzip + JSON),
    /// which have different fixes. Note `read` is uncompressed bytes over
    /// wall, so it understates the network rate for compressed inputs.
    fn report(&self, total: std::time::Duration) -> Option<String> {
        let n = self.per_sample.len();
        if n == 0 {
            return None;
        }
        let secs = total.as_secs_f64();
        let pct = |d: std::time::Duration| {
            if secs > 0.0 {
                100.0 * d.as_secs_f64() / secs
            } else {
                0.0
            }
        };
        let mb = self.bytes as f64 / 1_048_576.0;
        let mut sorted = self.per_sample.clone();
        sorted.sort_unstable();
        let ms = |d: std::time::Duration| d.as_secs_f64() * 1e3;

        let bytes_kind = match (self.any_fastq, self.n_json > 0) {
            (false, _) => "uncompressed",
            (true, false) => "on disk",
            (true, true) => "JSON uncompressed + FASTQ on disk",
        };
        let mut out = format!(
            "[dada-pooled] derep split (of {secs:.2}s over {n} sample(s), {mb:.0} MB {bytes_kind}):\n"
        );
        if !self.any_fastq {
            out.push_str(&format!(
                "[dada-pooled]   read+gunzip  {:7.2}s ({:4.1}%)   {:.0} MB/s\n",
                self.read.as_secs_f64(),
                pct(self.read),
                if self.read.as_secs_f64() > 0.0 {
                    mb / self.read.as_secs_f64()
                } else {
                    0.0
                },
            ));
            out.push_str(&format!(
                "[dada-pooled]   parse        {:7.2}s ({:4.1}%)   (serde_json)\n",
                self.parse.as_secs_f64(),
                pct(self.parse),
            ));
            out.push_str(&format!(
                "[dada-pooled]   build        {:7.2}s ({:4.1}%)   (sort + convert)\n",
                self.build.as_secs_f64(),
                pct(self.build),
            ));
        } else {
            out.push_str(&format!(
                "[dada-pooled]   dereplicate  {:7.2}s ({:4.1}%)   {:.0} MB/s  (FASTQ: read+decompress+derep in one pass)\n",
                self.build.as_secs_f64(),
                pct(self.build),
                if self.build.as_secs_f64() > 0.0 {
                    mb / self.build.as_secs_f64()
                } else {
                    0.0
                },
            ));
        }
        out.push_str(&format!(
            "[dada-pooled]   per-sample   min {:.0}ms  median {:.0}ms  max {:.0}ms",
            ms(sorted[0]),
            ms(sorted[n / 2]),
            ms(sorted[n - 1]),
        ));
        Some(out)
    }
}
