//! Error models: `learn-errors`, `errors-from-sample`, `kdist-calibrate`.

use std::{io, path::PathBuf};

use dada2_rs::{
    cli, cluster_trace, dada, error_models, kdist_calibrate, learn_errors, metrics, minimizers,
    misc, nwalign,
};

use learn_errors::{
    ErrFun, LearnDiagOptions, LearnedErrParams, learn_errors, load_derep_samples,
    load_fastq_samples,
};
use metrics::MeasureLevel;
use misc::Tagged;
use nwalign::AlignParams;
use serde::Serialize;

use super::common::{WFA_MAX_EDITS_DEFAULT, check_input_paths, note_homopolymer_gapping};
use error_models::{LoessConfig, LoessSurface};

pub(crate) fn run_learn_errors(args: cli::LearnErrorsArgs) -> io::Result<()> {
    let cli::LearnErrorsArgs {
        input,
        nbases,
        randomize,
        seed,
        phred_offset,
        fit,
        denoise,
        experimental,
        threads,
        output,
        compact,
        diag,
        verbose,
    } = args;
    check_input_paths("input", &input)?;
    let resolved = resolve_learn_params(&fit, denoise, experimental, threads, verbose)?;
    let max_consist = fit.max_consist;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;
    let all_inputs = load_fastq_samples(
        &input,
        nbases,
        randomize,
        seed,
        phred_offset,
        &pool,
        verbose,
    )?;

    learn_and_write(
        "learn-errors",
        all_inputs,
        resolved,
        max_consist,
        diag,
        &pool,
        verbose,
        compact,
        output,
    )
}

pub(crate) fn run_errors_from_sample(args: cli::ErrorsFromSampleArgs) -> io::Result<()> {
    let cli::ErrorsFromSampleArgs {
        input,
        fit,
        denoise,
        experimental,
        threads,
        output,
        compact,
        diag,
        verbose,
    } = args;
    check_input_paths("input", &input)?;
    let resolved = resolve_learn_params(&fit, denoise, experimental, threads, verbose)?;
    let max_consist = fit.max_consist;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?;

    let all_inputs = load_derep_samples(&input)?;

    if verbose {
        eprintln!(
            "[errors-from-sample] loaded {} sample(s) from JSON",
            all_inputs.len()
        );
    }

    learn_and_write(
        "errors-from-sample",
        all_inputs,
        resolved,
        max_consist,
        diag,
        &pool,
        verbose,
        compact,
        output,
    )
}

pub(crate) fn run_kdist_calibrate(args: cli::KdistCalibrateArgs) -> io::Result<()> {
    let cli::KdistCalibrateArgs {
        inputs,
        k,
        screen_backend,
        minimizer_k,
        minimizer_w,
        cutoff,
        leak_pct,
        band,
        max_pairs,
        max_uniques,
        per_sample,
        nearest_parent,
        from_dada,
        from_dada_pooled,
        derive_cutoff,
        derive_uniform_pairs,
        derep_dir,
        threads,
        seed,
        output,
        verbose,
    } = args;
    check_input_paths("input", &inputs)?;
    kdist_calibrate::run(
        &inputs,
        &kdist_calibrate::Params {
            k,
            screen_backend,
            minimizer_k,
            minimizer_w,
            cutoff,
            leak_pct,
            band,
            max_pairs,
            max_uniques,
            per_sample,
            nearest_parent,
            from_dada,
            from_dada_pooled,
            derive_cutoff,
            derive_uniform_pairs,
            derep_dir,
            threads,
            seed,
            output,
            verbose,
        },
    )?;
    Ok(())
}

/// Run the self-consistency loop on loaded samples and write the error-model
/// JSON: everything `learn-errors` and `errors-from-sample` do once their
/// inputs are loaded.
#[allow(clippy::too_many_arguments)]
fn learn_and_write(
    tag: &'static str,
    all_inputs: Vec<Vec<dada::RawInput>>,
    (err_fun, align_params, dada_params): (ErrFun, AlignParams, dada::DadaParams),
    max_consist: usize,
    diag: cli::LearnDiagArgs,
    pool: &rayon::ThreadPool,
    verbose: bool,
    compact: bool,
    output: Option<PathBuf>,
) -> io::Result<()> {
    if let Some(ref dir) = diag.diag_dir {
        std::fs::create_dir_all(dir)?;
    }

    let params_snapshot =
        build_learned_err_params(&err_fun, max_consist, &dada_params, &align_params);

    if let Some(ref dir) = diag.cluster_trace_dir {
        std::fs::create_dir_all(dir)?;
    }
    let trace_params = cluster_trace::TraceParams {
        no_members: diag.trace_no_members,
        min_abund: diag.trace_min_abund,
    };

    let result = pool.install(|| {
        learn_errors(
            all_inputs,
            &err_fun,
            dada_params,
            max_consist,
            LearnDiagOptions {
                verbose,
                diag_dir: diag.diag_dir.as_deref(),
                cluster_trace_dir: diag.cluster_trace_dir.as_deref(),
                trace_params,
            },
        )
    })?;

    // Serialize: represent the three matrices as Vec<Vec<T>> (16 rows × nq cols).
    #[derive(Serialize)]
    struct LearnErrorsOutput {
        nq: usize,
        converged: bool,
        stop_reason: learn_errors::StopReason,
        iterations: usize,
        /// Provenance: parameters used for the dada_uniques runs that
        /// produced this err model. Embedded so a downstream `dada`
        /// invocation can validate or inherit them. See
        /// `LearnedErrParams` for field details.
        params: LearnedErrParams,
        /// Transition counts: 16 rows (ref_nt*4+query_nt), nq columns.
        trans: Vec<Vec<u32>>,
        /// Error rates fed into the final DADA run: 16 × nq.
        err_in: Vec<Vec<f64>>,
        /// Error rates estimated from `trans`: 16 × nq.
        err_out: Vec<Vec<f64>>,
    }

    fn flat_to_rows_u32(flat: &[u32], nq: usize) -> Vec<Vec<u32>> {
        (0..16)
            .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
            .collect()
    }
    fn flat_to_rows_f64(flat: &[f64], nq: usize) -> Vec<Vec<f64>> {
        (0..16)
            .map(|r| flat[r * nq..(r + 1) * nq].to_vec())
            .collect()
    }

    let out = LearnErrorsOutput {
        nq: result.nq,
        converged: result.converged,
        stop_reason: result.stop_reason,
        iterations: result.iterations,
        params: params_snapshot,
        trans: flat_to_rows_u32(&result.trans, result.nq),
        err_in: flat_to_rows_f64(&result.err_in, result.nq),
        err_out: flat_to_rows_f64(&result.err_out, result.nq),
    };

    let tagged = Tagged::new(tag, out);
    let json = if compact {
        serde_json::to_string(&tagged)
    } else {
        serde_json::to_string_pretty(&tagged)
    }
    .map_err(io::Error::other)?;

    match output {
        Some(path) => misc::write_maybe_gz(&path, json.as_bytes())?,
        None => println!("{json}"),
    }
    Ok(())
}

/// Resolve a [`LoessConfig`] from CLI inputs: preset + per-knob overrides.
/// `--loess-surface`, `--loess-cell`, `--loess-max-rate`, and `--loess-min-rate`
/// each override the preset's value for that knob if supplied.  `--loess-cell`
/// is ignored unless the resolved surface is `Interpolate`.
///
/// `--loess-preset` is deprecated (#205). It survives as an alias because it
/// appears in shipped docs, the concordance runners' `ERRFUN_ARGS`, and users'
/// scripts -- silently ignoring it would look exactly like the preset having
/// no effect, which is a failure mode this project has already paid for. It
/// warns on use, and `--loess-surface` wins if both are given.
fn resolve_loess_config(
    preset: Option<&str>,
    surface: Option<&str>,
    cell: Option<f64>,
    max_rate: Option<f64>,
    min_rate: Option<f64>,
) -> LoessConfig {
    let base = match preset {
        Some(p) => {
            let mapped = if p == "r-dada2" {
                "interpolate"
            } else {
                "direct"
            };
            eprintln!(
                "[learn-errors] WARNING: --loess-preset is deprecated; use \
                 --loess-surface {mapped}. `interpolate` is now the default \
                 (it is the surface R's loess() uses), so \
                 `--loess-preset r-dada2` is redundant."
            );
            match p {
                "r-dada2" => LoessConfig::r_dada2(),
                // The old `default` preset selected Direct; preserve that
                // meaning for anyone who passed it explicitly.
                _ => LoessConfig {
                    surface: LoessSurface::Direct,
                    ..LoessConfig::default()
                },
            }
        }
        None => LoessConfig::default(),
    };
    let surface = match surface {
        Some("interpolate") => {
            let c = cell.unwrap_or(match base.surface {
                LoessSurface::Interpolate { cell } => cell,
                LoessSurface::Direct => 0.2,
            });
            LoessSurface::Interpolate { cell: c }
        }
        Some("direct") => LoessSurface::Direct,
        _ => match base.surface {
            LoessSurface::Interpolate { cell: base_cell } => LoessSurface::Interpolate {
                cell: cell.unwrap_or(base_cell),
            },
            LoessSurface::Direct => LoessSurface::Direct,
        },
    };
    LoessConfig {
        surface,
        max_error_rate: max_rate.unwrap_or(base.max_error_rate),
        min_error_rate: min_rate.unwrap_or(base.min_error_rate),
    }
}

/// Build a [`LearnedErrParams`] snapshot from the resolved errfun + dada/align
/// params, for embedding in the err-model JSON. Captures everything dada cares
/// about so a downstream invocation can validate or inherit.
fn build_learned_err_params(
    errfun: &ErrFun,
    max_consist: usize,
    dp: &dada::DadaParams,
    ap: &AlignParams,
) -> LearnedErrParams {
    let (errfun_name, errfun_pseudocount, errfun_bins, errfun_cmd, loess_cfg) = match errfun {
        ErrFun::Loess { config } => ("loess", None, None, None, Some(config)),
        ErrFun::Noqual {
            pseudocount,
            config,
        } => ("noqual", Some(*pseudocount), None, None, Some(config)),
        ErrFun::BinnedQual { bins, config } => {
            ("binned-qual", None, Some(bins.clone()), None, Some(config))
        }
        ErrFun::PacBio { config } => ("pacbio", None, None, None, Some(config)),
        ErrFun::External { command } => ("external", None, None, Some(command.clone()), None),
    };
    LearnedErrParams {
        errfun: errfun_name.to_string(),
        errfun_pseudocount,
        errfun_bins,
        errfun_cmd,
        loess: loess_cfg.map(Into::into),
        max_consist,
        omega_a: dp.omega_a,
        // Deliberately not embedded: learn-time `omega_c` (R default 0)
        // differs from dada-time (R default 1e-40) and must not transfer.
        omega_c: None,
        omega_p: dp.omega_p,
        min_fold: dp.min_fold,
        min_hamming: dp.min_hamming,
        min_abund: dp.min_abund,
        detect_singletons: dp.detect_singletons,
        use_quals: dp.use_quals,
        greedy: dp.greedy,
        max_clust: dp.max_clust,
        match_score: ap.match_score,
        mismatch: ap.mismatch,
        gap_p: ap.gap_p,
        homo_gap_p: ap.homo_gap_p,
        use_kmers: ap.use_kmers,
        kdist_cutoff: ap.kdist_cutoff,
        kmer_size: ap.kmer_size,
        band: ap.band,
        vectorized: ap.vectorized,
        gapless: ap.gapless,
        backend: ap.backend,
        wfa_max_edits: ap.wfa_max_edits,
        screen_backend: ap.screen_backend,
        minimizer_k: ap.minimizer_k,
        minimizer_w: ap.minimizer_w,
    }
}

/// Build the errfun and DADA parameters `learn-errors` and `errors-from-sample`
/// share. `dada`-side resolution, which can inherit from an error model, is
/// `commands::dada::resolve_dada_params`.
fn resolve_learn_params(
    fit: &cli::ErrModelFitArgs,
    denoise: cli::LearnDenoiseArgs,
    experimental: cli::ExperimentalArgs,
    threads: usize,
    verbose: bool,
) -> io::Result<(ErrFun, AlignParams, dada::DadaParams)> {
    let cli::LearnDenoiseArgs {
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
    // R's HOMOPOLYMER_GAP_PENALTY = NULL tracks GAP_PENALTY. R also
    // normalizes a positive penalty to negative before comparing them
    // (dada.R:223-227), so `--homo-gap-p 1` means the same as `-1`.
    let gap_p = gap_p.unwrap_or(-8);
    let gap_p = if gap_p > 0 { -gap_p } else { gap_p };
    let homo_gap_p = homo_gap_p.unwrap_or(gap_p);
    let homo_gap_p = if homo_gap_p > 0 {
        -homo_gap_p
    } else {
        homo_gap_p
    };
    note_homopolymer_gapping(verbose, gap_p, homo_gap_p);
    let greedy = greedy.unwrap_or(true);
    let use_quals = use_quals.unwrap_or(true);
    let loess_config = resolve_loess_config(
        fit.loess_preset.as_deref(),
        fit.loess_surface.as_deref(),
        fit.loess_cell,
        fit.loess_max_rate,
        fit.loess_min_rate,
    );
    let err_fun = match fit.errfun.as_str() {
        "loess" => ErrFun::Loess {
            config: loess_config,
        },
        "noqual" => ErrFun::Noqual {
            pseudocount: fit.pseudocount,
            config: loess_config,
        },
        "binned-qual" => {
            let bins = fit.binned_quals.clone().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--binned-quals is required when --errfun binned-qual is used",
                )
            })?;
            ErrFun::BinnedQual {
                bins,
                config: loess_config,
            }
        }
        "pacbio" => ErrFun::PacBio {
            config: loess_config,
        },
        "external" => {
            let command = fit.errfun_cmd.clone().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--errfun-cmd is required when --errfun external is used",
                )
            })?;
            if command.trim().is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--errfun-cmd cannot be empty",
                ));
            }
            ErrFun::External { command }
        }
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Unknown errfun '{other}'; expected one of: loess, noqual, binned-qual, pacbio, external"
                ),
            ));
        }
    };

    let align_params = AlignParams {
        backend: align_backend.unwrap_or_default(),
        wfa_max_edits: wfa_max_edits.unwrap_or(WFA_MAX_EDITS_DEFAULT),
        match_score,
        mismatch,
        gap_p,
        homo_gap_p,
        use_kmers: !no_kmer_screen,
        kdist_cutoff,
        kmer_size,
        screen_backend: screen_backend.unwrap_or_default(),
        minimizer_k: minimizer_k.unwrap_or(minimizers::MINIMIZER_K),
        minimizer_w: minimizer_w.unwrap_or(minimizers::MINIMIZER_W),
        screen_audit,
        band,
        vectorized: true,
        gapless: true,
    };

    let dada_params = dada::DadaParams::new(
        align_params,
        Vec::new(),
        0,
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
        // learn-errors has no --metrics-json yet; keep
        // --verbose's measurements exactly as they were.
        if verbose {
            MeasureLevel::Attribution
        } else {
            MeasureLevel::Off
        },
    );
    Ok((err_fun, align_params, dada_params))
}
