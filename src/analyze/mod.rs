use std::sync::Arc;

use crate::preflight;
use crate::scoring_support::print_scoring_rubric;
use crate::*;

mod bootstrap;
mod manifest;
mod pipeline;

use bootstrap::AnalyzeBootstrap;
use pipeline::{
    build_consensus_stage, finalize_outputs, prepare_emit_stage, run_diamond_stage,
    run_heavy_stage, run_scoring_and_render_stage, FinalizeInputs, ScoringRenderInputs,
    TimingInputs,
};

pub(crate) fn run_analyze(
    cli_config_path: Option<&str>,
    args: AnalyzeArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    let AnalyzeBootstrap {
        file_cfg,
        mut cfg,
        calibration,
        report_format,
        rhai_paths,
    } = bootstrap::bootstrap_analyze(cli_config_path, &args)?;

    if args.dry_run {
        print_scoring_rubric(&file_cfg.scoring, calibration, args.classify_no_data);
        return Ok(());
    }

    std::fs::create_dir_all(&cfg.out)?;

    let tools = preflight(
        &cfg.diamond_bin,
        args.mafft_bin.as_deref(),
        args.hmmscan_bin.as_deref(),
    );
    if tools.diamond.is_none() {
        return Err("DIAMOND not found or --version check failed".into());
    }
    if let Err(e) = preflight::check_diamond_db(&cfg.diamond_bin, &cfg.db) {
        return Err(e.into());
    }

    if args.nucleotide {
        cfg.fasta = orf::translate_nt_fasta_to_protein(&cfg.fasta)?;
    }

    let log_json = matches!(args.log_format, LogFormat::Json);
    let diamond = run_diamond_stage(&cfg, &file_cfg, &args, log_json)?;

    let t_ecs = step_start("ecs", log_json);
    let mut metrics: Vec<GeneMetrics> = run_scheduler(EcsConfig {
        fasta_path: cfg.fasta.clone(),
        diamond_tsv: diamond.diamond_tsv.to_string_lossy().to_string(),
        threads: cfg.threads,
        log_json,
    });
    let ecs_secs = step_finish("ecs", t_ecs, log_json);

    if args.resume {
        let rendered_ids = manifest::load_rendered_gene_ids(&cfg.out);
        if !rendered_ids.is_empty() {
            let before = metrics.len();
            metrics.retain(|m| !rendered_ids.contains(&m.gene_id));
            let skipped = before.saturating_sub(metrics.len());
            if skipped > 0 {
                log::info!("resume: skipping {} already-rendered genes", skipped);
            }
        }
    }
    if metrics.is_empty() {
        log::info!("resume: no pending genes to process");
        return Ok(());
    }

    let t_intrinsic = step_start("intrinsic", log_json);
    let intrinsic_map = Arc::new(compute_intrinsic_for_ids(&cfg.fasta, &metrics)?);
    let intrinsic_secs = step_finish("intrinsic", t_intrinsic, log_json);

    let consensus = build_consensus_stage(&cfg, &file_cfg, &args, &metrics, &diamond, log_json)?;
    let grouped = consensus.grouped;
    let cons_cfg = consensus.cons_cfg;
    let panel_map = consensus.panel_map;
    let len_map = consensus.len_map;
    let panel_prov_map = consensus.panel_prov_map;
    let panel_prov_rows = consensus.panel_prov_rows;
    let refprot_used_panels = consensus.refprot_used_panels;
    let taxonomy_hits_map = consensus.taxonomy_hits_map;
    let consensus_secs = consensus.consensus_secs;
    let backfill_used = consensus.backfill_used;
    let backfill_added_total = consensus.backfill_added_total;

    let config_sv = StructVarOverrideInputs {
        min_hsp_len: file_cfg.structvar.as_ref().and_then(|c| c.min_hsp_len),
        min_hsp_frac: file_cfg.structvar.as_ref().and_then(|c| c.min_hsp_frac),
        fusion_min_gap: file_cfg.structvar.as_ref().and_then(|c| c.fusion_min_gap),
        dup_max_gap: file_cfg.structvar.as_ref().and_then(|c| c.dup_max_gap),
        split_delta: file_cfg.structvar.as_ref().and_then(|c| c.split_delta),
        min_subject_cov: file_cfg.structvar.as_ref().and_then(|c| c.min_subject_cov),
        orient_majority: file_cfg.structvar.as_ref().and_then(|c| c.orient_majority),
    };
    let cli_sv = StructVarOverrideInputs {
        min_hsp_len: args.sv_min_hsp_len,
        min_hsp_frac: args.sv_min_hsp_frac,
        fusion_min_gap: args.sv_fusion_min_gap,
        dup_max_gap: args.sv_dup_max_gap,
        split_delta: args.sv_split_delta,
        min_subject_cov: args.sv_min_subject_cov,
        orient_majority: args.sv_orient_majority,
    };
    let sv_th = resolve_thresholds(&config_sv, &cli_sv);
    let structvar_map = analyze_structvar(&metrics, &grouped, &sv_th);

    dump_matches_fasta(
        &cfg.out,
        cfg.reference_fasta.as_deref(),
        args.dump_matches_gene.as_deref(),
        args.dump_matches_best,
        &metrics,
        &panel_map,
        intrinsic_map.as_ref(),
    )?;
    write_panel_sources_csv(&cfg.out, &panel_prov_rows)?;

    let heavy = run_heavy_stage(
        &cfg,
        &file_cfg,
        &args,
        &metrics,
        &panel_map,
        intrinsic_map.as_ref(),
        cons_cfg.min_hits,
        log_json,
    );
    let mafft_requested = heavy.mafft_requested;
    let hmmer_requested = heavy.hmmer_requested;
    let orphan_analysis_enabled = heavy.orphan_analysis_enabled;
    let alignment_map = Arc::new(heavy.alignment_map.clone());
    let hmmsum_map = Arc::new(heavy.hmmsum_map.clone());
    let domains_arch_map = Arc::new(heavy.domains_arch_map.clone());
    let len_map = Arc::new(len_map);
    let orphan_map = Arc::new(heavy.orphan_map.clone());
    let structvar_map = Arc::new(structvar_map);
    let taxonomy_hits_map = Arc::new(taxonomy_hits_map);

    let emit = prepare_emit_stage(
        &cfg,
        &file_cfg,
        &args,
        &metrics,
        &diamond,
        &tools,
        calibration,
        report_format,
        &rhai_paths,
        Arc::clone(&taxonomy_hits_map),
    )?;
    let scoring = run_scoring_and_render_stage(ScoringRenderInputs {
        cfg: &cfg,
        file_cfg: &file_cfg,
        args: &args,
        metrics: &metrics,
        report_format,
        tools: &tools,
        checksums: emit.checksums,
        stats: emit.stats,
        taxonomy_setup: emit.taxonomy_setup,
        intrinsic_map: Arc::clone(&intrinsic_map),
        alignment_map: Arc::clone(&alignment_map),
        hmmsum_map: Arc::clone(&hmmsum_map),
        domains_arch_map: Arc::clone(&domains_arch_map),
        len_map: Arc::clone(&len_map),
        orphan_map: Arc::clone(&orphan_map),
        structvar_map: Arc::clone(&structvar_map),
        panel_prov_map,
        domains_arch_dbg: &heavy.domains_arch_dbg,
        rhai_paths: &rhai_paths,
        taxonomy_enabled_effective: emit.taxonomy_enabled_effective,
        hmmer_requested,
        mafft_requested,
        orphan_analysis_enabled,
        calibration,
        log_json,
    })?;
    finalize_outputs(FinalizeInputs {
        cfg: &cfg,
        file_cfg: &file_cfg,
        args: &args,
        report_format,
        calibration,
        metrics: &metrics,
        render_summary: &scoring.render_summary,
        taxonomy_enabled_effective: scoring.taxonomy_enabled_effective,
        hmmer_requested,
        mafft_requested,
        plugin_count: scoring.plugin_count,
        rule_count: scoring.rule_count,
        taxsum_map: scoring.taxsum_map.as_ref(),
        resolver: scoring.resolver.as_deref(),
        hmmsum_map: &heavy.hmmsum_map,
        alignment_map: &heavy.alignment_map,
        alignment_secs: &heavy.alignment_secs,
        hmmer_secs: &heavy.hmmer_secs,
        backfill_used,
        backfill_added_total,
        refprot_used_panels,
        timings: TimingInputs {
            diamond_secs: diamond.diamond_secs,
            ecs_secs,
            intrinsic_secs,
            consensus_secs,
            emit_secs: scoring.emit_secs,
        },
    })?;
    Ok(())
}
