use crate::config_types::FileConfig;
use crate::*;
use std::collections::HashMap;
use std::sync::Arc;

use crate::analyze::manifest::{
    build_config_snapshot, collect_plugin_manifest, collect_rule_manifest, write_run_manifest,
};
use crate::analyze::pipeline::DiamondStageResult;

pub(in crate::analyze) struct EmitStageResult {
    pub(in crate::analyze) checksums: Checksums,
    pub(in crate::analyze) stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    pub(in crate::analyze) taxonomy_setup: taxonomy_support::TaxonomySetup,
    pub(in crate::analyze) taxonomy_enabled_effective: bool,
}

pub(in crate::analyze) fn prepare_emit_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    metrics: &[GeneMetrics],
    diamond: &DiamondStageResult,
    tools: &preflight::ToolVersions,
    calibration: CalibrationSettings,
    report_format: ReportFormat,
    rhai_paths: &[String],
    taxonomy_hits_map: Arc<HashMap<String, Vec<diamond::DiamondHitRow>>>,
) -> Result<EmitStageResult, Box<dyn std::error::Error>> {
    let qlen_map: HashMap<String, usize> = metrics
        .iter()
        .map(|m| (m.gene_id.clone(), m.length))
        .collect();
    let stats =
        Arc::new(parse_tsv_stats(&diamond.diamond_tsv, Some(&qlen_map)).unwrap_or_default());
    let checksums = collect_checksums(cfg)?;
    let snapshot = build_config_snapshot(cfg, args, file_cfg, &calibration, report_format);
    let plugin_manifest = collect_plugin_manifest(&args.plugin);
    let rule_manifest = collect_rule_manifest(rhai_paths);
    write_run_manifest(
        cfg,
        tools,
        &checksums,
        &snapshot,
        &plugin_manifest,
        &rule_manifest,
    )?;

    let mut taxonomy_enabled_effective = args.enable_taxonomy
        || file_cfg
            .taxonomy
            .as_ref()
            .and_then(|t| t.enabled)
            .unwrap_or(false);
    let cache_path = args
        .taxonomy_cache
        .as_deref()
        .or_else(|| {
            file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.cache_path.as_deref())
        })
        .or(file_cfg.taxonomy_cache.as_deref());
    let taxdump_dir = args
        .taxonomy_taxdump_dir
        .as_deref()
        .or_else(|| {
            file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.taxdump_dir.as_deref())
        })
        .or(file_cfg.taxonomy_taxdump_dir.as_deref());
    let tax_cfg = TaxonomyConsensusConfig {
        min_hits: args
            .taxonomy_min_consensus
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_consensus))
            .unwrap_or(5)
            .max(1),
        top_hits: args
            .taxonomy_top_hits
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.top_hits))
            .unwrap_or(20)
            .max(1),
        min_support: args
            .taxonomy_min_support
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_support))
            .unwrap_or(0.75),
        coarse_rank_index: args
            .taxonomy_coarse_rank_index
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.coarse_rank_index))
            .unwrap_or(1),
        coarse_min_support: args
            .taxonomy_coarse_min_support
            .or_else(|| {
                file_cfg
                    .taxonomy
                    .as_ref()
                    .and_then(|t| t.coarse_min_support)
            })
            .unwrap_or(0.6),
    };
    let taxonomy_setup = build_taxonomy_setup(
        metrics,
        &taxonomy_hits_map,
        stats.as_ref(),
        taxonomy_enabled_effective,
        cache_path,
        cfg.reference_fasta.as_deref(),
        taxdump_dir,
        tax_cfg,
        TaxonomyWarnings {
            expected_domain: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.expected_domain.clone()),
            warn_non_target_min_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_frac)
                .unwrap_or(0.05),
            warn_non_target_min_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_hits)
                .unwrap_or(200),
            warn_non_target_strong_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_frac)
                .unwrap_or(0.10),
            warn_non_target_strong_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_hits)
                .unwrap_or(500),
            warn_genus_min_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_frac)
                .unwrap_or(0.15),
            warn_genus_min_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_hits)
                .unwrap_or(300),
            low_coverage_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.low_coverage_frac)
                .unwrap_or(0.05),
        },
    )?;
    taxonomy_enabled_effective = taxonomy_setup.enabled;

    Ok(EmitStageResult {
        checksums,
        stats,
        taxonomy_setup,
        taxonomy_enabled_effective,
    })
}
