use crate::config_types::FileConfig;
use crate::scoring_support::{
    build_component_scores, build_scores_map, export_high_sequences, scoring_thresholds,
};
use crate::taxonomy::TaxonomyEvidence;
use crate::*;
use std::collections::HashMap;
use std::sync::Arc;

pub(in crate::analyze) struct ScoringRenderInputs<'a> {
    pub(in crate::analyze) cfg: &'a EffectiveConfig,
    pub(in crate::analyze) file_cfg: &'a FileConfig,
    pub(in crate::analyze) args: &'a AnalyzeArgs,
    pub(in crate::analyze) metrics: &'a [GeneMetrics],
    pub(in crate::analyze) report_format: ReportFormat,
    pub(in crate::analyze) tools: &'a preflight::ToolVersions,
    pub(in crate::analyze) checksums: Checksums,
    pub(in crate::analyze) stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    pub(in crate::analyze) taxonomy_setup: taxonomy_support::TaxonomySetup,
    pub(in crate::analyze) intrinsic_map:
        Arc<HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>>,
    pub(in crate::analyze) alignment_map: Arc<HashMap<String, AlignmentMetrics>>,
    pub(in crate::analyze) hmmsum_map: Arc<HashMap<String, HmmscanSummary>>,
    pub(in crate::analyze) domains_arch_map: Arc<HashMap<String, f64>>,
    pub(in crate::analyze) len_map: Arc<HashMap<String, LengthSummary>>,
    pub(in crate::analyze) orphan_map: Arc<HashMap<String, hmmer::OrphanAnalysis>>,
    pub(in crate::analyze) structvar_map: Arc<HashMap<String, structvar::StructVar>>,
    pub(in crate::analyze) panel_prov_map: HashMap<String, PanelProvenanceCounts>,
    pub(in crate::analyze) domains_arch_dbg: &'a [DomainsArchDebugRow],
    pub(in crate::analyze) rhai_paths: &'a [String],
    pub(in crate::analyze) taxonomy_enabled_effective: bool,
    pub(in crate::analyze) hmmer_requested: bool,
    pub(in crate::analyze) mafft_requested: bool,
    pub(in crate::analyze) orphan_analysis_enabled: bool,
    pub(in crate::analyze) calibration: CalibrationSettings,
    pub(in crate::analyze) log_json: bool,
}

pub(in crate::analyze) struct ScoringRenderResult {
    pub(in crate::analyze) render_summary: ecs::RenderSummary,
    pub(in crate::analyze) emit_secs: f64,
    pub(in crate::analyze) taxonomy_enabled_effective: bool,
    pub(in crate::analyze) plugin_count: usize,
    pub(in crate::analyze) rule_count: usize,
    pub(in crate::analyze) taxsum_map: Arc<HashMap<String, Option<TaxonomyEvidence>>>,
    pub(in crate::analyze) resolver: Option<Arc<taxonomy::TaxonomyResolver>>,
}

pub(in crate::analyze) fn run_scoring_and_render_stage(
    input: ScoringRenderInputs<'_>,
) -> Result<ScoringRenderResult, Box<dyn std::error::Error>> {
    let t_emit = step_start("emit_outputs", input.log_json);
    let resolver = input.taxonomy_setup.resolver.clone();
    let taxsum_map = Arc::clone(&input.taxonomy_setup.taxsum_map);
    let genomic_map = load_genomic_context(
        input.args.gff.as_deref(),
        input.args.genome.as_deref(),
        input.log_json,
    );
    let rnaseq_requested = input.args.rnaseq_file.is_some()
        || input
            .file_cfg
            .rnaseq
            .as_ref()
            .and_then(|r| r.enabled)
            .unwrap_or(false);
    let rnaseq_min_tpm = if input.args.rnaseq_file.is_some() {
        input.args.rnaseq_min_tpm
    } else {
        input
            .file_cfg
            .rnaseq
            .as_ref()
            .and_then(|r| r.min_tpm)
            .unwrap_or(input.args.rnaseq_min_tpm)
    };
    let (rnaseq_enabled, rnaseq_map) = load_rnaseq_data(
        rnaseq_requested,
        input
            .args
            .rnaseq_file
            .clone()
            .or_else(|| input.file_cfg.rnaseq.as_ref().and_then(|r| r.file.clone())),
        rnaseq_min_tpm,
    );
    let t_scoring = step_start("scoring", input.log_json);
    let (scores_map_inner, raw_scores_inner) = build_scores_map(
        input.metrics,
        input.stats.as_ref(),
        input.intrinsic_map.as_ref(),
        &input.file_cfg.scoring,
        input.taxonomy_enabled_effective,
        input.hmmer_requested,
        if input.hmmer_requested {
            Some(input.hmmsum_map.as_ref())
        } else {
            None
        },
        Some(input.domains_arch_map.as_ref()),
        Some(input.len_map.as_ref()),
        if input.orphan_analysis_enabled {
            Some(input.orphan_map.as_ref())
        } else {
            None
        },
        if input.taxonomy_enabled_effective {
            Some(taxsum_map.as_ref())
        } else {
            None
        },
        Some(input.alignment_map.as_ref()),
        Some(input.structvar_map.as_ref()),
        genomic_map.as_ref().map(|gm| gm.as_ref()),
        Some(rnaseq_map.as_ref()),
        rnaseq_enabled,
        input.calibration,
        input.args.classify_no_data,
    );
    let comp_map = build_component_scores(
        input.metrics,
        input.stats.as_ref(),
        input.intrinsic_map.as_ref(),
        input.taxonomy_enabled_effective,
        if input.orphan_analysis_enabled {
            Some(input.orphan_map.as_ref())
        } else {
            None
        },
        if input.args.enable_taxonomy {
            Some(taxsum_map.as_ref())
        } else {
            None
        },
        Some(input.alignment_map.as_ref()),
        Some(input.len_map.as_ref()),
        genomic_map.as_ref().map(|gm| gm.as_ref()),
        Some(rnaseq_map.as_ref()),
        rnaseq_enabled,
        &input.file_cfg.scoring,
    );
    let _ = step_finish("scoring", t_scoring, input.log_json);
    if input.args.export_high || input.file_cfg.export_high.unwrap_or(false) {
        let export_high_path = input
            .args
            .export_high_path
            .clone()
            .or_else(|| input.file_cfg.export_high_path.clone());
        if let Some((path, count)) = export_high_sequences(
            &input.cfg.out,
            export_high_path.as_deref(),
            input.metrics,
            input.intrinsic_map.as_ref(),
            &scores_map_inner,
        )? {
            log::info!(
                "exported {} high-scoring genes to {}",
                count,
                path.display()
            );
        } else {
            log::info!("no high-scoring genes to export");
        }
    }
    let scores_map = Arc::new(scores_map_inner);
    let raw_scores_map = Arc::new(raw_scores_inner);
    let comp_map = Arc::new(comp_map);
    let panel_prov_map = Arc::new(input.panel_prov_map);
    let extensions = load_extensions(&input.args.plugin, input.rhai_paths);
    let plugin_count = extensions.plugin_count;
    let rule_count = extensions.rule_count;
    let (th_high, th_med) = scoring_thresholds(&input.file_cfg.scoring);
    let render_ctx = build_render_context(RenderContextInputs {
        stats: Arc::clone(&input.stats),
        intrinsic_map: Arc::clone(&input.intrinsic_map),
        alignment_map: Arc::clone(&input.alignment_map),
        hmmsum_map: Arc::clone(&input.hmmsum_map),
        taxsum_map: Arc::clone(&taxsum_map),
        taxonomy_resolver: resolver.clone(),
        genomic_map: genomic_map.clone(),
        scores_map: Arc::clone(&scores_map),
        raw_scores_map: Arc::clone(&raw_scores_map),
        comp_map: Arc::clone(&comp_map),
        arch_map: Arc::clone(&input.domains_arch_map),
        len_map: Arc::clone(&input.len_map),
        orphan_map: Arc::clone(&input.orphan_map),
        structvar_map: Arc::clone(&input.structvar_map),
        panel_prov_map,
        cov_delta_thresh: input.args.coverage_delta_threshold,
        taxonomy_enabled: input.taxonomy_enabled_effective,
        taxonomy_warnings: input.taxonomy_setup.warnings.clone(),
        orphan_analysis_enabled: input.orphan_analysis_enabled,
        csv_verbose: input.args.csv_verbose,
        mafft_missing_exon_thresh: input
            .args
            .alignment_missing_exon
            .or(input.file_cfg.alignment_missing_exon)
            .unwrap_or(30),
        mafft_retained_intron_thresh: input
            .args
            .alignment_retained_intron
            .or(input.file_cfg.alignment_retained_intron)
            .unwrap_or(30),
        features_string: build_features_string(&FeatureSummary {
            profile_name: &format!("{:?}", input.args.profile),
            has_genomic_input: input.args.gff.is_some() || input.args.genome.is_some(),
            nucleotide_input: input.args.nucleotide,
            taxonomy_enabled: input.taxonomy_enabled_effective,
            genomic_enabled: genomic_map.is_some(),
            plugins_enabled: plugin_count > 0,
            rules_enabled: rule_count > 0,
            hmmer_enabled: input.hmmer_requested,
            orphan_enabled: input.orphan_analysis_enabled,
            alignment_enabled: input.mafft_requested,
        }),
        extensions,
        rnaseq_map: Arc::clone(&rnaseq_map),
        rnaseq_enabled,
        th_high,
        th_med,
    });
    let render_summary = run_render_pipeline(
        &input.cfg.out,
        input.metrics,
        render_ctx,
        resolve_render_max_jobs(
            input.args.render_max_jobs,
            input.file_cfg.render_max_jobs,
            input.cfg.threads,
        ),
        input.tools,
        &input.checksums,
        input.report_format.output_config(input.args.resume),
    )?;
    write_domains_arch_debug(&input.cfg.out, input.domains_arch_dbg)?;
    write_structvar_summary(&input.cfg.out, input.structvar_map.as_ref())?;
    let emit_secs = step_finish("emit_outputs", t_emit, input.log_json);

    Ok(ScoringRenderResult {
        render_summary,
        emit_secs,
        taxonomy_enabled_effective: input.taxonomy_enabled_effective,
        plugin_count,
        rule_count,
        taxsum_map,
        resolver,
    })
}
