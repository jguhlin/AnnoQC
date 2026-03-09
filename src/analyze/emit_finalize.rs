use crate::config_types::FileConfig;
use crate::taxonomy::TaxonomyEvidence;
use crate::*;
use std::collections::HashMap;

pub(in crate::analyze) struct TimingInputs {
    pub(in crate::analyze) diamond_secs: f64,
    pub(in crate::analyze) ecs_secs: f64,
    pub(in crate::analyze) intrinsic_secs: f64,
    pub(in crate::analyze) consensus_secs: f64,
    pub(in crate::analyze) emit_secs: f64,
}

pub(in crate::analyze) struct FinalizeInputs<'a> {
    pub(in crate::analyze) cfg: &'a EffectiveConfig,
    pub(in crate::analyze) file_cfg: &'a FileConfig,
    pub(in crate::analyze) args: &'a AnalyzeArgs,
    pub(in crate::analyze) report_format: ReportFormat,
    pub(in crate::analyze) calibration: CalibrationSettings,
    pub(in crate::analyze) metrics: &'a [GeneMetrics],
    pub(in crate::analyze) render_summary: &'a ecs::RenderSummary,
    pub(in crate::analyze) taxonomy_enabled_effective: bool,
    pub(in crate::analyze) hmmer_requested: bool,
    pub(in crate::analyze) mafft_requested: bool,
    pub(in crate::analyze) plugin_count: usize,
    pub(in crate::analyze) rule_count: usize,
    pub(in crate::analyze) taxsum_map: &'a HashMap<String, Option<TaxonomyEvidence>>,
    pub(in crate::analyze) resolver: Option<&'a taxonomy::TaxonomyResolver>,
    pub(in crate::analyze) hmmsum_map: &'a HashMap<String, HmmscanSummary>,
    pub(in crate::analyze) alignment_map: &'a HashMap<String, AlignmentMetrics>,
    pub(in crate::analyze) alignment_secs: &'a HashMap<String, f64>,
    pub(in crate::analyze) hmmer_secs: &'a HashMap<String, f64>,
    pub(in crate::analyze) backfill_used: u64,
    pub(in crate::analyze) backfill_added_total: u64,
    pub(in crate::analyze) refprot_used_panels: u64,
    pub(in crate::analyze) timings: TimingInputs,
}

pub(in crate::analyze) fn finalize_outputs(
    input: FinalizeInputs<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let hmmer_genes = input.hmmsum_map.len() as u64;
    let hmmer_hits_total: u64 = input.hmmsum_map.values().map(|s| s.hits_count as u64).sum();
    let mafft_alignments = input.alignment_map.len() as u64;
    let mut totals = build_base_totals(
        input.metrics.len(),
        hmmer_genes,
        hmmer_hits_total,
        mafft_alignments,
        &format!("{:?}", input.args.diamond_mode),
        input.backfill_used,
        input.backfill_added_total,
        input.refprot_used_panels,
    );
    let steps = build_steps(
        input.timings.diamond_secs,
        input.timings.ecs_secs,
        input.timings.intrinsic_secs,
        input.timings.consensus_secs,
        input.timings.emit_secs,
    );
    if input.taxonomy_enabled_effective {
        augment_totals_from_taxonomy(
            &input.cfg.out,
            &mut totals,
            input.metrics,
            input.taxsum_map,
            input.resolver,
            input
                .file_cfg
                .refprot
                .as_ref()
                .map(|cfg_rp| RefprotSelectionConfig {
                    enabled: cfg_rp.enabled.unwrap_or(false),
                    readme_path: cfg_rp.readme_path.as_deref(),
                    max_scopes: cfg_rp.max_scopes.unwrap_or(2),
                    max_proteomes: cfg_rp.max_proteomes.unwrap_or(10),
                }),
        );
    }
    write_slowest_genes(&input.cfg.out, input.alignment_secs, input.hmmer_secs)?;
    let run_summary = build_run_summary_obj(RunSummaryInputs {
        schema_version: OUTPUT_SCHEMA_VERSION,
        total_genes: input.metrics.len(),
        render_summary: input.render_summary,
        aligner: format!("{:?}", input.args.aligner).to_lowercase(),
        report_format: format!("{:?}", input.report_format).to_lowercase(),
        resume: input.args.resume,
        taxonomy_enabled: input.taxonomy_enabled_effective,
        hmmer_enabled: input.hmmer_requested,
        alignment_enabled: input.mafft_requested,
        genomic_enabled: input.args.gff.is_some(),
        plugins: input.plugin_count,
        rules: input.rule_count,
        nucleotide: input.args.nucleotide,
        calibration: format!("{:?}", input.calibration.mode).to_lowercase(),
        timings: steps.clone(),
    });
    write_run_summary(&input.cfg.out, &run_summary)?;
    let runm = build_run_metrics_obj(OUTPUT_SCHEMA_VERSION, steps, totals);
    write_run_metrics(&input.cfg.out, &runm)?;
    Ok(())
}
