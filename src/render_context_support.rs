use std::collections::HashMap;
use std::sync::Arc;

use crate::analysis_support::LoadedExtensions;
use crate::consensus_support::{LengthSummary, PanelProvenanceCounts};
use crate::diamond;
use crate::genomic;
use crate::hmmer;
use crate::mafft;
use crate::scoring_support::ComponentScores;
use crate::structvar;
use crate::taxonomy;
use crate::taxonomy_support::TaxonomyWarnings;
use crate::RenderContext;

pub struct RenderContextInputs {
    pub stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    pub intrinsic_map: Arc<HashMap<String, (crate::metrics::IntrinsicMetrics, Vec<u8>)>>,
    pub alignment_map: Arc<HashMap<String, mafft::AlignmentMetrics>>,
    pub hmmsum_map: Arc<HashMap<String, hmmer::HmmscanSummary>>,
    pub taxsum_map: Arc<HashMap<String, Option<taxonomy::TaxonomyEvidence>>>,
    pub taxonomy_resolver: Option<Arc<taxonomy::TaxonomyResolver>>,
    pub genomic_map: Option<Arc<HashMap<String, genomic::GenomicMetrics>>>,
    pub scores_map: Arc<HashMap<String, (f64, String)>>,
    pub raw_scores_map: Arc<HashMap<String, f64>>,
    pub comp_map: Arc<HashMap<String, ComponentScores>>,
    pub arch_map: Arc<HashMap<String, f64>>,
    pub len_map: Arc<HashMap<String, LengthSummary>>,
    pub orphan_map: Arc<HashMap<String, hmmer::OrphanAnalysis>>,
    pub structvar_map: Arc<HashMap<String, structvar::StructVar>>,
    pub panel_prov_map: Arc<HashMap<String, PanelProvenanceCounts>>,
    pub cov_delta_thresh: f64,
    pub taxonomy_enabled: bool,
    pub taxonomy_warnings: TaxonomyWarnings,
    pub orphan_analysis_enabled: bool,
    pub csv_verbose: bool,
    pub mafft_missing_exon_thresh: usize,
    pub mafft_retained_intron_thresh: usize,
    pub features_string: String,
    pub extensions: LoadedExtensions,
    pub rnaseq_map: Arc<HashMap<String, crate::rnaseq::RnaseqMetrics>>,
    pub rnaseq_enabled: bool,
    pub th_high: f64,
    pub th_med: f64,
}

pub fn build_render_context(inputs: RenderContextInputs) -> Arc<RenderContext> {
    Arc::new(RenderContext {
        stats: inputs.stats,
        intrinsic_map: inputs.intrinsic_map,
        alignment_map: inputs.alignment_map,
        hmmsum_map: inputs.hmmsum_map,
        taxsum_map: inputs.taxsum_map,
        taxonomy_resolver: inputs.taxonomy_resolver,
        genomic_map: inputs.genomic_map,
        scores_map: inputs.scores_map,
        raw_scores_map: inputs.raw_scores_map,
        comp_map: inputs.comp_map,
        arch_map: inputs.arch_map,
        len_map: inputs.len_map,
        orphan_map: inputs.orphan_map,
        structvar_map: inputs.structvar_map,
        panel_prov_map: inputs.panel_prov_map,
        cov_delta_thresh: inputs.cov_delta_thresh,
        taxonomy_enabled: inputs.taxonomy_enabled,
        taxonomy_expected_domain: inputs.taxonomy_warnings.expected_domain,
        taxonomy_warn_non_target_min_frac: inputs.taxonomy_warnings.warn_non_target_min_frac,
        taxonomy_warn_non_target_min_hits: inputs.taxonomy_warnings.warn_non_target_min_hits,
        taxonomy_warn_non_target_strong_frac: inputs.taxonomy_warnings.warn_non_target_strong_frac,
        taxonomy_warn_non_target_strong_hits: inputs.taxonomy_warnings.warn_non_target_strong_hits,
        taxonomy_warn_genus_min_frac: inputs.taxonomy_warnings.warn_genus_min_frac,
        taxonomy_warn_genus_min_hits: inputs.taxonomy_warnings.warn_genus_min_hits,
        taxonomy_low_coverage_frac: inputs.taxonomy_warnings.low_coverage_frac,
        orphan_analysis_enabled: inputs.orphan_analysis_enabled,
        csv_verbose: inputs.csv_verbose,
        mafft_missing_exon_thresh: inputs.mafft_missing_exon_thresh,
        mafft_retained_intron_thresh: inputs.mafft_retained_intron_thresh,
        features_string: inputs.features_string,
        plugins: inputs.extensions.plugins,
        rhai_runtime: inputs.extensions.rhai_runtime,
        rnaseq_map: inputs.rnaseq_map,
        rnaseq_enabled: inputs.rnaseq_enabled,
        th_high: inputs.th_high,
        th_med: inputs.th_med,
    })
}
