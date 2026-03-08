use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use clap::{Args, Parser, Subcommand, ValueEnum};
use mimalloc::MiMalloc;
use serde::{Deserialize, Serialize};
mod alignment_support;
mod analyze;
mod debug_outputs;
use alignment_support::{build_alignment_setup, AlignmentOverrideInputs};
use analysis_support::{
    build_features_string, load_extensions, load_genomic_context, load_rnaseq_data, FeatureSummary,
};
mod analysis_support;
mod checkpoint;
mod clusters;
mod consensus;
mod consensus_support;
mod diamond;
mod ecs;
mod explain;
mod genomic;
mod hmmer;
mod hmmer_support;
mod length;
mod mafft;
mod metrics;
mod orchestration;
mod orf;
mod plugins;
mod preflight;
mod prepare;
mod profiles;
mod provenance;
mod refprot;
mod render_support;
mod reporting_support;
mod rhai_rules;
mod rnaseq;
mod scoring;
mod scoring_support;
mod stage_config_support;
mod structvar;
mod structvar_support;
mod taxonomy;
mod taxonomy_support;
use consensus_support::{
    annotate_refprot_hits, build_consensus_outputs, filter_refprot_hits, load_refprot_proteome_map,
    LengthSummary, PanelProvenanceCounts,
};
use debug_outputs::{dump_matches_fasta, write_panel_sources_csv};
use diamond::{
    blastp_once, cluster as diamond_cluster, linclust as diamond_linclust, parse_tsv_stats,
    recluster as diamond_recluster, DiamondConfig, HitSource,
};
use ecs::{
    run_heavy_pipelines, run_render_pipeline, run_scheduler, EcsConfig, GeneMetrics,
    HeavyPipelineConfig, RenderOutputConfig,
};
use hmmer::HmmscanSummary;
use hmmer_support::{
    build_hmmer_setup, postprocess_hmmer, DomainsArchDebugRow, HmmerOverrideInputs,
    HmmerPostProcessInputs,
};
use mafft::{AlignerBackend, AlignmentMetrics};
use metrics::compute_intrinsic;
use orchestration::{
    step_finish, step_start, write_run_metrics, write_run_summary, write_slowest_genes,
};
use plugins::{run_plugin, PluginDefinition, PluginInput};
use preflight::preflight;
use profiles::Profile;
use provenance::filehash_xx64;
use render_support::{build_render_context, resolve_render_max_jobs, RenderContextInputs};
use reporting_support::{
    augment_totals_from_taxonomy, build_base_totals, build_run_metrics as build_run_metrics_obj,
    build_run_summary as build_run_summary_obj, build_steps, write_domains_arch_debug,
    write_structvar_summary, RefprotSelectionConfig, RunSummaryInputs,
};
use rnaseq::RnaseqMetrics;
use scoring::{
    compute_conserved_regions_score, compute_domains_strength_score, compute_genomic_score,
    compute_genomic_score_with_cfg, compute_homology_score, compute_intrinsic_score,
    compute_rnaseq_score, compute_structvar_multiplier, compute_subject_cov_penalty,
    compute_subject_cov_score, compute_taxonomy_score,
};
use scoring_support::ComponentScores;
use stage_config_support::{
    resolve_consensus_and_refprot_config, ConsensusOverrideInputs, RefProtOverrideInputs,
};
use structvar_support::{analyze_structvar, resolve_thresholds, StructVarOverrideInputs};
use taxonomy::{TaxonomyConsensusConfig, TaxonomyEvidence};
use taxonomy_support::{build_taxonomy_setup, TaxonomyWarnings};

pub(crate) const OUTPUT_SCHEMA_VERSION: &str = "1.1";
const PLUGIN_SCHEMA_VERSION: &str = "v1";
pub(crate) const TOOL_NAME: &str = env!("CARGO_PKG_NAME");
pub(crate) const TOOL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Clone)]
pub(crate) struct RenderContext {
    stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    intrinsic_map: Arc<HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>>,
    alignment_map: Arc<HashMap<String, mafft::AlignmentMetrics>>,
    hmmsum_map: Arc<HashMap<String, hmmer::HmmscanSummary>>,
    taxsum_map: Arc<HashMap<String, Option<TaxonomyEvidence>>>,
    taxonomy_resolver: Option<Arc<taxonomy::TaxonomyResolver>>,
    genomic_map: Option<Arc<HashMap<String, genomic::GenomicMetrics>>>,
    scores_map: Arc<HashMap<String, (f64, String)>>,
    raw_scores_map: Arc<HashMap<String, f64>>,
    comp_map: Arc<HashMap<String, ComponentScores>>,
    arch_map: Arc<HashMap<String, f64>>,
    len_map: Arc<HashMap<String, LengthSummary>>,
    orphan_map: Arc<HashMap<String, hmmer::OrphanAnalysis>>,
    structvar_map: Arc<HashMap<String, structvar::StructVar>>,
    panel_prov_map: Arc<HashMap<String, PanelProvenanceCounts>>,
    cov_delta_thresh: f64,
    taxonomy_enabled: bool,
    taxonomy_expected_domain: Option<String>,
    taxonomy_warn_non_target_min_frac: f64,
    taxonomy_warn_non_target_min_hits: usize,
    taxonomy_warn_non_target_strong_frac: f64,
    taxonomy_warn_non_target_strong_hits: usize,
    taxonomy_warn_genus_min_frac: f64,
    taxonomy_warn_genus_min_hits: usize,
    taxonomy_low_coverage_frac: f64,
    orphan_analysis_enabled: bool,
    csv_verbose: bool,
    mafft_missing_exon_thresh: usize,
    mafft_retained_intron_thresh: usize,
    features_string: String,
    plugins: Vec<PluginDefinition>,
    rhai_runtime: Option<rhai_rules::RhaiRuntime>,
    rnaseq_map: Arc<HashMap<String, RnaseqMetrics>>,
    rnaseq_enabled: bool,
    th_high: f64,
    th_med: f64,
}

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[derive(Clone, Debug, Default)]
pub struct ScoreCard {
    pub gene_id: String,
    pub hits_count: usize,
    pub panel_swissprot: usize,
    pub panel_refprot: usize,
    pub panel_cluster: usize,
    pub top_hit: String,
    pub top_bitscore: f64,
    pub top_evalue: f64,
    pub top_qcov: f64,
    pub top_scov: f64,
    pub bitscore_density: f64,
    pub coverage_delta: f64,
    pub coverage_ratio: f64,
    pub subject_cov_score: f64,
    pub subject_cov_penalty: f64,
    pub fusion_split: bool,
    pub structvar_multiplier: f64,
    pub final_score: f64,
    pub homology_score: f64,
    pub intrinsic_score: f64,
    pub taxonomy_score: Option<f64>,
    pub domains_score: Option<f64>,
    pub domains_arch_score: f64,
    pub orphan_domain_score: f64,
    pub length_score: f64,
    pub length_ratio: f64,
    pub length_class: String,
    pub expected_len_min: Option<f64>,
    pub expected_len_max: Option<f64>,
    pub length_in_expected_range: Option<bool>,
    pub length_panel_n: Option<usize>,
    pub conserved_regions_score: f64,
    pub termini_score: f64,
    pub divergence_score: f64,
    pub mafft_enabled: bool,
    pub conserved_fraction: f64,
    pub pairwise_identity: f64,
    pub panel_pairwise_identity: f64,
    pub divergence_ratio: f64,
    pub sequences_aligned: usize,
    pub query_gap_fraction: f64,
    pub gap_run_count: usize,
    pub max_gap_run: usize,
    pub missing_exon_run: usize,
    pub retained_intron_run: usize,
    pub start_concordance: f64,
    pub start_class: String,
    pub end_concordance: f64,
    pub end_class: String,
    pub structvar_class: String,
    pub structvar_gap: Option<usize>,
    pub orphan_status: String,
    pub taxonomy_contamination: Option<f64>,
    pub taxonomy_support: usize,
    pub taxonomy_considered: usize,
    pub taxonomy_support_frac: f64,
    pub consensus_taxon: String,
    pub taxonomy_rank: String,
    pub taxonomy_status: String,
    pub taxonomy_domain: String,
    pub taxonomy_genus: String,
    pub genomic_introns: Option<usize>,
    pub genomic_splice_canonical: Option<usize>,
    pub genomic_splice_noncanonical: Option<usize>,
    pub genomic_splice_weird: Option<usize>,
    pub genomic_score: f64,
    pub rnaseq_score: f64,
    pub rnaseq_tpm: String,
    pub rnaseq_num_reads: String,
    pub classification_base: String,
    pub classification_final: String,
    pub plugin_penalty: f64,
    pub plugin_count: usize,
    pub plugin_names: String,
    pub plugin_scores: String,
    pub plugin_penalties: String,
    pub plugin_metadata: String,
    pub warnings: String,
}

pub(crate) struct RenderedRecord {
    pub(crate) index: usize,
    pub(crate) json_line: String,
    pub(crate) csv_line: String,
    pub(crate) card: Option<ScoreCard>,
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[arg(long)]
    config: Option<String>,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Prepare(PrepareArgs),
    Analyze(Box<AnalyzeArgs>),
    Explain(ExplainArgs),
    TaxonomyCache(TaxonomyCacheArgs),
    RefprotIndex(RefProtIndexArgs),
    TaxonomyCount(TaxonomyCountArgs),
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize, Default)]
enum CalibrationMode {
    #[default]
    Off,
    Percentile,
    Isotonic,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
enum Mode {
    Centroid,
    Members,
    Auto,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
enum AlignmentStrategy {
    Auto,
    Add,
    AddFragments,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
enum AlignerCliBackend {
    #[serde(alias = "mafft")]
    Mafft,
    #[serde(alias = "spoa")]
    Spoa,
}

impl From<AlignerCliBackend> for AlignerBackend {
    fn from(v: AlignerCliBackend) -> Self {
        match v {
            AlignerCliBackend::Mafft => AlignerBackend::Mafft,
            AlignerCliBackend::Spoa => AlignerBackend::Spoa,
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
enum DiamondMode {
    Auto,
    Batch,
    Single,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
enum LogFormat {
    Text,
    Json,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
enum ReportFormat {
    Jsonl,
    Csv,
    Parquet,
    All,
}

impl ReportFormat {
    fn output_config(self, resume: bool) -> RenderOutputConfig {
        match self {
            ReportFormat::Jsonl => RenderOutputConfig {
                jsonl: true,
                csv: false,
                parquet: false,
                resume,
            },
            ReportFormat::Csv => RenderOutputConfig {
                jsonl: false,
                csv: true,
                parquet: false,
                resume,
            },
            ReportFormat::Parquet => RenderOutputConfig {
                jsonl: false,
                csv: false,
                parquet: true,
                resume,
            },
            ReportFormat::All => RenderOutputConfig {
                jsonl: true,
                csv: true,
                parquet: true,
                resume,
            },
        }
    }
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct AnalyzeArgs {
    #[arg(long)]
    fasta: Option<String>,
    #[arg(long)]
    db: Option<String>,
    #[arg(long)]
    config: Option<String>,
    #[arg(long, default_value_t = 16)]
    threads: usize,
    #[arg(long, default_value_t = 50)]
    approx_id: u32,
    #[arg(long, default_value_t = 80)]
    member_cover: u32,
    #[arg(long)]
    out: Option<String>,
    #[arg(long, value_enum)]
    report_format: Option<ReportFormat>,
    #[arg(long, default_value_t = false)]
    resume: bool,
    #[arg(long, value_enum, default_value_t = Mode::Auto)]
    mode: Mode,
    #[arg(long)]
    top: Option<usize>,
    #[arg(long)]
    diamond_bin: Option<String>,
    #[arg(long)]
    hmmscan_bin: Option<String>,
    #[arg(long)]
    mafft_bin: Option<String>,
    #[arg(long)]
    mafft_threads_per_job: Option<usize>,
    #[arg(long)]
    mafft_max_jobs: Option<usize>,
    #[arg(long, default_value_t = false)]
    mafft_fast: bool,
    /// Alignment backend (mafft or spoa)
    #[arg(long, value_enum, default_value_t = AlignerCliBackend::Mafft)]
    aligner: AlignerCliBackend,
    #[arg(long)]
    render_max_jobs: Option<usize>,
    #[arg(long)]
    reference_fasta: Option<String>,
    #[arg(long, default_value_t = 5)]
    alignment_top_hits: usize,
    #[arg(long, value_enum, default_value_t = AlignmentStrategy::Auto)]
    alignment_strategy: AlignmentStrategy,
    #[arg(long, default_value_t = 0.8)]
    conserved_identity_min: f64,
    #[arg(long)]
    alignment_missing_exon: Option<usize>,
    #[arg(long)]
    alignment_retained_intron: Option<usize>,
    #[arg(long, default_value_t = 8)]
    batch_size: usize,
    #[arg(long)]
    enable_taxonomy: bool,
    #[arg(long)]
    taxonomy_min_support: Option<f64>,
    #[arg(long)]
    taxonomy_top_hits: Option<usize>,
    #[arg(long)]
    taxonomy_min_consensus: Option<usize>,
    #[arg(long)]
    taxonomy_coarse_rank_index: Option<usize>,
    #[arg(long)]
    taxonomy_coarse_min_support: Option<f64>,
    #[arg(long)]
    taxonomy_profile_db: Option<String>,
    #[arg(long)]
    taxonomy_cache: Option<String>,
    #[arg(long)]
    taxonomy_taxdump_dir: Option<String>,
    #[arg(long)]
    pfam_metadata: Option<String>,
    #[arg(long)]
    pfam_clans: Option<String>,
    #[arg(long)]
    pfam_db: Option<String>,
    #[arg(long)]
    refprot_db: Option<String>,
    #[arg(long)]
    refprot_trigger_k: Option<usize>,
    #[arg(long)]
    refprot_min_qcov: Option<f64>,
    #[arg(long)]
    refprot_min_scov: Option<f64>,
    #[arg(long)]
    refprot_max_evalue: Option<f64>,
    #[arg(long)]
    refprot_min_pident: Option<f64>,
    #[arg(long)]
    refprot_max_hits: Option<usize>,
    #[arg(long)]
    refprot_proteome_cap: Option<usize>,
    #[arg(long, default_value_t = true)]
    classify_no_data: bool,
    #[arg(long, default_value_t = false)]
    export_high: bool,
    #[arg(long)]
    export_high_path: Option<String>,
    #[arg(long, default_value_t = false)]
    csv_verbose: bool,
    #[arg(long)]
    hmmer_top_n: Option<usize>,
    #[arg(long)]
    hmmer_threads: Option<usize>,
    #[arg(long)]
    hmmer_ievalue: Option<f64>,
    #[arg(long)]
    hmmer_ref_ievalue: Option<f64>,
    #[arg(long, default_value_t = false)]
    disable_orphan_analysis: bool,
    #[arg(long, value_enum, default_value_t = DiamondMode::Auto)]
    diamond_mode: DiamondMode,
    #[arg(long, value_enum, default_value_t = LogFormat::Text)]
    log_format: LogFormat,
    #[arg(long, default_value_t = 0.25)]
    coverage_delta_threshold: f64,
    #[arg(long)]
    diamond_auto_threshold: Option<usize>,
    #[arg(long)]
    dump_matches_best: bool,
    #[arg(long)]
    dump_matches_gene: Option<String>,
    #[arg(long, default_value_t = false)]
    nucleotide: bool,
    #[arg(long, default_value_t = false)]
    dry_run: bool,
    #[arg(long)]
    diamond_max_hsps: Option<usize>,
    #[arg(long, value_enum)]
    calibration_mode: Option<CalibrationMode>,
    // StructVar thresholds (optional overrides)
    #[arg(long)]
    sv_min_hsp_len: Option<usize>,
    #[arg(long)]
    sv_min_hsp_frac: Option<f64>,
    #[arg(long)]
    sv_fusion_min_gap: Option<usize>,
    #[arg(long)]
    sv_dup_max_gap: Option<usize>,
    #[arg(long)]
    sv_split_delta: Option<f64>,
    #[arg(long)]
    sv_min_subject_cov: Option<f64>,
    #[arg(long)]
    sv_orient_majority: Option<f64>,
    #[arg(long, value_enum, default_value_t = Profile::Standard)]
    profile: Profile,
    #[arg(long)]
    gff: Option<String>,
    #[arg(long)]
    genome: Option<String>,
    #[arg(long)]
    plugin: Vec<String>,
    #[arg(long)]
    rhai: Vec<String>,
    /// RNA-seq expression data file (quant table or expression support)
    #[arg(long, value_name = "FILE")]
    rnaseq_file: Option<String>,
    /// Minimum TPM threshold for expression support
    #[arg(long, value_name = "FLOAT", default_value = "1.0")]
    rnaseq_min_tpm: f64,
}

#[derive(Args, Debug, Clone)]
struct ExplainArgs {
    gene_id: String,
    #[arg(long, default_value_t = String::from("results"))]
    out: String,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct PrepareArgs {
    #[arg(long, default_value_t = String::from("uniprot_sprot.fasta.gz"))]
    fasta: String,
    #[arg(long, default_value_t = String::from("uniprot_sprot.dmnd"))]
    db_out: String,
    #[arg(long, default_value_t = 16)]
    threads: usize,
    #[arg(long, default_value_t = 50)]
    approx_id: u32,
    #[arg(long, default_value_t = 80)]
    member_cover: u32,
    #[arg(long)]
    diamond_bin: Option<String>,
    #[arg(long, default_value_t = false)]
    resume: bool,
    #[arg(long, value_enum, default_value_t = LogFormat::Text)]
    log_format: LogFormat,
    #[arg(long, default_value_t = 4)]
    refprot_workers: usize,
    #[arg(long, default_value_t = 3)]
    download_retries: usize,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct TaxonomyCacheArgs {
    #[arg(long)]
    input: String,
    #[arg(long)]
    output: String,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct RefProtIndexArgs {
    #[arg(long, default_value_t = String::from("share/uniprot/reference_proteomes/README"))]
    readme: String,
    #[arg(long, default_value_t = String::from("share/taxonomy/new_taxdump"))]
    taxdump_dir: String,
    #[arg(long, default_value_t = String::from("Aves"))]
    rank_name: String,
    #[arg(long, default_value_t = 8782)]
    rank_taxid: u32,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct TaxonomyCountArgs {
    #[arg(long, default_value_t = String::from("share/uniprot/reference_proteomes/README"))]
    readme: String,
    #[arg(long, default_value_t = String::from("share/taxonomy/new_taxdump"))]
    taxdump_dir: String,
    #[arg(long)]
    taxid: Option<u32>,
    #[arg(long)]
    name: Option<String>,
    #[arg(long)]
    rank: Option<String>,
    #[arg(long, default_value_t = 20)]
    top: usize,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(dead_code)]
struct FileConfig {
    fasta: Option<String>,
    db: Option<String>,
    threads: Option<usize>,
    approx_id: Option<u32>,
    member_cover: Option<u32>,
    out: Option<String>,
    report_format: Option<ReportFormat>,
    mode: Option<Mode>,
    top: Option<usize>,
    diamond_bin: Option<String>,
    hmmscan_bin: Option<String>,
    scoring: Option<ScoringConfigOverride>,
    reference_fasta: Option<String>,
    mafft_bin: Option<String>,
    alignment_top_hits: Option<usize>,
    alignment_strategy: Option<AlignmentStrategy>,
    conserved_identity_min: Option<f64>,
    alignment_missing_exon: Option<usize>,
    alignment_retained_intron: Option<usize>,
    mafft_threads_per_job: Option<usize>,
    mafft_max_jobs: Option<usize>,
    mafft_fast: Option<bool>,
    mafft_backend: Option<String>,
    render_max_jobs: Option<usize>,
    batch_size: Option<usize>,
    taxonomy: Option<TaxonomyConfigOverride>,
    taxonomy_cache: Option<String>,
    taxonomy_taxdump_dir: Option<String>,
    pfam_metadata: Option<String>,
    pfam_clans: Option<String>,
    pfam_db: Option<String>,
    diamond_mode: Option<DiamondMode>,
    hmmer: Option<HmmerConfigOverride>,
    diamond: Option<DiamondConfigOverride>,
    consensus: Option<ConsensusConfigOverride>,
    refprot: Option<RefProtConfigOverride>,
    structvar: Option<StructVarConfigOverride>,
    export_high: Option<bool>,
    export_high_path: Option<String>,
    calibration: Option<CalibrationConfigOverride>,
    rhai: Option<Vec<String>>,
    rnaseq: Option<RnaseqConfig>,
}

#[derive(Debug, Clone, Default, Deserialize)]
struct RnaseqConfig {
    enabled: Option<bool>,
    file: Option<String>,
    min_tpm: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct TaxonomyConfigOverride {
    enabled: Option<bool>,
    min_support: Option<f64>,
    top_hits: Option<usize>,
    min_consensus: Option<usize>,
    coarse_rank_index: Option<usize>,
    coarse_min_support: Option<f64>,
    expected_domain: Option<String>,
    warn_non_target_min_frac: Option<f64>,
    warn_non_target_min_hits: Option<usize>,
    warn_non_target_strong_frac: Option<f64>,
    warn_non_target_strong_hits: Option<usize>,
    warn_genus_min_frac: Option<f64>,
    warn_genus_min_hits: Option<usize>,
    low_coverage_frac: Option<f64>,
    profile_db: Option<String>,
    #[serde(alias = "cache_path")]
    cache_path: Option<String>,
    #[serde(alias = "taxdump_dir")]
    taxdump_dir: Option<String>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct HmmerConfigOverride {
    top_n: Option<usize>,
    threads: Option<usize>,
    ievalue: Option<f64>,
    ref_ievalue: Option<f64>,
    orphan_analysis: Option<bool>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct DiamondConfigOverride {
    auto_threshold: Option<usize>,
    max_hsps: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct StructVarConfigOverride {
    min_hsp_len: Option<usize>,
    min_hsp_frac: Option<f64>,
    fusion_min_gap: Option<usize>,
    dup_max_gap: Option<usize>,
    split_delta: Option<f64>,
    min_subject_cov: Option<f64>,
    orient_majority: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ConsensusConfigOverride {
    min_hits: Option<usize>,
    max_panel: Option<usize>,
    filt_qcov: Option<f64>,
    filt_scov: Option<f64>,
    filt_evalue: Option<f64>,
    filt_pident: Option<f64>,
    redundancy_pident: Option<f64>,
    max_high_identity: Option<usize>,
    len_ratio_tolerance: Option<f64>,
    backfill_enabled: Option<bool>,
    backfill_min_primary_hits: Option<usize>,
    backfill_max_added: Option<usize>,
    refprot_proteome_cap: Option<usize>,
    diversity_rank_index: Option<usize>,
    diversity_rank_cap: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct CalibrationConfigOverride {
    mode: Option<CalibrationMode>,
    min_samples: Option<usize>,
    min_unique: Option<usize>,
}

#[derive(Debug, Clone, Copy)]
struct CalibrationSettings {
    mode: CalibrationMode,
    min_samples: usize,
    min_unique: usize,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct RefProtConfigOverride {
    enabled: Option<bool>,
    base_dir: Option<String>,
    readme_path: Option<String>,
    taxon_scope_rank: Option<String>, // e.g., "class", "order"
    max_scopes: Option<usize>,        // how many dominant taxa to consider
    max_proteomes: Option<usize>,     // cap selected proteomes
    trigger_k: Option<usize>,
    min_qcov: Option<f64>,
    min_scov: Option<f64>,
    max_evalue: Option<f64>,
    min_pident: Option<f64>,
    max_hits: Option<usize>,
    proteome_cap: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringConfigOverride {
    #[serde(default)]
    weights: HashMap<String, f64>,
    thresholds: Option<ScoringThresholds>,
    caps: Option<ScoringCapsConfigOverride>,
    genomic: Option<ScoringGenomicConfigOverride>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringThresholds {
    high: Option<f64>,
    medium: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringGenomicConfigOverride {
    min_canonical: Option<f64>,
    max_noncanonical: Option<f64>,
    max_weird: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringCapsConfigOverride {
    /// Maximum raw score allowed when `structvar.classification == FusionPossible`.
    structvar_fusion_max: Option<f64>,
    /// Maximum raw score allowed when `structvar.classification == SplitPossible`.
    structvar_split_max: Option<f64>,
    /// Maximum raw score allowed when `structvar.classification == InternalDuplicationPossible`.
    structvar_dup_max: Option<f64>,
}

#[derive(Debug, Clone)]
struct EffectiveConfig {
    fasta: String,
    db: String,
    out: String,
    threads: usize,
    diamond_bin: String,
    reference_fasta: Option<String>,
    refprot_db: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Prepare(p) => prepare::run_prepare(p),
        Commands::Analyze(args) => analyze::run_analyze(cli.config.as_deref(), *args),
        Commands::Explain(args) => {
            explain::explain_gene(&args.out, &args.gene_id)?;
            Ok(())
        }
        Commands::TaxonomyCache(t) => {
            let count = taxonomy::write_cache_from_fasta(&t.input, &t.output)
                .map_err(|e| format!("taxonomy cache failed: {}", e))?;
            eprintln!("wrote {} taxonomy entries to {}", count, t.output);
            Ok(())
        }
        Commands::RefprotIndex(p) => {
            // Build resolver from taxdump only
            let resolver =
                taxonomy::TaxonomyResolver::from_sources(None, None, Some(&p.taxdump_dir))
                    .map_err(|e| format!("taxonomy setup failed: {}", e))?
                    .ok_or("failed to build taxonomy resolver from taxdump")?;
            let entries = refprot::parse_readme(&p.readme)
                .map_err(|e| format!("refprot parse failed: {}", e))?;
            // If rank_name provided and resolves, override rank_taxid
            let target_tid = resolver
                .find_taxid_by_name_exact(&p.rank_name)
                .unwrap_or(p.rank_taxid);
            let aves = refprot::select_by_taxon(&entries, &[target_tid], &resolver, usize::MAX);
            println!(
                "reference_proteomes_total\t{}\nreference_proteomes_{}\t{}",
                entries.len(),
                p.rank_name,
                aves.len()
            );
            // Write a TSV cache for quick reuse
            let outp = std::path::Path::new("share/uniprot/reference_proteomes/refprot_index.tsv");
            if let Some(parent) = outp.parent() {
                std::fs::create_dir_all(parent).ok();
            }
            let mut w = std::fs::File::create(outp)?;
            writeln!(w, "proteome_id\ttaxid\torganism\tis_{}", p.rank_name)?;
            let ave_set: std::collections::HashSet<String> =
                aves.iter().map(|e| e.proteome_id.clone()).collect();
            for e in entries {
                let is = ave_set.contains(&e.proteome_id);
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}",
                    e.proteome_id,
                    e.taxid,
                    e.organism,
                    if is { 1 } else { 0 }
                )?;
            }
            Ok(())
        }
        Commands::TaxonomyCount(args) => {
            run_taxonomy_count(args)?;
            Ok(())
        }
    }
}

fn resolve_effective_config(
    file: &FileConfig,
    args: &AnalyzeArgs,
) -> Result<EffectiveConfig, Box<dyn std::error::Error>> {
    let fasta = args
        .fasta
        .clone()
        .or_else(|| file.fasta.clone())
        .ok_or("--fasta or config fasta required")?;
    let db = args
        .db
        .clone()
        .or_else(|| file.db.clone())
        .ok_or("--db or config db required")?;
    let out = args
        .out
        .clone()
        .or_else(|| file.out.clone())
        .unwrap_or_else(|| "results".to_string());
    let threads = args.threads;
    let diamond_bin = args
        .diamond_bin
        .clone()
        .or_else(|| file.diamond_bin.clone())
        .unwrap_or_else(|| "diamond".to_string());
    let reference_fasta = args
        .reference_fasta
        .clone()
        .or_else(|| file.reference_fasta.clone());
    let refprot_db = if let Some(db) = args.refprot_db.clone() {
        Some(db)
    } else {
        let default = std::path::Path::new("share/refprot/aves/aves_refprot.dmnd");
        if default.exists() {
            Some(default.to_string_lossy().to_string())
        } else {
            None
        }
    };

    Ok(EffectiveConfig {
        fasta,
        db,
        out,
        threads,
        diamond_bin,
        reference_fasta,
        refprot_db,
    })
}

#[derive(Serialize, Clone)]
pub struct Checksums {
    pub fasta_xx64: Option<String>,
    pub db_xx64: Option<String>,
}

fn collect_checksums(cfg: &EffectiveConfig) -> Result<Checksums, Box<dyn std::error::Error>> {
    let fasta = Path::new(&cfg.fasta);
    let db = Path::new(&cfg.db);
    Ok(Checksums {
        fasta_xx64: if fasta.exists() {
            Some(filehash_xx64(fasta)?)
        } else {
            None
        },
        db_xx64: if db.exists() {
            Some(filehash_xx64(db)?)
        } else {
            None
        },
    })
}

fn run_taxonomy_count(args: TaxonomyCountArgs) -> Result<(), String> {
    let entries = refprot::parse_readme(&args.readme)
        .map_err(|e| format!("refprot README parse failed: {}", e))?;
    if entries.is_empty() {
        return Err(format!("no proteomes found in {}", args.readme));
    }
    let resolver = taxonomy::TaxonomyResolver::from_sources(None, None, Some(&args.taxdump_dir))
        .map_err(|e| format!("taxonomy resolver setup failed: {}", e))?
        .ok_or_else(|| {
            format!(
                "taxonomy data missing; ensure {} contains nodes.dmp/names.dmp",
                args.taxdump_dir
            )
        })?;
    let target_taxid = if let Some(tid) = args.taxid {
        tid
    } else if let Some(name) = args.name.as_deref() {
        resolver
            .find_taxid_by_name_exact(name)
            .ok_or_else(|| format!("taxonomy name '{}' not found in taxdump", name))?
    } else {
        return Err("taxonomy-count requires --taxid or --name".into());
    };
    let (_, _, target_name) = resolver.reconstruct_lineage_public(target_taxid);
    let selected = refprot::select_by_taxon(&entries, &[target_taxid], &resolver, usize::MAX);
    println!(
        "Taxon: {} ({})",
        target_taxid,
        target_name.unwrap_or_else(|| "unknown".to_string())
    );
    println!("Reference proteomes total: {}", entries.len());
    println!("Proteomes at/under target: {}", selected.len());
    if let Some(rank) = args.rank.as_deref() {
        let buckets = group_proteomes_by_rank(&resolver, &selected, rank);
        if buckets.is_empty() {
            println!(
                "No descendants expose rank '{}' under taxid {}",
                rank, target_taxid
            );
        } else {
            println!("Top {} {} descendants (by proteome count):", args.top, rank);
            for (idx, (taxid, name, count)) in buckets.iter().enumerate() {
                if idx >= args.top {
                    break;
                }
                println!("  {} (taxid {})\t{}", name, taxid, count);
            }
        }
    }
    Ok(())
}

fn group_proteomes_by_rank(
    resolver: &taxonomy::TaxonomyResolver,
    entries: &[refprot::ProteomeEntry],
    rank: &str,
) -> Vec<(u32, String, usize)> {
    let mut counts: HashMap<u32, (usize, String)> = HashMap::new();
    let rank_lower = rank.to_ascii_lowercase();
    for entry in entries {
        let (lineage_names, lineage_ids, _) = resolver.reconstruct_lineage_public(entry.taxid);
        for (tid, name) in lineage_ids.iter().zip(lineage_names.iter()) {
            if let Some(r) = resolver.rank_of(*tid) {
                if r.eq_ignore_ascii_case(&rank_lower) {
                    let entry = counts.entry(*tid).or_insert((0, name.clone()));
                    entry.0 += 1;
                    break;
                }
            }
        }
    }
    let mut out: Vec<(u32, String, usize)> = counts
        .into_iter()
        .map(|(tid, (count, name))| (tid, name, count))
        .collect();
    out.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| a.1.cmp(&b.1)));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecs::GeneMetrics;
    use crate::taxonomy::TaxonomyEvidence;
    use tempfile::TempDir;

    #[test]
    fn export_high_sequences_writes_only_high() {
        let tmp = TempDir::new().unwrap();
        let out_dir = tmp.path().to_str().unwrap();
        let metrics = vec![
            GeneMetrics {
                gene_id: "gene_high".into(),
                length: 10,
                hits: 0,
            },
            GeneMetrics {
                gene_id: "gene_low".into(),
                length: 10,
                hits: 0,
            },
        ];
        let mut seqs = HashMap::new();
        seqs.insert(
            "gene_high".into(),
            (metrics::IntrinsicMetrics::default(), b"MKTAA".to_vec()),
        );
        seqs.insert(
            "gene_low".into(),
            (metrics::IntrinsicMetrics::default(), b"AAAAA".to_vec()),
        );
        let mut scores = HashMap::new();
        scores.insert("gene_high".into(), (0.9, "High".into()));
        scores.insert("gene_low".into(), (0.4, "Low".into()));
        let res = scoring_support::export_high_sequences(out_dir, None, &metrics, &seqs, &scores)
            .expect("export succeeds");
        assert!(res.is_some());
        let (path, count) = res.unwrap();
        assert_eq!(count, 1);
        let fasta = std::fs::read_to_string(path).unwrap();
        assert!(fasta.contains("gene_high"));
        assert!(!fasta.contains("gene_low"));
    }

    #[test]
    fn propagate_transcript_taxonomy_borrows_best() {
        use crate::taxonomy::TaxonomyDetail;
        use crate::taxonomy_support::propagate_transcript_taxonomy;
        let metrics = vec![
            GeneMetrics {
                gene_id: "gene1.t1".into(),
                length: 100,
                hits: 0,
            },
            GeneMetrics {
                gene_id: "gene1.t2".into(),
                length: 100,
                hits: 0,
            },
            GeneMetrics {
                gene_id: "gene2".into(),
                length: 100,
                hits: 0,
            },
        ];
        let mut map: HashMap<String, Option<TaxonomyEvidence>> = HashMap::new();
        let mut ev = TaxonomyEvidence::default();
        ev.detail = TaxonomyDetail::Consensus;
        ev.congruence_score = 0.9;
        ev.contamination_score = 0.1;
        ev.support = 5;
        ev.considered = 5;
        ev.support_fraction = 1.0;
        map.insert("gene1.t1".into(), Some(ev));
        map.insert("gene1.t2".into(), None);
        map.insert("gene2".into(), None);
        propagate_transcript_taxonomy(&mut map, &metrics);
        let borrowed = map.get("gene1.t2").and_then(|v| v.clone()).unwrap();
        assert_eq!(borrowed.detail, TaxonomyDetail::Borrowed);
        assert_eq!(borrowed.support, 0);
        assert!(map.get("gene2").unwrap().is_none());
    }

    #[test]
    fn dbinfo_text_detection_handles_positive_and_negative() {
        let ok = "Database sequences\nTaxon count: 100";
        assert!(prepare::dbinfo_text_has_taxonomy(ok));
        let bad = "Database sequences\nNo taxonomy found";
        assert!(!prepare::dbinfo_text_has_taxonomy(bad));
    }
}

fn resolve_calibration_settings(args: &AnalyzeArgs, file_cfg: &FileConfig) -> CalibrationSettings {
    let mode = args
        .calibration_mode
        .or_else(|| file_cfg.calibration.as_ref().and_then(|c| c.mode))
        .unwrap_or(CalibrationMode::Off);
    let min_samples = file_cfg
        .calibration
        .as_ref()
        .and_then(|c| c.min_samples)
        .unwrap_or(50)
        .max(1);
    let min_unique = file_cfg
        .calibration
        .as_ref()
        .and_then(|c| c.min_unique)
        .unwrap_or(5)
        .max(1);
    CalibrationSettings {
        mode,
        min_samples,
        min_unique,
    }
}

#[allow(clippy::type_complexity)]
fn compute_intrinsic_for_ids(
    fasta_path: &str,
    metrics: &[GeneMetrics],
) -> Result<HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>, Box<dyn std::error::Error>> {
    let wanted: std::collections::HashSet<String> =
        metrics.iter().map(|m| m.gene_id.clone()).collect();
    let mut map: HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)> = HashMap::new();
    let mut reader = needletail::parse_fastx_file(fasta_path)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let gid = id.split_whitespace().next().unwrap_or("").to_string();
        if !wanted.contains(&gid) {
            continue;
        }
        let seq = rec.seq().to_vec();
        let im = compute_intrinsic(&seq);
        map.insert(gid, (im, seq));
    }
    Ok(map)
}

pub(crate) fn render_gene_record(
    index: usize,
    m: &ecs::GeneMetrics,
    ctx: &RenderContext,
) -> Result<RenderedRecord, String> {
    let summary = ctx.stats.get(&m.gene_id);
    let (intrinsic, seq_bytes) = ctx
        .intrinsic_map
        .get(&m.gene_id)
        .map(|t| (t.0.clone(), t.1.clone()))
        .unwrap_or_default();
    let genomic_metrics = ctx.genomic_map.as_ref().and_then(|map| map.get(&m.gene_id));
    let taxonomy_entry = ctx.taxsum_map.get(&m.gene_id).and_then(|x| x.as_ref());
    let panel_prov = ctx
        .panel_prov_map
        .get(&m.gene_id)
        .cloned()
        .unwrap_or_default();

    // Plugins
    let mut plugin_results = Vec::new();
    let mut plugin_penalty = 0.0;
    if !ctx.plugins.is_empty() || ctx.rhai_runtime.is_some() {
        let homology = summary.map(|s| plugins::PluginHomology {
            hits_count: s.count,
            top_hit: s.top_sseqid.clone(),
            top_bitscore: s.top_bitscore,
            top_evalue: s.top_evalue.clone(),
            top_qcov: s.top_qcov,
            top_scov: s.top_scov,
            bitscore_density: if s.top_len > 0 {
                s.top_bitscore / s.top_len as f64
            } else {
                0.0
            },
            coverage_delta: s.coverage_delta,
            coverage_ratio: s.coverage_ratio,
        });
        let intrinsic_snapshot = plugins::PluginIntrinsic {
            ambiguous_fraction: intrinsic.ambiguous_fraction,
            max_homopolymer: intrinsic.max_homopolymer,
            low_complexity_fraction: intrinsic.low_complexity_fraction,
            low_complexity_windows: intrinsic.low_complexity_windows,
            orf_start_score: intrinsic.orf_start_score,
        };
        let taxonomy_snapshot = if ctx.taxonomy_enabled {
            taxonomy_entry.map(|ev| plugins::PluginTaxonomy {
                detail: ev.detail.to_string(),
                congruence_score: ev.congruence_score,
                contamination_score: ev.contamination_score,
                support_fraction: ev.support_fraction,
                support: ev.support,
                considered: ev.considered,
                consensus_rank: ev.consensus_rank.clone(),
                consensus_taxid: ev.consensus.as_ref().map(|c| c.taxid),
                consensus_name: ev.consensus.as_ref().and_then(|c| c.name.clone()),
            })
        } else {
            None
        };
        let panel_snapshot = Some(plugins::PluginPanel {
            swissprot: panel_prov.swissprot,
            refprot: panel_prov.refprot,
            cluster: panel_prov.cluster,
        });
        let genomic_snapshot = genomic_metrics.map(|g| plugins::PluginGenomic {
            introns_total: g.introns_total,
            splice_canonical: g.splice_canonical,
            splice_noncanonical: g.splice_major_noncan + g.splice_minor,
            splice_weird: g.splice_weird,
            intron_len_min: g.intron_len_min,
            intron_len_max: g.intron_len_max,
            intron_len_avg: g.intron_len_avg,
        });
        let plugin_input = PluginInput {
            gene_id: m.gene_id.clone(),
            sequence: String::from_utf8_lossy(&seq_bytes).to_string(),
            homology,
            intrinsic: intrinsic_snapshot,
            taxonomy: taxonomy_snapshot,
            panel: panel_snapshot,
            genomic: genomic_snapshot,
        };
        for p in &ctx.plugins {
            match run_plugin(p, &plugin_input) {
                Ok(res) => {
                    if let Some(pen) = res.penalty {
                        plugin_penalty += pen;
                    }
                    plugin_results.push(res);
                }
                Err(e) => {
                    // Log debug to avoid spam
                    log::debug!("plugin {} failed for {}: {}", p.name, m.gene_id, e);
                }
            }
        }
        if let Some(rt) = &ctx.rhai_runtime {
            for res in rt.run(&plugin_input) {
                if let Some(pen) = res.penalty {
                    plugin_penalty += pen;
                }
                plugin_results.push(res);
            }
        }
    }
    let sanitize_csv_field = |value: String| value.replace([',', '\n', '\r'], " ");
    let plugin_names_raw = plugin_results
        .iter()
        .map(|r| r.name.clone())
        .collect::<Vec<_>>()
        .join("|");
    let plugin_scores_raw = plugin_results
        .iter()
        .map(|r| {
            let score = r.score.unwrap_or(0.0);
            format!("{}={:.4}", r.name, score)
        })
        .collect::<Vec<_>>()
        .join("|");
    let plugin_penalties_raw = plugin_results
        .iter()
        .map(|r| {
            let pen = r.penalty.unwrap_or(0.0);
            format!("{}={:.4}", r.name, pen)
        })
        .collect::<Vec<_>>()
        .join("|");
    let plugin_metadata_raw = plugin_results
        .iter()
        .filter_map(|r| r.metadata.as_ref().map(|m| (r.name.as_str(), m)))
        .map(|(name, meta)| format!("{}={}", name, meta))
        .collect::<Vec<_>>()
        .join("|");
    let plugin_names = sanitize_csv_field(plugin_names_raw);
    let plugin_scores = sanitize_csv_field(plugin_scores_raw);
    let plugin_penalties = sanitize_csv_field(plugin_penalties_raw);
    let plugin_metadata = sanitize_csv_field(plugin_metadata_raw);
    let plugin_count = plugin_results.len();

    let aln = ctx.alignment_map.get(&m.gene_id);
    let hmmsum = ctx.hmmsum_map.get(&m.gene_id);
    let comp_entry = ctx.comp_map.get(&m.gene_id);
    let homology_score = comp_entry
        .map(|c| c.homology)
        .unwrap_or_else(|| compute_homology_score(summary));
    let intrinsic_score = comp_entry
        .map(|c| c.intrinsic)
        .unwrap_or_else(|| compute_intrinsic_score(&intrinsic));
    let taxonomy_score = if ctx.taxonomy_enabled {
        if let Some(val) = comp_entry.and_then(|c| c.taxonomy) {
            Some(val)
        } else {
            taxonomy_entry.and_then(|ev| {
                if ev.considered > 0 || ev.top_hit.is_some() {
                    Some(compute_taxonomy_score(Some(ev)))
                } else {
                    None
                }
            })
        }
    } else {
        None
    };
    let domains_arch_value = ctx.arch_map.get(&m.gene_id).copied();
    let domains_arch_score = domains_arch_value.unwrap_or(0.0);
    let (length_score, _len_z, len_ratio, len_class, len_min, len_max, len_in_range, len_panel_n) =
        ctx.len_map.get(&m.gene_id).cloned().unwrap_or((
            0.0,
            0.0,
            0.0,
            String::new(),
            0.0,
            0.0,
            false,
            0,
        ));
    let orphan_score = if ctx.orphan_analysis_enabled {
        comp_entry.map(|c| c.orphan).unwrap_or_else(|| {
            ctx.orphan_map
                .get(&m.gene_id)
                .map(|oa| oa.score)
                .unwrap_or(1.0)
        })
    } else {
        1.0
    };
    let subject_cov_score = comp_entry
        .map(|c| c.subject_cov)
        .unwrap_or_else(|| compute_subject_cov_score(summary));
    let subject_cov_penalty = compute_subject_cov_penalty(summary);
    let rnaseq_score = if ctx.rnaseq_enabled {
        comp_entry.map(|c| c.rnaseq).unwrap_or(0.0)
    } else {
        0.0
    };
    let rnaseq_metrics = ctx.rnaseq_map.get(&m.gene_id);
    let rnaseq_tpm = rnaseq_metrics
        .and_then(|r| r.tpm)
        .map_or(String::new(), |v| v.to_string());
    let rnaseq_num_reads = rnaseq_metrics
        .and_then(|r| r.num_reads)
        .map_or(String::new(), |v| v.to_string());
    let termini_score = comp_entry
        .and_then(|c| c.termini)
        .or_else(|| {
            aln.and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            })
        })
        .unwrap_or(0.0);
    let conserved_regions_score = comp_entry
        .and_then(|c| c.conserved_regions)
        .or_else(|| {
            aln.and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
                {
                    Some(compute_conserved_regions_score(Some(a)))
                } else {
                    None
                }
            })
        })
        .unwrap_or(0.0);
    let block_conservation_score = aln
        .and_then(|a| {
            if a.mafft_enabled && !a.conserved_blocks.is_empty() {
                Some(scoring::compute_block_conservation_score(Some(a)))
            } else {
                None
            }
        })
        .unwrap_or(1.0);
    let genomic_metrics = ctx.genomic_map.as_ref().and_then(|map| map.get(&m.gene_id));
    let genomic_score = comp_entry
        .map(|c| c.genomic)
        .unwrap_or_else(|| compute_genomic_score(genomic_metrics));
    let (base_final_score, classif_base) = ctx
        .scores_map
        .get(&m.gene_id)
        .cloned()
        .unwrap_or((0.0, "Low".to_string()));

    // Apply plugin penalty
    let final_score = (base_final_score - plugin_penalty).clamp(0.0, 1.0);

    // Compute final classification based on final_score (post-plugin-penalty)
    let classif_final = if final_score >= ctx.th_high {
        "High"
    } else if final_score >= ctx.th_med {
        "Medium"
    } else {
        "Low"
    }
    .to_string();

    let raw_final_score = ctx
        .raw_scores_map
        .get(&m.gene_id)
        .copied()
        .unwrap_or(base_final_score);
    let fusion_split_flag = summary
        .map(|s| s.coverage_delta > ctx.cov_delta_thresh)
        .unwrap_or(false);
    let mut base_warnings: Vec<String> = Vec::new();
    if m.hits == 0 {
        base_warnings.push("No DIAMOND hits".to_string());
    }
    let mut warn_extra: Vec<String> = Vec::new();
    let mut warning_msgs = base_warnings.clone();
    if let Some(a) = aln {
        if a.missing_exon_run >= ctx.mafft_missing_exon_thresh {
            warn_extra.push("MissingExonPossible".into());
            warning_msgs.push("MissingExonPossible".into());
        }
        if a.retained_intron_run >= ctx.mafft_retained_intron_thresh {
            warn_extra.push("RetainedIntronPossible".into());
            warning_msgs.push("RetainedIntronPossible".into());
        }
        // Block conservation warnings
        if !a.missing_blocks.is_empty() {
            let missing_count = a.missing_blocks.len();
            warn_extra.push(format!("MissingConservedBlocks(count={})", missing_count));
            warning_msgs.push(format!("MissingConservedBlocks(count={})", missing_count));
        }
        if !a.extra_blocks.is_empty() {
            let extra_count = a.extra_blocks.len();
            warn_extra.push(format!("ExtraConservedBlocks(count={})", extra_count));
            warning_msgs.push(format!("ExtraConservedBlocks(count={})", extra_count));
        }
    }
    let sv_obj = ctx.structvar_map.get(&m.gene_id);
    let structvar_multiplier = if m.hits > 0 {
        compute_structvar_multiplier(sv_obj)
    } else {
        1.0
    };
    if let Some(sv) = sv_obj {
        match sv.classification.as_str() {
            "FusionPossible" => {
                warn_extra.push("FusionPossible".into());
                warning_msgs.push("FusionPossible".into());
            }
            "SplitPossible" => {
                warn_extra.push("SplitPossible".into());
                warning_msgs.push("SplitPossible".into());
            }
            "InternalDuplicationPossible" => {
                warn_extra.push("InternalDuplicationPossible".into());
                warning_msgs.push("InternalDuplicationPossible".into());
            }
            _ => {}
        }
        for w in &sv.warnings {
            warn_extra.push(w.clone());
            warning_msgs.push(w.clone());
        }
    }
    let taxonomy_json = if ctx.taxonomy_enabled {
        if let Some(Some(ev)) = ctx.taxsum_map.get(&m.gene_id) {
            let mut obj = serde_json::Map::new();
            obj.insert("status".into(), serde_json::Value::String("enabled".into()));
            obj.insert(
                "detail".into(),
                serde_json::Value::String(ev.detail.to_string()),
            );
            obj.insert(
                "congruence_score".into(),
                serde_json::json!(ev.congruence_score),
            );
            obj.insert(
                "contamination_score".into(),
                serde_json::json!(ev.contamination_score),
            );
            obj.insert(
                "support_fraction".into(),
                serde_json::json!(ev.support_fraction),
            );
            obj.insert(
                "support_hits".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.support as u64)),
            );
            obj.insert(
                "considered_hits".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.considered as u64)),
            );
            obj.insert(
                "consensus_depth".into(),
                serde_json::Value::Number(serde_json::Number::from(ev.consensus_depth as u64)),
            );
            if let Some(rank) = &ev.consensus_rank {
                obj.insert(
                    "consensus_rank".into(),
                    serde_json::Value::String(rank.clone()),
                );
            }
            if let Some(top) = &ev.top_hit {
                obj.insert(
                    "taxid".into(),
                    serde_json::Value::Number(serde_json::Number::from(top.taxid)),
                );
                if let Some(name) = &top.name {
                    obj.insert("name".into(), serde_json::Value::String(name.clone()));
                }
                obj.insert("lineage".into(), serde_json::json!(top.lineage));
                obj.insert(
                    "top_hit".into(),
                    serde_json::json!({
                        "taxid": top.taxid,
                        "name": top.name,
                        "lineage": top.lineage,
                    }),
                );
            }
            if let Some(consensus) = &ev.consensus {
                obj.insert(
                    "consensus_taxid".into(),
                    serde_json::Value::Number(serde_json::Number::from(consensus.taxid)),
                );
                if let Some(name) = &consensus.name {
                    obj.insert(
                        "consensus_name".into(),
                        serde_json::Value::String(name.clone()),
                    );
                }
                obj.insert(
                    "consensus_lineage".into(),
                    serde_json::json!(consensus.lineage.clone()),
                );
            }
            serde_json::Value::Object(obj)
        } else {
            serde_json::json!({"status":"enabled","detail":"NoResolver"})
        }
    } else {
        serde_json::json!({"status":"disabled"})
    };
    let domains = hmmsum.map(|d| {
        let score = if let Some(ev) = d.top_evalue {
            let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
            (le / 20.0).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let arch_score = ctx.arch_map.get(&m.gene_id).copied().unwrap_or(0.0);
        serde_json::json!({
            "hits_count": d.hits_count,
            "top_accession": d.top_accession,
            "top_evalue": d.top_evalue,
            "domains_score": score,
            "domains_arch_score": arch_score,
            "hits": d.hits.iter().map(|h| serde_json::json!({
                "target_name": h.target_name,
                "accession": h.accession,
                "evalue": h.evalue,
                "score": h.score,
                "bias": h.bias
            })).collect::<Vec<_>>()
        })
    });
    let length_block = ctx
        .len_map
        .get(&m.gene_id)
        .map(|(s, z, r, c, mn, mx, in_rng, n)| {
            serde_json::json!({
                "length_score": s,
                "length_z": z,
                "length_ratio": r,
                "length_class": c,
                "expected_len_min": mn,
                "expected_len_max": mx,
                "in_expected_range": in_rng,
                "panel_n": n,
            })
        });
    let orphan_entry = if ctx.orphan_analysis_enabled {
        ctx.orphan_map.get(&m.gene_id)
    } else {
        None
    };
    let divergence_score = comp_entry
        .and_then(|c| c.divergence)
        .unwrap_or_else(|| scoring::compute_divergence_score(aln));
    let record = serde_json::json!({
        "gene_id": m.gene_id,
        "taxonomy": taxonomy_json,
        "score_components": {
            "taxonomy": taxonomy_score,
            "homology": homology_score,
            "intrinsic": intrinsic_score,
            "domains": domains_arch_score,
            "domains_strength": compute_domains_strength_score(hmmsum),
            "subject_coverage": subject_cov_score,
            "subject_cov_penalty": subject_cov_penalty,
            "length": length_score,
            "orphan": orphan_score,
            "conserved_regions": conserved_regions_score,
            "block_conservation": block_conservation_score,
            "structvar_multiplier": structvar_multiplier,
            "termini": termini_score,
            "genomic": genomic_score,
            "rnaseq": rnaseq_score,
            "divergence": divergence_score
        },
        "final_score_raw": raw_final_score,
        "base_score": base_final_score,
        "final_score": final_score,
        "classification_base": classif_base.clone(),
        "classification_final": classif_final.clone(),
        "homology": {
            "hits_count": m.hits,
            "top_hit": summary.and_then(|s| s.top_sseqid.clone()),
            "top_bitscore": summary.map(|s| s.top_bitscore),
            "top_evalue": summary.as_ref().map(|s| s.top_evalue.clone()),
            "top_qcov": summary.map(|s| s.top_qcov),
            "top_scov": summary.map(|s| s.top_scov),
            "bitscore_density": summary.map(|s| if s.top_len > 0 {
                s.top_bitscore / s.top_len as f64
            } else { 0.0 }),
            "coverage_delta": summary.map(|s| s.coverage_delta),
            "coverage_ratio": summary.map(|s| s.coverage_ratio),
            "fusion_split_flag": fusion_split_flag,
        },
        "panel_provenance": {
            "swissprot": panel_prov.swissprot,
            "refprot": panel_prov.refprot,
            "cluster": panel_prov.cluster,
        },
        "intrinsic": {
            "ambiguous_fraction": intrinsic.ambiguous_fraction,
            "max_homopolymer": intrinsic.max_homopolymer,
            "low_complexity_fraction": intrinsic.low_complexity_fraction,
            "low_complexity_windows": intrinsic.low_complexity_windows,
            "orf_diagnostics": {
                "start_methionine": intrinsic.start_methionine,
                "alt_start_pos": intrinsic.alt_start_pos,
                "internal_stop_count": intrinsic.internal_stop_count,
                "terminal_stop": intrinsic.terminal_stop,
                "orf_start_score": intrinsic.orf_start_score,
                "orf_score": scoring::compute_orf_score(&intrinsic),
            },
        },
        "alignment": aln.as_ref().map(|a| serde_json::json!({
            "mafft_enabled": a.mafft_enabled,
            "strategy_used": a.strategy_used,
            "conserved_fraction": a.conserved_fraction,
            "pairwise_identity": a.pairwise_identity,
            "sequences_aligned": a.sequences_aligned,
            "query_gap_fraction": a.query_gap_fraction,
            "gap_run_count": a.gap_run_count,
            "max_gap_run": a.max_gap_run,
            "motif_mismatch_fraction": a.motif_mismatch_fraction,
            "start_concordance": a.start_concordance,
            "start_class": a.start_class,
            "end_concordance": a.end_concordance,
            "end_class": a.end_class,
            "missing_exon_run": a.missing_exon_run,
            "retained_intron_run": a.retained_intron_run,
        })),
        "domains": domains,
        "length": length_block,
        "orphan_analysis": orphan_entry.map(|oa| serde_json::json!({
            "status": oa.status.as_str(),
            "score": oa.score,
            "details": oa.details.iter().map(|d| serde_json::json!({
                "accession": d.accession,
                "domain_index": d.domain_index,
                "total_domains": d.total_domains,
                "completeness": d.completeness,
                "hmm_from": d.hmm_from,
                "hmm_to": d.hmm_to,
                "hmm_len": d.hmm_len,
            })).collect::<Vec<_>>()
        })),
        "structvar": sv_obj.map(|sv| serde_json::json!({
            "classification": sv.classification,
            "fusion_possible": sv.fusion_possible,
            "split_possible": sv.split_possible,
            "duplication_possible": sv.duplication_possible,
            "spans": sv.spans,
            "fusion_subjects": sv.fusion_subjects,
            "fusion_gap": sv.fusion_gap,
            "fusion_left_len": sv.fusion_left_len,
            "fusion_right_len": sv.fusion_right_len,
            "fusion_cover_fracs": sv.fusion_cover_fracs,
            "subjects": sv.subjects,
            "subject_warnings": sv.warnings,
        })),
        "genomic": genomic_metrics.map(|g| serde_json::json!({
            "introns_total": g.introns_total,
            "splice_canonical": g.splice_canonical,
            "splice_noncanonical": g.splice_major_noncan + g.splice_minor,
            "splice_weird": g.splice_weird,
            "intron_len_min": g.intron_len_min,
            "intron_len_max": g.intron_len_max,
            "intron_len_avg": g.intron_len_avg,
            "genomic_score": genomic_score,
        })),
        "rnaseq": ctx.rnaseq_map.get(&m.gene_id).map(|r| serde_json::json!({
            "tpm": r.tpm,
            "num_reads": r.num_reads,
            "expression_score": r.expression_score,
        })),
        "warnings": if warn_extra.is_empty() { base_warnings.clone() } else { warn_extra.clone() },
        "plugins": serde_json::json!({
            "total_penalty": plugin_penalty,
            "results": plugin_results.clone(),
        }),
    });

    let json_line = serde_json::to_string(&record).map_err(|e| e.to_string())?;

    let (top_hit, top_bitscore, top_evalue, top_qcov, top_scov, bsd, cov_delta, cov_ratio) =
        if let Some(s) = summary {
            (
                s.top_sseqid.clone().unwrap_or_default(),
                s.top_bitscore,
                s.top_evalue.clone(),
                s.top_qcov,
                s.top_scov,
                if s.top_len > 0 {
                    s.top_bitscore / s.top_len as f64
                } else {
                    0.0
                },
                s.coverage_delta,
                s.coverage_ratio,
            )
        } else {
            (String::new(), 0.0, String::new(), 0.0, 0.0, 0.0, 0.0, 0.0)
        };
    let subject_cov_score_str = format!("{:.3}", subject_cov_score);
    let subject_cov_penalty_str = format!("{:.3}", subject_cov_penalty);

    let domains_score_field = if let Some(d) = hmmsum {
        let score = if let Some(ev) = d.top_evalue {
            let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
            (le / 20.0).clamp(0.0, 1.0)
        } else {
            0.0
        };
        format!("{:.4}", score)
    } else {
        String::new()
    };
    let domains_arch_field = domains_arch_value
        .map(|v| format!("{:.4}", v))
        .unwrap_or_default();
    let orphan_status_str = orphan_entry
        .map(|oa| oa.status.as_str().to_string())
        .unwrap_or_default();
    let orphan_score_field = if ctx.orphan_analysis_enabled {
        orphan_entry
            .map(|oa| format!("{:.4}", oa.score))
            .unwrap_or_default()
    } else {
        String::new()
    };
    let warnings_field = if warning_msgs.is_empty() {
        String::new()
    } else {
        warning_msgs.join(";")
    };
    let (
        mafft_enabled,
        conserved,
        pid,
        panel_pid,
        div_ratio,
        seqs_aln,
        qgap,
        gap_runs,
        max_gap,
        missing_run,
        intron_run,
        start_conc,
        start_class,
        end_conc,
        end_class,
    ) = if let Some(a) = aln {
        (
            a.mafft_enabled,
            a.conserved_fraction,
            a.pairwise_identity,
            a.panel_pairwise_identity,
            a.divergence_ratio,
            a.sequences_aligned,
            a.query_gap_fraction,
            a.gap_run_count,
            a.max_gap_run,
            a.missing_exon_run,
            a.retained_intron_run,
            a.start_concordance,
            a.start_class.clone(),
            a.end_concordance,
            a.end_class.clone(),
        )
    } else {
        (
            false,
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            0.0,
            0,
            0,
            0,
            0,
            0.0,
            String::new(),
            0.0,
            String::new(),
        )
    };

    let (
        structvar_class,
        structvar_gap,
        structvar_left_len,
        structvar_right_len,
        structvar_cov_left,
        structvar_cov_right,
    ) = if let Some(sv) = sv_obj {
        (
            sv.classification.clone(),
            sv.fusion_gap.map(|g| g.to_string()).unwrap_or_default(),
            sv.fusion_left_len
                .map(|g| g.to_string())
                .unwrap_or_default(),
            sv.fusion_right_len
                .map(|g| g.to_string())
                .unwrap_or_default(),
            sv.fusion_cover_fracs
                .map(|(a, _)| format!("{:.3}", a))
                .unwrap_or_default(),
            sv.fusion_cover_fracs
                .map(|(_, b)| format!("{:.3}", b))
                .unwrap_or_default(),
        )
    } else {
        (
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        )
    };
    let taxonomy_score_value = if ctx.taxonomy_enabled {
        taxonomy_entry
            .and_then(|ev| {
                if ev.considered > 0 || ev.top_hit.is_some() {
                    Some(compute_taxonomy_score(Some(ev)))
                } else {
                    None
                }
            })
            .unwrap_or(0.0)
    } else {
        0.0
    };
    let top_evalue_val = summary
        .and_then(|s| s.top_evalue.parse::<f64>().ok())
        .unwrap_or(0.0);
    let domains_strength_score = compute_domains_strength_score(hmmsum);
    let domains_score = hmmsum.map(|_| domains_strength_score);

    let (
        taxonomy_contamination_field,
        taxonomy_support_field,
        taxonomy_considered_field,
        taxonomy_support_frac_field,
        taxonomy_status_str,
        consensus_label,
        taxonomy_rank_field,
    ) = if ctx.taxonomy_enabled {
        if let Some(ev) = taxonomy_entry {
            let label = ev
                .consensus
                .as_ref()
                .map(|c| {
                    if let Some(name) = &c.name {
                        format!("{} ({})", name, c.taxid)
                    } else {
                        c.taxid.to_string()
                    }
                })
                .unwrap_or_default();
            (
                format!("{:.4}", ev.contamination_score),
                ev.support.to_string(),
                ev.considered.to_string(),
                format!("{:.4}", ev.support_fraction),
                ev.detail.to_string(),
                label,
                ev.consensus_rank.clone().unwrap_or_default(),
            )
        } else {
            (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                "NoResolver".to_string(),
                String::new(),
                String::new(),
            )
        }
    } else {
        (
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            "disabled".to_string(),
            String::new(),
            String::new(),
        )
    };
    let (taxonomy_domain, taxonomy_genus) = if ctx.taxonomy_enabled {
        if let (Some(ev), Some(resolver)) = (taxonomy_entry, ctx.taxonomy_resolver.as_ref()) {
            if let Some(consensus) = &ev.consensus {
                let domain = resolver
                    .ancestor_at_rank(consensus.taxid, "domain")
                    .and_then(|tid| resolver.name_of(tid).map(|s| s.to_string()))
                    .unwrap_or_else(|| "Unknown".to_string());
                let genus = resolver
                    .ancestor_at_rank(consensus.taxid, "genus")
                    .and_then(|tid| resolver.name_of(tid).map(|s| s.to_string()))
                    .unwrap_or_default();
                (domain, genus)
            } else {
                (String::new(), String::new())
            }
        } else {
            (String::new(), String::new())
        }
    } else {
        (String::new(), String::new())
    };

    let len_ratio_str = if len_ratio > 0.0 {
        format!("{:.3}", len_ratio)
    } else {
        String::new()
    };
    let len_expected_min_str = if len_panel_n > 0 && len_min > 0.0 {
        format!("{:.0}", len_min)
    } else {
        String::new()
    };
    let len_expected_max_str = if len_panel_n > 0 && len_max > 0.0 {
        format!("{:.0}", len_max)
    } else {
        String::new()
    };
    let len_in_range_str = if len_panel_n > 0 {
        len_in_range.to_string()
    } else {
        String::new()
    };
    let len_panel_n_str = if len_panel_n > 0 {
        len_panel_n.to_string()
    } else {
        String::new()
    };
    let panel_swissprot_field = panel_prov.swissprot.to_string();
    let panel_refprot_field = panel_prov.refprot.to_string();
    let panel_cluster_field = panel_prov.cluster.to_string();
    let orphan_component_str = if ctx.orphan_analysis_enabled {
        if orphan_score_field.is_empty() {
            let orphan_component = comp_entry.map(|c| c.orphan).unwrap_or(orphan_score);
            format!("{:.4}", orphan_component)
        } else {
            orphan_score_field.clone()
        }
    } else {
        String::new()
    };
    let genomic_score_str = if ctx.genomic_map.is_some() {
        format!("{:.4}", genomic_score)
    } else {
        String::new()
    };

    let rnaseq_score_str = if ctx.rnaseq_enabled {
        format!("{:.4}", rnaseq_score)
    } else {
        String::new()
    };

    let warnings_field_clone = warnings_field.clone();
    let csv_line = if ctx.csv_verbose {
        let taxonomy_component_str = if ctx.taxonomy_enabled {
            format!("{:.4}", taxonomy_score.unwrap_or(0.0))
        } else {
            String::new()
        };
        let row = vec![
            m.gene_id.clone(),
            m.hits.to_string(),
            panel_swissprot_field.clone(),
            panel_refprot_field.clone(),
            panel_cluster_field.clone(),
            top_hit.clone(),
            format!("{:.3}", top_bitscore),
            top_evalue.clone(),
            format!("{:.3}", top_qcov),
            format!("{:.3}", top_scov),
            format!("{:.3}", bsd),
            format!("{:.3}", cov_delta),
            format!("{:.3}", cov_ratio),
            subject_cov_score_str.clone(),
            subject_cov_penalty_str.clone(),
            fusion_split_flag.to_string(),
            format!("{:.4}", structvar_multiplier),
            format!("{:.3}", final_score),
            classif_base.clone(),
            classif_final.clone(),
            format!("{:.4}", homology_score),
            format!("{:.4}", intrinsic_score),
            genomic_score_str.clone(),
            taxonomy_component_str,
            rnaseq_score_str.clone(),
            rnaseq_tpm.clone(),
            rnaseq_num_reads.clone(),
            domains_score_field.clone(),
            domains_arch_field.clone(),
            orphan_component_str,
            format!("{:.4}", length_score),
            len_ratio_str,
            len_class.clone(),
            len_expected_min_str,
            len_expected_max_str,
            len_in_range_str,
            len_panel_n_str,
            format!("{:.4}", conserved_regions_score),
            format!("{:.4}", termini_score),
            format!("{:.4}", divergence_score),
            intrinsic.start_methionine.to_string(),
            intrinsic.internal_stop_count.to_string(),
            intrinsic.terminal_stop.to_string(),
            format!("{:.4}", scoring::compute_orf_score(&intrinsic)),
            mafft_enabled.to_string(),
            format!("{:.3}", conserved),
            format!("{:.3}", pid),
            format!("{:.3}", panel_pid),
            format!("{:.3}", div_ratio),
            seqs_aln.to_string(),
            format!("{:.3}", qgap),
            gap_runs.to_string(),
            max_gap.to_string(),
            missing_run.to_string(),
            intron_run.to_string(),
            format!("{:.3}", start_conc),
            start_class.clone(),
            format!("{:.3}", end_conc),
            end_class.clone(),
            structvar_class.clone(),
            structvar_gap.clone(),
            structvar_left_len.clone(),
            structvar_right_len.clone(),
            structvar_cov_left.clone(),
            structvar_cov_right.clone(),
            orphan_status_str.clone(),
            taxonomy_contamination_field.clone(),
            taxonomy_support_field.clone(),
            taxonomy_considered_field.clone(),
            taxonomy_support_frac_field.clone(),
            consensus_label.clone(),
            taxonomy_rank_field.clone(),
            taxonomy_status_str.clone(),
            format!("{:.4}", plugin_penalty),
            plugin_names.clone(),
            plugin_scores.clone(),
            plugin_penalties.clone(),
            plugin_metadata.clone(),
            warnings_field,
        ]
        .join(",");
        row
    } else {
        let row = vec![
            m.gene_id.clone(),
            m.hits.to_string(),
            panel_swissprot_field.clone(),
            panel_refprot_field.clone(),
            panel_cluster_field.clone(),
            top_hit.clone(),
            format!("{:.3}", top_bitscore),
            top_evalue.clone(),
            format!("{:.3}", top_qcov),
            format!("{:.3}", top_scov),
            format!("{:.3}", bsd),
            format!("{:.3}", cov_delta),
            format!("{:.3}", cov_ratio),
            subject_cov_score_str.clone(),
            subject_cov_penalty_str.clone(),
            fusion_split_flag.to_string(),
            format!("{:.4}", structvar_multiplier),
            format!("{:.3}", final_score),
            classif_base.clone(),
            classif_final.clone(),
            mafft_enabled.to_string(),
            format!("{:.3}", conserved),
            format!("{:.3}", pid),
            format!("{:.3}", panel_pid),
            format!("{:.3}", div_ratio),
            seqs_aln.to_string(),
            format!("{:.3}", qgap),
            gap_runs.to_string(),
            max_gap.to_string(),
            domains_score_field.clone(),
            domains_arch_field.clone(),
            orphan_score_field.clone(),
            structvar_class.clone(),
            structvar_gap.clone(),
            structvar_left_len.clone(),
            structvar_right_len.clone(),
            structvar_cov_left.clone(),
            structvar_cov_right.clone(),
            orphan_status_str.clone(),
            genomic_score_str.clone(),
            format!("{:.4}", taxonomy_score_value),
            taxonomy_contamination_field.clone(),
            taxonomy_support_field.clone(),
            taxonomy_considered_field.clone(),
            taxonomy_support_frac_field.clone(),
            consensus_label.clone(),
            taxonomy_rank_field.clone(),
            taxonomy_status_str.clone(),
            format!("{:.4}", plugin_penalty),
            plugin_names.clone(),
            plugin_scores.clone(),
            plugin_penalties.clone(),
            plugin_metadata.clone(),
            warnings_field_clone.clone(),
        ]
        .join(",");
        row
    };

    let card = ScoreCard {
        gene_id: m.gene_id.clone(),
        hits_count: m.hits,
        panel_swissprot: panel_prov.swissprot,
        panel_refprot: panel_prov.refprot,
        panel_cluster: panel_prov.cluster,
        top_hit: top_hit.clone(),
        top_bitscore,
        top_evalue: top_evalue_val,
        top_qcov,
        top_scov,
        bitscore_density: bsd,
        coverage_delta: cov_delta,
        coverage_ratio: cov_ratio,
        subject_cov_score,
        subject_cov_penalty,
        fusion_split: fusion_split_flag,
        structvar_multiplier,
        final_score,
        classification_base: classif_base.clone(),
        classification_final: classif_final.clone(),
        homology_score,
        intrinsic_score,
        taxonomy_score,
        domains_score,
        domains_arch_score,
        orphan_domain_score: comp_entry.map(|c| c.orphan).unwrap_or(orphan_score),
        length_score,
        length_ratio: len_ratio,
        length_class: len_class.clone(),
        expected_len_min: if len_panel_n > 0 { Some(len_min) } else { None },
        expected_len_max: if len_panel_n > 0 { Some(len_max) } else { None },
        length_in_expected_range: if len_panel_n > 0 {
            Some(len_in_range)
        } else {
            None
        },
        length_panel_n: if len_panel_n > 0 {
            Some(len_panel_n)
        } else {
            None
        },
        conserved_regions_score,
        termini_score,
        divergence_score,
        mafft_enabled,
        conserved_fraction: conserved,
        pairwise_identity: pid,
        panel_pairwise_identity: panel_pid,
        divergence_ratio: div_ratio,
        sequences_aligned: seqs_aln,
        query_gap_fraction: qgap,
        gap_run_count: gap_runs,
        max_gap_run: max_gap,
        missing_exon_run: missing_run,
        retained_intron_run: intron_run,
        start_concordance: start_conc,
        start_class: start_class.clone(),
        end_concordance: end_conc,
        end_class: end_class.clone(),
        structvar_class: structvar_class.clone(),
        structvar_gap: sv_obj.and_then(|sv| sv.fusion_gap),
        orphan_status: orphan_status_str.clone(),
        taxonomy_contamination: if ctx.taxonomy_enabled {
            taxonomy_entry.map(|ev| ev.contamination_score)
        } else {
            None
        },
        taxonomy_support: if let Some(ev) = taxonomy_entry {
            ev.support
        } else {
            0
        },
        taxonomy_considered: if let Some(ev) = taxonomy_entry {
            ev.considered
        } else {
            0
        },
        taxonomy_support_frac: if let Some(ev) = taxonomy_entry {
            ev.support_fraction
        } else {
            0.0
        },
        consensus_taxon: consensus_label.clone(),
        taxonomy_rank: taxonomy_rank_field.clone(),
        taxonomy_status: taxonomy_status_str.clone(),
        taxonomy_domain,
        taxonomy_genus,
        genomic_introns: genomic_metrics.map(|g| g.introns_total),
        genomic_splice_canonical: genomic_metrics.map(|g| g.splice_canonical),
        genomic_splice_noncanonical: genomic_metrics
            .map(|g| g.splice_major_noncan + g.splice_minor),
        genomic_splice_weird: genomic_metrics.map(|g| g.splice_weird),
        genomic_score,
        rnaseq_score,
        rnaseq_tpm,
        rnaseq_num_reads,
        plugin_penalty,
        plugin_count,
        plugin_names,
        plugin_scores,
        plugin_penalties,
        plugin_metadata,
        warnings: warnings_field_clone,
    };

    Ok(RenderedRecord {
        index,
        json_line,
        csv_line,
        card: Some(card),
    })
}

// ... (rest of file unchanged for brevity in plan)
