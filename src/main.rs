use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use clap::{Args, Parser, Subcommand, ValueEnum};
use mimalloc::MiMalloc;
use serde::{Deserialize, Serialize};
use std::time::{Duration, Instant};

mod checkpoint;
mod clusters;
mod consensus;
mod diamond;
mod ecs;
mod explain;
mod genomic;
mod hmmer;
mod length;
mod mafft;
mod metrics;
mod orf;
mod plugins;
mod preflight;
mod profiles;
mod provenance;
mod refprot;
mod rhai_rules;
mod scoring;
mod structvar;
mod taxonomy;
use diamond::{
    blastp_once, cluster as diamond_cluster, linclust as diamond_linclust, parse_tsv_stats,
    recluster as diamond_recluster, DiamondConfig, HitSource,
};
use ecs::{
    run_heavy_pipelines, run_render_pipeline, run_scheduler, AlignmentPipelineConfig, EcsConfig,
    GeneMetrics, HeavyPipelineConfig, HmmerPipelineConfig, RenderOutputConfig,
};
use hmmer::HmmscanSummary;
use mafft::{load_sequences_by_ids, AlignerBackend, AlignerConfig, AlignmentMetrics};
use metrics::compute_intrinsic;
use plugins::{run_plugin, PluginDefinition, PluginInput};
use preflight::preflight;
use profiles::Profile;
use provenance::filehash_xx64;
use scoring::{
    compute_conserved_regions_score, compute_domains_strength_score, compute_genomic_score,
    compute_genomic_score_with_cfg, compute_homology_score, compute_intrinsic_score,
    compute_structvar_multiplier, compute_subject_cov_penalty, compute_subject_cov_score,
    compute_taxonomy_score,
};
use taxonomy::{TaxonomyConsensusConfig, TaxonomyDetail, TaxonomyEvidence, TaxonomyResolver};

type DomainsArchDebugRow = (
    String,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
);
type LengthSummary = (f64, f64, f64, String, f64, f64, bool, usize);

pub(crate) const OUTPUT_SCHEMA_VERSION: &str = "1.1";
const PLUGIN_SCHEMA_VERSION: &str = "v1";
pub(crate) const TOOL_NAME: &str = env!("CARGO_PKG_NAME");
pub(crate) const TOOL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Clone)]
struct PanelAggDebugRow {
    gene_id: String,
    subject_id: String,
    qcov: f64,
    scov: f64,
    len_ratio: f64,
    bitscore: f64,
    hit_count: usize,
    selected: bool,
    source: HitSource,
    quality: f64,
    diversity_key: Option<String>,
}

#[derive(Clone, Default)]
struct PanelProvenanceCounts {
    swissprot: usize,
    refprot: usize,
    cluster: usize,
}

#[derive(Clone)]
struct RefProtFallbackConfig {
    trigger_k: usize,
    min_qcov: f64,
    min_scov: f64,
    max_evalue: f64,
    min_pident: f64,
    max_hits: usize,
}

impl Default for RefProtFallbackConfig {
    fn default() -> Self {
        Self {
            trigger_k: 5,
            min_qcov: 0.50,
            min_scov: 0.25,
            max_evalue: 1e-10,
            min_pident: 30.0,
            max_hits: 100,
        }
    }
}

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
    pub classification: String,
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

#[derive(Debug, Clone, Serialize)]
struct StepDuration<'a> {
    name: &'a str,
    seconds: f64,
}

#[derive(Debug, Clone, Serialize)]
struct RunMetrics<'a> {
    schema_version: &'a str,
    steps: Vec<StepDuration<'a>>,
    totals: std::collections::HashMap<&'a str, serde_json::Value>,
}

#[derive(Debug, Clone, Serialize)]
struct RunCounts {
    total_genes: usize,
    rendered_genes: usize,
    high: usize,
    low: usize,
    x: usize,
    high_complete: usize,
    high_fragmented: usize,
    low_novel: usize,
    low_artifact: usize,
}

#[derive(Debug, Clone, Serialize)]
struct RunThroughput {
    total_seconds: f64,
    genes_per_second: f64,
}

#[derive(Debug, Clone, Serialize)]
struct RunFeatures {
    aligner: String,
    report_format: String,
    resume: bool,
    taxonomy_enabled: bool,
    hmmer_enabled: bool,
    alignment_enabled: bool,
    genomic_enabled: bool,
    plugins: usize,
    rules: usize,
    nucleotide: bool,
    calibration: String,
}

#[derive(Debug, Clone, Serialize)]
struct RunSummary<'a> {
    schema_version: &'a str,
    counts: RunCounts,
    throughput: RunThroughput,
    features: RunFeatures,
    timings: Vec<StepDuration<'a>>,
    errors: Vec<String>,
}

fn step_start(name: &str, json_logs: bool) -> Instant {
    if json_logs {
        log::info!("{}", serde_json::json!({"event":"step_start","name":name}));
    } else {
        log::info!("step_start: {}", name);
    }
    Instant::now()
}

fn step_finish(name: &str, start: Instant, json_logs: bool) -> f64 {
    let secs = start.elapsed().as_secs_f64();
    if json_logs {
        log::info!(
            "{}",
            serde_json::json!({"event":"step_finish","name":name,"seconds":format!("{:.3}",secs)})
        );
    } else {
        log::info!("step_finish: {} ({:.3}s)", name, secs);
    }
    secs
}

fn write_run_metrics(
    out_dir: &str,
    metrics: &RunMetrics,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = std::path::Path::new(out_dir).join("run_metrics.json");
    let text = serde_json::to_string_pretty(metrics)?;
    std::fs::write(path, text)?;
    Ok(())
}

fn write_run_summary(
    out_dir: &str,
    summary: &RunSummary,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = std::path::Path::new(out_dir).join("run_summary.json");
    let text = serde_json::to_string_pretty(summary)?;
    std::fs::write(path, text)?;
    Ok(())
}

#[derive(Debug, Clone, Serialize)]
struct SlowGene {
    gene_id: String,
    seconds: f64,
}

#[derive(Debug, Clone, Serialize)]
struct SlowestGenes {
    alignment: Vec<SlowGene>,
    hmmer: Vec<SlowGene>,
}

fn write_slowest_genes(
    out_dir: &str,
    alignment_secs: &HashMap<String, f64>,
    hmmer_secs: &HashMap<String, f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut align: Vec<SlowGene> = alignment_secs
        .iter()
        .map(|(gid, secs)| SlowGene {
            gene_id: gid.clone(),
            seconds: *secs,
        })
        .collect();
    let mut hmmer: Vec<SlowGene> = hmmer_secs
        .iter()
        .map(|(gid, secs)| SlowGene {
            gene_id: gid.clone(),
            seconds: *secs,
        })
        .collect();
    align.sort_by(|a, b| b.seconds.total_cmp(&a.seconds));
    hmmer.sort_by(|a, b| b.seconds.total_cmp(&a.seconds));
    let payload = SlowestGenes {
        alignment: align.into_iter().take(50).collect(),
        hmmer: hmmer.into_iter().take(50).collect(),
    };
    let path = std::path::Path::new(out_dir).join("slowest_genes.json");
    let text = serde_json::to_string_pretty(&payload)?;
    std::fs::write(path, text)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Prepare(p) => {
            prepare_cmd(p)?;
            Ok(())
        }
        Commands::Analyze(args) => {
            let args = *args;
            let config_path = args.config.as_deref().or(cli.config.as_deref());
            let mut file_cfg = if let Some(path) = config_path {
                let text = fs::read_to_string(path)?;
                toml::from_str::<FileConfig>(&text)?
            } else {
                FileConfig::default()
            };

            // Apply profile overrides to scoring config
            let mut scoring = file_cfg.scoring.unwrap_or_default();
            args.profile.apply_to(&mut scoring);
            file_cfg.scoring = Some(scoring);

            let mut cfg = resolve_effective_config(&file_cfg, &args)?;
            let calibration = resolve_calibration_settings(&args, &file_cfg);

            if args.dry_run {
                print_scoring_rubric(&file_cfg.scoring, calibration, args.classify_no_data);
                return Ok(());
            }

            fs::create_dir_all(&cfg.out)?;

            let report_format = args
                .report_format
                .or(file_cfg.report_format)
                .unwrap_or(ReportFormat::All);
            let mut rhai_paths = file_cfg.rhai.clone().unwrap_or_default();
            if !args.rhai.is_empty() {
                rhai_paths.extend(args.rhai.clone());
            }

            // Pre-flight: tool versions
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

            // If nucleotide mode, translate to protein first
            if args.nucleotide {
                let prot_fa = orf::translate_nt_fasta_to_protein(&cfg.fasta)?;
                cfg.fasta = prot_fa;
            }

            // DIAMOND pre-run
            let dia_cfg = DiamondConfig {
                bin: cfg.diamond_bin.clone(),
                db: cfg.db.clone(),
                query_fasta: cfg.fasta.clone(),
                threads: cfg.threads,
                out_dir: cfg.out.clone(),
                out_name: "diamond.blastp.tsv".to_string(),
                retries: 2,
                max_hsps: args
                    .diamond_max_hsps
                    .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.max_hsps))
                    .unwrap_or(5),
            };
            log::info!(
                "diamond: started mode={:?} out_dir={} out_name={}",
                args.diamond_mode,
                cfg.out,
                dia_cfg.out_name
            );
            let log_json = matches!(args.log_format, LogFormat::Json);
            let t_diamond = step_start("diamond", log_json);
            let diamond_tsv = match args.diamond_mode {
                DiamondMode::Single => blastp_once(&dia_cfg)?,
                DiamondMode::Batch => {
                    diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?
                }
                DiamondMode::Auto => {
                    let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
                    // Heuristic: use config/flag threshold, default 200k
                    let threshold = args
                        .diamond_auto_threshold
                        .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold))
                        .unwrap_or(200_000usize);
                    if n > threshold {
                        diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?
                    } else {
                        blastp_once(&dia_cfg)?
                    }
                }
            };
            let diamond_secs = step_finish("diamond", t_diamond, log_json);
            log::info!(
                "diamond: finished out={} seconds={:.3}",
                diamond_tsv.display(),
                diamond_secs
            );
            let refprot_map_path = Path::new("share/refprot/aves/refprot_proteome_map.tsv");
            let refprot_proteome_map = load_refprot_proteome_map(refprot_map_path);

            // Optional refprot fallback DIAMOND pass (best-effort; never fatal)
            let mut refprot_grouped: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
            if let Some(ref_db) = cfg.refprot_db.as_ref() {
                if Path::new(ref_db).exists() {
                    let refprot_ok =
                        if let Err(e) = preflight::check_diamond_db(&cfg.diamond_bin, ref_db) {
                            log::warn!("refprot db preflight failed; skipping refprot: {}", e);
                            false
                        } else {
                            true
                        };
                    if refprot_ok {
                        let ref_db_size = std::fs::metadata(ref_db).map(|m| m.len()).unwrap_or(0);
                        log::info!(
                            "refprot: db={} size_bytes={} fasta={} mode={:?} out_dir={}",
                            ref_db,
                            ref_db_size,
                            cfg.fasta,
                            args.diamond_mode,
                            cfg.out
                        );
                        let t_ref = step_start("diamond_refprot", log_json);
                        let mut ref_cfg = dia_cfg.clone();
                        ref_cfg.out_name = "diamond.refprot.tsv".into();
                        ref_cfg.db = ref_db.clone();
                        let mut used_chunked = false;
                        let ref_tsv_result = match args.diamond_mode {
                            DiamondMode::Batch => {
                                used_chunked = true;
                                diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json)
                            }
                            DiamondMode::Auto => {
                                let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
                                let threshold = args
                                    .diamond_auto_threshold
                                    .or_else(|| {
                                        file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold)
                                    })
                                    .unwrap_or(200_000usize);
                                if n > threshold {
                                    used_chunked = true;
                                    diamond::blastp_chunked(
                                        &ref_cfg,
                                        args.batch_size.max(1),
                                        log_json,
                                    )
                                } else {
                                    blastp_once(&ref_cfg)
                                }
                            }
                            DiamondMode::Single => blastp_once(&ref_cfg),
                        };
                        if let Ok(ref_tsv) = ref_tsv_result {
                            let mut ref_bytes =
                                std::fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                            if ref_bytes == 0 && !used_chunked {
                                log::warn!(
                                    "refprot blastp produced empty output; retrying in chunked mode"
                                );
                                let _ = diamond::blastp_chunked(
                                    &ref_cfg,
                                    args.batch_size.max(1),
                                    log_json,
                                );
                                ref_bytes =
                                    std::fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                            }
                            if ref_bytes == 0 {
                                log::warn!(
                                    "refprot blastp output is empty; continuing without refprot"
                                );
                            } else {
                                // parse grouped (no qlen map needed here)
                                refprot_grouped =
                                    diamond::parse_tsv_grouped(&ref_tsv, None, Some(100))
                                        .unwrap_or_default();
                                annotate_refprot_hits(&mut refprot_grouped, &refprot_proteome_map);
                            }
                        } else if let Err(e) = ref_tsv_result {
                            log::warn!("refprot blastp failed; continuing without refprot: {}", e);
                        }
                        let _ = step_finish("diamond_refprot", t_ref, log_json);
                    }
                } else {
                    log::warn!("refprot db '{}' not found; skipping", ref_db);
                }
            }

            // ECS: compute per-gene metrics from FASTA and DIAMOND
            let t_ecs = step_start("ecs", log_json);
            let mut metrics: Vec<GeneMetrics> = run_scheduler(EcsConfig {
                fasta_path: cfg.fasta.clone(),
                diamond_tsv: diamond_tsv.to_string_lossy().to_string(),
                threads: cfg.threads,
                log_json: matches!(args.log_format, LogFormat::Json),
            });
            let ecs_secs = step_finish("ecs", t_ecs, log_json);

            if args.resume {
                let rendered_ids = load_rendered_gene_ids(&cfg.out);
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

            // Intrinsic metrics per gene
            let t_intrinsic = step_start("intrinsic", log_json);
            let intrinsic_map = Arc::new(compute_intrinsic_for_ids(&cfg.fasta, &metrics)?);
            let intrinsic_secs = step_finish("intrinsic", t_intrinsic, log_json);

            // Build consensus panels per gene from DIAMOND hits
            let t_consensus = step_start("consensus", log_json);
            let qlen_map: HashMap<String, usize> = metrics
                .iter()
                .map(|m| (m.gene_id.clone(), m.length))
                .collect();
            let grouped = diamond::parse_tsv_grouped(&diamond_tsv, Some(&qlen_map), Some(1000))
                .unwrap_or_default();
            let mut cons_cfg = consensus::ConsensusConfig::default();
            if let Some(cfg_override) = file_cfg.consensus.as_ref() {
                if let Some(v) = cfg_override.min_hits {
                    cons_cfg.min_hits = v;
                }
                if let Some(v) = cfg_override.max_panel {
                    cons_cfg.max_panel = v;
                }
                if let Some(v) = cfg_override.filt_qcov {
                    cons_cfg.filt_qcov = v;
                }
                if let Some(v) = cfg_override.filt_scov {
                    cons_cfg.filt_scov = v;
                }
                if let Some(v) = cfg_override.filt_evalue {
                    cons_cfg.filt_evalue = v;
                }
                if let Some(v) = cfg_override.filt_pident {
                    cons_cfg.filt_pident = v;
                }
                if let Some(v) = cfg_override.redundancy_pident {
                    cons_cfg.redundancy_pident = v;
                }
                if let Some(v) = cfg_override.max_high_identity {
                    cons_cfg.max_high_identity = v;
                }
                if let Some(v) = cfg_override.len_ratio_tolerance {
                    cons_cfg.len_ratio_tolerance = v;
                }
                if let Some(v) = cfg_override.backfill_enabled {
                    cons_cfg.backfill_enabled = v;
                }
                if let Some(v) = cfg_override.backfill_min_primary_hits {
                    cons_cfg.backfill_min_primary_hits = v;
                }
                if let Some(v) = cfg_override.backfill_max_added {
                    cons_cfg.backfill_max_added = v;
                }
                if let Some(v) = cfg_override.refprot_proteome_cap {
                    cons_cfg.refprot_proteome_cap = v;
                }
                if let Some(v) = cfg_override.diversity_rank_index {
                    cons_cfg.diversity_rank_index = v;
                }
                if let Some(v) = cfg_override.diversity_rank_cap {
                    cons_cfg.diversity_rank_cap = v;
                }
            }
            let mut refprot_fallback_cfg = RefProtFallbackConfig {
                trigger_k: cons_cfg.min_hits,
                ..Default::default()
            };
            if let Some(cfg_override) = file_cfg.refprot.as_ref() {
                if let Some(v) = cfg_override.trigger_k {
                    refprot_fallback_cfg.trigger_k = v;
                }
                if let Some(v) = cfg_override.min_qcov {
                    refprot_fallback_cfg.min_qcov = v;
                }
                if let Some(v) = cfg_override.min_scov {
                    refprot_fallback_cfg.min_scov = v;
                }
                if let Some(v) = cfg_override.max_evalue {
                    refprot_fallback_cfg.max_evalue = v;
                }
                if let Some(v) = cfg_override.min_pident {
                    refprot_fallback_cfg.min_pident = v;
                }
                if let Some(v) = cfg_override.max_hits {
                    refprot_fallback_cfg.max_hits = v.max(1);
                }
                if let Some(v) = cfg_override.proteome_cap {
                    cons_cfg.refprot_proteome_cap = v;
                }
            }
            if let Some(v) = args.refprot_trigger_k {
                refprot_fallback_cfg.trigger_k = v;
            }
            if let Some(v) = args.refprot_min_qcov {
                refprot_fallback_cfg.min_qcov = v;
            }
            if let Some(v) = args.refprot_min_scov {
                refprot_fallback_cfg.min_scov = v;
            }
            if let Some(v) = args.refprot_max_evalue {
                refprot_fallback_cfg.max_evalue = v;
            }
            if let Some(v) = args.refprot_min_pident {
                refprot_fallback_cfg.min_pident = v;
            }
            if let Some(v) = args.refprot_max_hits {
                refprot_fallback_cfg.max_hits = v.max(1);
            }
            if let Some(v) = args.refprot_proteome_cap {
                cons_cfg.refprot_proteome_cap = v;
            }
            if !refprot_grouped.is_empty() {
                filter_refprot_hits(&mut refprot_grouped, &refprot_fallback_cfg);
            }
            // Optional cluster backfill: try default clusters file
            let cl_path = Path::new("clusters.recluster");
            if cl_path.exists() {
                match clusters::load_cluster_map(cl_path.to_str().unwrap()) {
                    Ok(map) => {
                        cons_cfg.clusters = Some(std::sync::Arc::new(map));
                        log::info!("consensus: loaded cluster map from {}", cl_path.display());
                    }
                    Err(e) => log::warn!("consensus: failed to load cluster map: {}", e),
                }
            }
            let mut panel_map: HashMap<String, Vec<String>> = HashMap::new();
            let mut panel_stats_rows: Vec<(String, consensus::PanelStats)> = Vec::new();
            let mut panel_len_stats: HashMap<String, consensus::LenStats> = HashMap::new();
            let mut panel_agg_rows: Vec<PanelAggDebugRow> = Vec::new();
            let mut backfill_used: u64 = 0;
            let mut backfill_added_total: u64 = 0;
            let mut len_map: HashMap<String, LengthSummary> = HashMap::new();
            let mut panel_prov_map: HashMap<String, PanelProvenanceCounts> = HashMap::new();
            let mut panel_prov_rows: Vec<(String, PanelProvenanceCounts)> = Vec::new();
            let mut refprot_used_panels: u64 = 0;
            let mut taxonomy_hits_map: HashMap<String, Vec<diamond::DiamondHitRow>> =
                HashMap::new();
            // Use trigger_k to decide when to consider refprot fallback
            let trigger_k = refprot_fallback_cfg.trigger_k.max(0);
            for m in &metrics {
                let primary_hits = grouped.get(&m.gene_id);
                let fallback_hits = refprot_grouped.get(&m.gene_id);
                if primary_hits.is_some() || fallback_hits.is_some() {
                    let mut panel_input: Vec<diamond::DiamondHitRow> =
                        primary_hits.cloned().unwrap_or_default();
                    if panel_input.len() < trigger_k {
                        if let Some(extra) = fallback_hits {
                            panel_input.extend_from_slice(extra);
                        }
                    }
                    let mut taxonomy_hits: Vec<diamond::DiamondHitRow> =
                        primary_hits.cloned().unwrap_or_default();
                    if let Some(extra) = fallback_hits {
                        taxonomy_hits.extend_from_slice(extra);
                    }
                    taxonomy_hits.sort_by(|a, b| {
                        b.bitscore
                            .partial_cmp(&a.bitscore)
                            .unwrap_or(Ordering::Equal)
                    });
                    taxonomy_hits_map.insert(m.gene_id.clone(), taxonomy_hits);
                    if panel_input.is_empty() {
                        panel_stats_rows.push((
                            m.gene_id.clone(),
                            consensus::PanelStats {
                                total_hits: 0,
                                ..Default::default()
                            },
                        ));
                        panel_prov_map.insert(m.gene_id.clone(), PanelProvenanceCounts::default());
                        panel_prov_rows.push((m.gene_id.clone(), PanelProvenanceCounts::default()));
                        continue;
                    }
                    let panel_res = consensus::select_panel_with_result(&panel_input, &cons_cfg);
                    let selection = panel_res.selection.clone();
                    panel_len_stats.insert(m.gene_id.clone(), panel_res.len_stats.clone());
                    let selected_set: HashSet<_> = selection.ids.iter().cloned().collect();
                    let mut source_by_id: HashMap<String, HitSource> = HashMap::new();
                    for agg in &panel_res.aggregated_hits {
                        source_by_id
                            .entry(agg.sseqid.clone())
                            .or_insert_with(|| agg.source.clone());
                        panel_agg_rows.push(PanelAggDebugRow {
                            gene_id: m.gene_id.clone(),
                            subject_id: agg.sseqid.clone(),
                            qcov: agg.qcov,
                            scov: agg.scov,
                            len_ratio: agg.len_ratio,
                            bitscore: agg.bitscore,
                            hit_count: agg.hit_count,
                            selected: selected_set.contains(&agg.sseqid),
                            source: agg.source.clone(),
                            quality: agg.quality,
                            diversity_key: agg.diversity_key(cons_cfg.diversity_rank_index),
                        });
                    }
                    let mut prov_counts = PanelProvenanceCounts::default();
                    for sid in &selection.ids {
                        match source_by_id.get(sid) {
                            Some(HitSource::SwissProt) => prov_counts.swissprot += 1,
                            Some(HitSource::RefProt(_)) => prov_counts.refprot += 1,
                            Some(HitSource::Cluster) => prov_counts.cluster += 1,
                            // IDs added by cluster backfill have no DIAMOND row, so they aren't in
                            // `aggregated_hits` and are counted as cluster provenance here.
                            None => prov_counts.cluster += 1,
                        }
                    }
                    let panel_ids = selection.ids.clone();
                    panel_stats_rows.push((m.gene_id.clone(), selection.stats.clone()));
                    if !panel_ids.is_empty() {
                        panel_map.insert(m.gene_id.clone(), panel_ids.clone());
                        // Length consistency against available subject lengths from the selected panel only
                        let id_set: HashSet<String> = panel_ids.iter().cloned().collect();
                        let slens: Vec<usize> = panel_input
                            .iter()
                            .filter(|h| id_set.contains(&h.sseqid))
                            .map(|h| h.slen)
                            .filter(|&x| x > 0)
                            .collect();
                        if slens.len() >= cons_cfg.min_hits {
                            if let Some(lc) = length::compute_length_consistency(m.length, &slens) {
                                len_map.insert(
                                    m.gene_id.clone(),
                                    (
                                        lc.score,
                                        lc.z,
                                        lc.ratio,
                                        lc.class_,
                                        lc.expected_min,
                                        lc.expected_max,
                                        lc.in_expected_range,
                                        lc.panel_n,
                                    ),
                                );
                            }
                        }
                    }
                    if prov_counts.refprot > 0 {
                        refprot_used_panels += 1;
                    }
                    panel_prov_map.insert(m.gene_id.clone(), prov_counts.clone());
                    panel_prov_rows.push((m.gene_id.clone(), prov_counts));
                } else {
                    taxonomy_hits_map.insert(m.gene_id.clone(), Vec::new());
                    panel_stats_rows.push((
                        m.gene_id.clone(),
                        consensus::PanelStats {
                            total_hits: 0,
                            ..Default::default()
                        },
                    ));
                    panel_prov_map.insert(m.gene_id.clone(), PanelProvenanceCounts::default());
                    panel_prov_rows.push((m.gene_id.clone(), PanelProvenanceCounts::default()));
                }
            }
            let consensus_secs = step_finish("consensus", t_consensus, log_json);
            if !panel_stats_rows.is_empty() {
                let panel_dbg = Path::new(&cfg.out).join("panel_debug.csv");
                let mut f = File::create(panel_dbg)?;
                writeln!(f, "gene_id,total_hits,phase,selected,filtered_hits,len_ratio_window_min,len_ratio_window_max,median_len_ratio,len_ratio_min,len_ratio_max,median_pident,high_identity_dropped,backfill_from_clusters,diversity_rank_index,diversity_cap,diversity_keys_used,diversity_skipped")?;
                for (gid, st) in &panel_stats_rows {
                    writeln!(
                        f,
                        "{},{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{},{},{},{}",
                        gid,
                        st.total_hits,
                        st.phase,
                        st.selected,
                        st.filtered_hits,
                        st.len_ratio_window_min,
                        st.len_ratio_window_max,
                        st.median_len_ratio,
                        st.len_ratio_min,
                        st.len_ratio_max,
                        st.median_pident,
                        st.high_identity_dropped,
                        st.backfill_from_clusters,
                        st.diversity_rank_index,
                        st.diversity_cap,
                        st.diversity_keys_used,
                        st.diversity_skipped,
                    )?;
                }
                if !panel_agg_rows.is_empty() {
                    let agg_path = Path::new(&cfg.out).join("panel_agg_debug.csv");
                    let mut af = File::create(agg_path)?;
                    writeln!(
                        af,
                        "gene_id,subject_id,qcov,scov,len_ratio,bitscore,hit_count,selected,source,quality,diversity_key,len_median,len_mad,expected_len,qlen"
                    )?;
                    for row in &panel_agg_rows {
                        let stats = panel_len_stats
                            .get(&row.gene_id)
                            .cloned()
                            .unwrap_or_default();
                        writeln!(
                            af,
                            "{},{},{:.3},{:.3},{:.3},{:.1},{},{},{:?},{:.3},{},{:.3},{},{},{}",
                            row.gene_id,
                            row.subject_id,
                            row.qcov,
                            row.scov,
                            row.len_ratio,
                            row.bitscore,
                            row.hit_count,
                            row.selected,
                            row.source,
                            row.quality,
                            row.diversity_key.clone().unwrap_or_else(|| "".to_string()),
                            stats.median,
                            stats.mad,
                            stats.expected_len,
                            stats.qlen,
                        )?;
                    }
                }
                // Summarize backfill usage for run_metrics
                let mut bf_used = 0u64;
                let mut bf_added = 0u64;
                for (_gid, st) in &panel_stats_rows {
                    if st.backfill_from_clusters > 0 {
                        bf_used += 1;
                        bf_added += st.backfill_from_clusters as u64;
                    }
                }
                backfill_used = bf_used;
                backfill_added_total = bf_added;
            }

            // Structural variation analysis (heuristic)
            let mut structvar_map: HashMap<String, structvar::StructVar> = HashMap::new();
            // Build thresholds from flags/config
            let mut sv_th = structvar::StructVarThresholds::default();
            if let Some(ref svcfg) = file_cfg.structvar {
                if let Some(v) = svcfg.min_hsp_len {
                    sv_th.min_strong_len = v;
                }
                if let Some(v) = svcfg.min_hsp_frac {
                    sv_th.min_strong_frac = v;
                }
                if let Some(v) = svcfg.fusion_min_gap {
                    sv_th.fusion_min_gap = v;
                }
                if let Some(v) = svcfg.dup_max_gap {
                    sv_th.dup_max_gap = v;
                }
                if let Some(v) = svcfg.split_delta {
                    sv_th.split_delta = v;
                }
                if let Some(v) = svcfg.min_subject_cov {
                    sv_th.min_subject_cov = v;
                }
                if let Some(v) = svcfg.orient_majority {
                    sv_th.orient_majority = v;
                }
            }
            if let Some(v) = args.sv_min_hsp_len {
                sv_th.min_strong_len = v;
            }
            if let Some(v) = args.sv_min_hsp_frac {
                sv_th.min_strong_frac = v;
            }
            if let Some(v) = args.sv_fusion_min_gap {
                sv_th.fusion_min_gap = v;
            }
            if let Some(v) = args.sv_dup_max_gap {
                sv_th.dup_max_gap = v;
            }
            if let Some(v) = args.sv_split_delta {
                sv_th.split_delta = v;
            }
            if let Some(v) = args.sv_min_subject_cov {
                sv_th.min_subject_cov = v;
            }
            if let Some(v) = args.sv_orient_majority {
                sv_th.orient_majority = v;
            }

            for m in &metrics {
                if let Some(hits) = grouped.get(&m.gene_id) {
                    let sv = structvar::analyze(hits, &sv_th);
                    structvar_map.insert(m.gene_id.clone(), sv);
                }
            }

            // Optional MAFFT alignment metrics
            let mut alignment_map: HashMap<String, AlignmentMetrics> = HashMap::new();
            let mut alignment_secs: HashMap<String, f64> = HashMap::new();
            let mut align_cfg_opt: Option<AlignmentPipelineConfig> = None;
            let mut t_mafft_opt: Option<Instant> = None;
            if let (Some(ref_fasta), Some(mafft_bin)) =
                (cfg.reference_fasta.as_ref(), args.mafft_bin.as_ref())
            {
                t_mafft_opt = Some(step_start("mafft", log_json));
                let total_threads = cfg.threads.max(1);
                let default_per_job = if total_threads >= 16 {
                    8
                } else if total_threads >= 8 {
                    4
                } else if total_threads >= 4 {
                    2
                } else {
                    1
                };
                let mafft_threads_per_job = args
                    .mafft_threads_per_job
                    .or(file_cfg.mafft_threads_per_job)
                    .unwrap_or(default_per_job)
                    .clamp(1, total_threads);
                let default_workers = (total_threads / mafft_threads_per_job).max(1);
                let requested_workers = args
                    .mafft_max_jobs
                    .or(file_cfg.mafft_max_jobs)
                    .unwrap_or(default_workers)
                    .max(1);
                let mafft_workers = requested_workers.min(default_workers).max(1);
                std::env::set_var("MAFFT_THREADS", mafft_threads_per_job.to_string());
                let mafft_fast = if args.mafft_fast {
                    true
                } else {
                    file_cfg.mafft_fast.unwrap_or(false)
                };
                if mafft_fast {
                    std::env::set_var("MAFFT_FAST", "1");
                } else {
                    std::env::remove_var("MAFFT_FAST");
                }

                let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
                let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                let query_seq_map: HashMap<String, Vec<u8>> = intrinsic_map
                    .iter()
                    .map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone()))
                    .collect();
                let jobs: Vec<(String, Vec<String>)> = metrics
                    .iter()
                    .filter_map(|g| {
                        let ids = panel_map.get(&g.gene_id)?.clone();
                        if ids.len() < cons_cfg.min_hits {
                            return None;
                        }
                        if !query_seq_map.contains_key(&g.gene_id) {
                            return None;
                        }
                        Some((g.gene_id.clone(), ids))
                    })
                    .collect();

                let backend = if let Some(b) = file_cfg.mafft_backend.as_deref() {
                    match b.to_ascii_lowercase().as_str() {
                        "spoa" => AlignerBackend::Spoa,
                        _ => AlignerBackend::Mafft,
                    }
                } else {
                    args.aligner.into()
                };

                if !jobs.is_empty() {
                    log::info!("mafft jobs queued: {}", jobs.len());
                    let aligner_cfg = AlignerConfig {
                        backend,
                        mafft_bin: mafft_bin.to_string(),
                        mafft_fast,
                        mafft_threads_per_job,
                        mafft_max_jobs: mafft_workers,
                    };
                    align_cfg_opt = Some(AlignmentPipelineConfig {
                        aligner: aligner_cfg,
                        jobs,
                        query_map: query_seq_map,
                        ref_map: ref_seqs,
                    });
                } else {
                    log::info!("mafft jobs queued: 0");
                }
            }

            // Optional: dump matches for inspection
            if let Some(ref_fasta) = cfg.reference_fasta.as_ref() {
                let mut target_gene: Option<String> = None;
                if let Some(g) = args.dump_matches_gene.as_ref() {
                    target_gene = Some(g.clone());
                } else if args.dump_matches_best {
                    // pick gene with largest panel size
                    let mut best: Option<(String, usize)> = None;
                    for m in &metrics {
                        let n = panel_map.get(&m.gene_id).map(|v| v.len()).unwrap_or(0);
                        if best.as_ref().map(|b| n > b.1).unwrap_or(true) {
                            best = Some((m.gene_id.clone(), n));
                        }
                    }
                    if let Some((gid, _)) = best {
                        target_gene = Some(gid);
                    }
                }
                if let Some(gid) = target_gene {
                    let ids = panel_map.get(&gid).cloned().unwrap_or_default();
                    if !ids.is_empty() {
                        let seqs = load_sequences_by_ids(ref_fasta, &ids).unwrap_or_default();
                        let path = std::path::Path::new(&cfg.out).join("matches.fasta");
                        let mut f = std::fs::File::create(path)?;
                        use std::io::Write as _;
                        // write query first
                        if let Some((_im, qseq)) = intrinsic_map.get(&gid) {
                            writeln!(f, ">{}", gid)?;
                            writeln!(f, "{}", String::from_utf8_lossy(qseq))?;
                        }
                        for id in &ids {
                            let key = taxonomy::canonical_accession(id);
                            if let Some(s) = seqs.get(&key) {
                                writeln!(f, ">{}", key)?;
                                writeln!(f, "{}", String::from_utf8_lossy(s))?;
                            }
                        }
                    }
                }
            }
            if !panel_prov_rows.is_empty() {
                let prov_path = Path::new(&cfg.out).join("panel_sources.csv");
                let mut f = File::create(prov_path)?;
                writeln!(f, "gene_id,swissprot_count,refprot_count,cluster_count")?;
                for (gid, prov) in &panel_prov_rows {
                    writeln!(
                        f,
                        "{},{},{},{}",
                        gid, prov.swissprot, prov.refprot, prov.cluster
                    )?;
                }
            }

            // Optional HMMER/Pfam domain summary (JSONL)
            let mut hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut hmmer_secs: HashMap<String, f64> = HashMap::new();
            let mut _ref_hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut domains_arch_map: HashMap<String, f64> = HashMap::new();
            let mut domains_arch_dbg: Vec<DomainsArchDebugRow> = Vec::new();
            let mut orphan_map: HashMap<String, hmmer::OrphanAnalysis> = HashMap::new();
            let mut orphan_analysis_enabled = false;
            let mut hmmer_cfg_opt: Option<HmmerPipelineConfig> = None;
            let mut hmmer_timer: Option<Instant> = None;
            let mut hmmer_bin_owned: Option<String> = None;
            let mut hmmer_db_owned: Option<String> = None;
            let mut hmmer_top_n = 5usize;
            let mut hmmer_thread_cap = cfg.threads;
            let mut hmmer_ref_ievalue = None;
            if let (Some(hmm), Some(pfam_db)) = (
                args.hmmscan_bin.as_ref(),
                args.pfam_db
                    .as_ref()
                    .or(file_cfg.pfam_db.as_ref())
                    .or(file_cfg.pfam_db.as_ref()),
            ) {
                let items: Vec<(String, Vec<u8>)> = intrinsic_map
                    .iter()
                    .map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone()))
                    .collect();
                hmmer_top_n = args
                    .hmmer_top_n
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.top_n))
                    .unwrap_or(5);
                hmmer_thread_cap = args
                    .hmmer_threads
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.threads))
                    .unwrap_or(cfg.threads);
                let query_ievalue = args
                    .hmmer_ievalue
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.ievalue));
                hmmer_ref_ievalue = args
                    .hmmer_ref_ievalue
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.ref_ievalue));
                hmmer_bin_owned = Some(hmm.clone());
                hmmer_db_owned = Some(pfam_db.to_string());
                if !items.is_empty() {
                    log::info!("hmmscan jobs queued: {}", items.len());
                    hmmer_timer = Some(step_start("hmmer", log_json));
                    hmmer_cfg_opt = Some(HmmerPipelineConfig {
                        hmmscan_bin: hmm.to_string(),
                        db_path: pfam_db.to_string(),
                        items,
                        threads_per_job: 1,
                        max_jobs: hmmer_thread_cap.max(1),
                        top_n: hmmer_top_n,
                        max_ievalue: query_ievalue,
                    });
                } else {
                    log::info!("hmmscan jobs queued: 0");
                }
                let orphan_cfg = file_cfg
                    .hmmer
                    .as_ref()
                    .and_then(|h| h.orphan_analysis)
                    .unwrap_or(true);
                orphan_analysis_enabled = orphan_cfg && !args.disable_orphan_analysis;
            }

            let mafft_requested = align_cfg_opt.is_some();
            let hmmer_requested = hmmer_cfg_opt.is_some();
            if mafft_requested || hmmer_requested {
                let reserve_threads = cfg.threads.min(2);
                let heavy_results = run_heavy_pipelines(HeavyPipelineConfig {
                    cpu_threads: cfg.threads,
                    reserve_threads,
                    log_json,
                    alignment: align_cfg_opt.take(),
                    hmmer: hmmer_cfg_opt.take(),
                });
                alignment_map = heavy_results.alignment_map;
                hmmsum_map = heavy_results.hmmer_map;
                alignment_secs = heavy_results.alignment_secs;
                hmmer_secs = heavy_results.hmmer_secs;
            }
            if let Some(t) = t_mafft_opt {
                let _ = step_finish("mafft", t, log_json);
            }
            if let Some(t) = hmmer_timer {
                let _ = step_finish("hmmer", t, log_json);
            }

            if let Some(hmm_bin) = hmmer_bin_owned.clone() {
                if orphan_analysis_enabled {
                    orphan_map = hmmsum_map
                        .iter()
                        .filter(|(_, summary)| !summary.hits.is_empty())
                        .map(|(gid, summary)| (gid.clone(), hmmer::analyze_orphan_domains(summary)))
                        .collect();
                } else {
                    orphan_map.clear();
                }
                if let (Some(ref_fasta), Some(db_path)) =
                    (cfg.reference_fasta.as_ref(), hmmer_db_owned.clone())
                {
                    let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
                    let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                    if !ref_seqs.is_empty() {
                        let ref_items: Vec<(String, Vec<u8>)> = ref_seqs.into_iter().collect();
                        if !ref_items.is_empty() {
                            log::info!("hmmscan reference jobs queued: {}", ref_items.len());
                            let reserve_threads = cfg.threads.min(2);
                            let ref_cfg = HmmerPipelineConfig {
                                hmmscan_bin: hmm_bin.clone(),
                                db_path,
                                items: ref_items,
                                threads_per_job: 1,
                                max_jobs: hmmer_thread_cap.max(1),
                                top_n: hmmer_top_n,
                                max_ievalue: hmmer_ref_ievalue,
                            };
                            let ref_results = run_heavy_pipelines(HeavyPipelineConfig {
                                cpu_threads: cfg.threads,
                                reserve_threads,
                                log_json,
                                alignment: None,
                                hmmer: Some(ref_cfg),
                            });
                            _ref_hmmsum_map = ref_results.hmmer_map;
                        }
                    }
                }
                // Optional clan collapse and architecture scoring
                let clan_map = args
                    .pfam_clans
                    .as_deref()
                    .or(file_cfg.pfam_clans.as_deref())
                    .and_then(|p| hmmer::load_pfam_clans(p).ok());
                for g in &metrics {
                    if let Some(qsum) = hmmsum_map.get(&g.gene_id) {
                        let qsum_c = if let Some(ref clans) = clan_map {
                            hmmer::collapse_by_clan(qsum, clans)
                        } else {
                            qsum.clone()
                        };
                        let ref_ids = panel_map.get(&g.gene_id).cloned().unwrap_or_default();
                        // Canonicalize panel IDs to match keys used in reference hmmscan map
                        let ref_ids_canonical: Vec<String> = ref_ids
                            .iter()
                            .map(|rid| taxonomy::canonical_accession(rid))
                            .collect();
                        let mut ref_map_c: HashMap<String, HmmscanSummary> = HashMap::new();
                        for key in &ref_ids_canonical {
                            if let Some(s) = _ref_hmmsum_map.get(key) {
                                let val = if let Some(ref clans) = clan_map {
                                    hmmer::collapse_by_clan(s, clans)
                                } else {
                                    s.clone()
                                };
                                ref_map_c.insert(key.clone(), val);
                            }
                        }
                        // Light diagnostics if nothing matched (helps debugging join issues)
                        if ref_map_c.is_empty() && !ref_ids_canonical.is_empty() {
                            log::debug!(
                                "ref_join_empty: gene={} panel_n={} keys_example={:?} ref_map_total={}",
                                g.gene_id,
                                ref_ids_canonical.len(),
                                &ref_ids_canonical.iter().take(5).collect::<Vec<_>>(),
                                _ref_hmmsum_map.len()
                            );
                        }
                        let dbg = hmmer::domains_architecture_diagnostics(
                            &qsum_c,
                            &ref_ids_canonical,
                            &ref_map_c,
                        );
                        let score = dbg.score;
                        domains_arch_map.insert(g.gene_id.clone(), score);
                        domains_arch_dbg.push((
                            g.gene_id.clone(),
                            dbg.panel_size,
                            dbg.refs_with_domains,
                            dbg.query_domains,
                            dbg.core_count,
                            dbg.accessory_count,
                            dbg.overlap_core,
                            dbg.overlap_accessory,
                            dbg.recall_core,
                            dbg.precision_acc,
                            dbg.extras_pen,
                            dbg.score,
                        ));
                    }
                }
            }

            let alignment_map = Arc::new(alignment_map);
            let hmmsum_map = Arc::new(hmmsum_map);
            let domains_arch_map = Arc::new(domains_arch_map);
            let len_map = Arc::new(len_map);
            let orphan_map = Arc::new(orphan_map);
            let structvar_map = Arc::new(structvar_map);
            let taxonomy_hits_map = Arc::new(taxonomy_hits_map);

            // Emit outputs
            let t_emit = step_start("emit_outputs", log_json);
            let qlen_map: HashMap<String, usize> = metrics
                .iter()
                .map(|m| (m.gene_id.clone(), m.length))
                .collect();
            let stats =
                Arc::new(parse_tsv_stats(&diamond_tsv, Some(&qlen_map)).unwrap_or_default());
            let checksums = collect_checksums(&cfg)?;
            let snapshot =
                build_config_snapshot(&cfg, &args, &file_cfg, &calibration, report_format);
            let plugin_manifest = collect_plugin_manifest(&args.plugin);
            let rule_manifest = collect_rule_manifest(&rhai_paths);
            write_run_manifest(
                &cfg,
                &tools,
                &checksums,
                &snapshot,
                &plugin_manifest,
                &rule_manifest,
            )?;
            // Taxonomy auto-enable: if resolver can be built, enable unless explicitly disabled.
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
            let taxonomy_top_hits = args
                .taxonomy_top_hits
                .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.top_hits))
                .unwrap_or(20)
                .max(1);
            let taxonomy_min_consensus = args
                .taxonomy_min_consensus
                .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_consensus))
                .unwrap_or(5)
                .max(1);
            let taxonomy_min_support = args
                .taxonomy_min_support
                .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_support))
                .unwrap_or(0.75);
            let taxonomy_coarse_rank = args
                .taxonomy_coarse_rank_index
                .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.coarse_rank_index))
                .unwrap_or(1);
            let taxonomy_coarse_support = args
                .taxonomy_coarse_min_support
                .or_else(|| {
                    file_cfg
                        .taxonomy
                        .as_ref()
                        .and_then(|t| t.coarse_min_support)
                })
                .unwrap_or(0.6);
            let tax_cfg = TaxonomyConsensusConfig {
                min_hits: taxonomy_min_consensus,
                top_hits: taxonomy_top_hits,
                min_support: taxonomy_min_support,
                coarse_rank_index: taxonomy_coarse_rank,
                coarse_min_support: taxonomy_coarse_support,
            };
            let resolver = TaxonomyResolver::from_sources(
                cache_path,
                cfg.reference_fasta.as_deref(),
                taxdump_dir,
            )
            .map_err(|e| format!("taxonomy setup failed: {}", e))?;
            let resolver = resolver.map(Arc::new);
            if resolver.is_some() {
                taxonomy_enabled_effective = true;
            }
            let mut taxsum_local: HashMap<String, Option<TaxonomyEvidence>> = HashMap::new();
            if let Some(ref resolver) = resolver {
                for m in &metrics {
                    let hit_ids: Vec<String> = taxonomy_hits_map
                        .get(&m.gene_id)
                        .map(|rows| {
                            rows.iter()
                                .take(tax_cfg.top_hits)
                                .map(|r| r.sseqid.clone())
                                .collect()
                        })
                        .unwrap_or_default();
                    let mut evidence = resolver.summarize_panel(&hit_ids, &tax_cfg);
                    if evidence.top_hit.is_none() {
                        if let Some(acc) = stats.get(&m.gene_id).and_then(|s| s.top_sseqid.clone())
                        {
                            evidence.top_hit = resolver.lookup(&acc);
                        }
                    }
                    taxsum_local.insert(m.gene_id.clone(), Some(evidence));
                }
                propagate_transcript_taxonomy(&mut taxsum_local, &metrics);
            } else {
                for m in &metrics {
                    taxsum_local.insert(m.gene_id.clone(), None);
                }
            }
            let taxsum_map = Arc::new(taxsum_local);
            let taxonomy_expected_domain = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.expected_domain.clone());
            let taxonomy_warn_non_target_min_frac = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_frac)
                .unwrap_or(0.05);
            let taxonomy_warn_non_target_min_hits = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_hits)
                .unwrap_or(200);
            let taxonomy_warn_non_target_strong_frac = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_frac)
                .unwrap_or(0.10);
            let taxonomy_warn_non_target_strong_hits = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_hits)
                .unwrap_or(500);
            let taxonomy_warn_genus_min_frac = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_frac)
                .unwrap_or(0.15);
            let taxonomy_warn_genus_min_hits = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_hits)
                .unwrap_or(300);
            let taxonomy_low_coverage_frac = file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.low_coverage_frac)
                .unwrap_or(0.05);

            let mut genomic_map: Option<Arc<HashMap<String, genomic::GenomicMetrics>>> = None;
            if let (Some(gff), Some(genome)) = (&args.gff, &args.genome) {
                let t_gff = step_start("genomic_context", log_json);
                match genomic::analyze_gff_context(gff, genome) {
                    Ok(map) => {
                        log::info!("genomic context: loaded {} genes", map.len());
                        genomic_map = Some(Arc::new(map));
                    }
                    Err(e) => {
                        log::warn!("genomic context analysis failed: {}", e);
                    }
                }
                let _ = step_finish("genomic_context", t_gff, log_json);
            }

            let t_scoring = step_start("scoring", log_json);
            let (scores_map_inner, raw_scores_inner) = build_scores_map(
                &metrics,
                stats.as_ref(),
                intrinsic_map.as_ref(),
                &file_cfg.scoring,
                taxonomy_enabled_effective,
                hmmer_requested,
                if hmmer_requested {
                    Some(hmmsum_map.as_ref())
                } else {
                    None
                },
                Some(domains_arch_map.as_ref()),
                Some(len_map.as_ref()),
                if orphan_analysis_enabled {
                    Some(orphan_map.as_ref())
                } else {
                    None
                },
                if taxonomy_enabled_effective {
                    Some(taxsum_map.as_ref())
                } else {
                    None
                },
                Some(alignment_map.as_ref()),
                Some(structvar_map.as_ref()),
                genomic_map.as_ref().map(|gm| gm.as_ref()),
                calibration,
                args.classify_no_data,
            );
            let export_high_enabled = if args.export_high {
                true
            } else {
                file_cfg.export_high.unwrap_or(false)
            };
            let export_high_path = args
                .export_high_path
                .clone()
                .or_else(|| file_cfg.export_high_path.clone());
            let comp_map = build_component_scores(
                &metrics,
                stats.as_ref(),
                intrinsic_map.as_ref(),
                taxonomy_enabled_effective,
                if orphan_analysis_enabled {
                    Some(orphan_map.as_ref())
                } else {
                    None
                },
                if args.enable_taxonomy {
                    Some(taxsum_map.as_ref())
                } else {
                    None
                },
                Some(alignment_map.as_ref()),
                Some(len_map.as_ref()),
                genomic_map.as_ref().map(|gm| gm.as_ref()),
                &file_cfg.scoring,
            );
            let _ = step_finish("scoring", t_scoring, log_json);
            if export_high_enabled {
                if let Some((path, count)) = export_high_sequences(
                    &cfg.out,
                    export_high_path.as_deref(),
                    &metrics,
                    intrinsic_map.as_ref(),
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
            let panel_prov_map = Arc::new(panel_prov_map);
            let mafft_missing_exon_threshold = args
                .alignment_missing_exon
                .or(file_cfg.alignment_missing_exon)
                .unwrap_or(30);
            let mafft_retained_intron_threshold = args
                .alignment_retained_intron
                .or(file_cfg.alignment_retained_intron)
                .unwrap_or(30);

            let mut loaded_plugins = Vec::new();
            for p_path in &args.plugin {
                match PluginDefinition::load(p_path) {
                    Ok(p) => {
                        log::info!("loaded plugin: {}", p.name);
                        loaded_plugins.push(p);
                    }
                    Err(e) => {
                        log::warn!("failed to load plugin {}: {}", p_path, e);
                    }
                }
            }
            let mut rhai_runtime = None;
            if !rhai_paths.is_empty() {
                match rhai_rules::RhaiRuntime::load(&rhai_paths) {
                    Ok(rt) => {
                        let names: Vec<_> = rt.rules().iter().map(|r| r.name.clone()).collect();
                        log::info!("loaded rhai rules: {}", names.join(", "));
                        rhai_runtime = Some(rt);
                    }
                    Err(e) => {
                        log::warn!("failed to load rhai rules: {}", e);
                    }
                }
            }

            let mut features_list = vec!["Homology", "Intrinsic"];
            if taxonomy_enabled_effective {
                features_list.push("Taxonomy");
            }
            if genomic_map.is_some() {
                features_list.push("Genomic");
            }
            if !loaded_plugins.is_empty() {
                features_list.push("Plugins");
            }
            if rhai_runtime.is_some() {
                features_list.push("Rules");
            }
            if hmmer_requested {
                features_list.push("Domains");
            }
            if orphan_analysis_enabled {
                features_list.push("Orphan");
            }
            if mafft_requested {
                features_list.push("Alignment");
                features_list.push("Divergence");
            }
            // Add Profile info if possible, but args.profile is available
            let profile_name = format!("{:?}", args.profile);
            let input_label = if args.gff.is_some() || args.genome.is_some() {
                "GFF+Genome"
            } else if args.nucleotide {
                "Nucleotide FASTA"
            } else {
                "Protein FASTA"
            };
            let features_string = format!(
                "Input: {}, Profile: {}, Features: [{}]",
                input_label,
                profile_name,
                features_list.join("+")
            );

            let plugin_count = loaded_plugins.len();
            let rule_count = rhai_paths.len();
            let render_ctx = Arc::new(RenderContext {
                stats: Arc::clone(&stats),
                intrinsic_map: Arc::clone(&intrinsic_map),
                alignment_map: Arc::clone(&alignment_map),
                hmmsum_map: Arc::clone(&hmmsum_map),
                taxsum_map: Arc::clone(&taxsum_map),
                taxonomy_resolver: resolver.clone(),
                genomic_map: genomic_map.clone(),
                scores_map: Arc::clone(&scores_map),
                raw_scores_map: Arc::clone(&raw_scores_map),
                comp_map: Arc::clone(&comp_map),
                arch_map: Arc::clone(&domains_arch_map),
                len_map: Arc::clone(&len_map),
                orphan_map: Arc::clone(&orphan_map),
                structvar_map: Arc::clone(&structvar_map),
                panel_prov_map: Arc::clone(&panel_prov_map),
                cov_delta_thresh: args.coverage_delta_threshold,
                taxonomy_enabled: taxonomy_enabled_effective,
                taxonomy_expected_domain,
                taxonomy_warn_non_target_min_frac,
                taxonomy_warn_non_target_min_hits,
                taxonomy_warn_non_target_strong_frac,
                taxonomy_warn_non_target_strong_hits,
                taxonomy_warn_genus_min_frac,
                taxonomy_warn_genus_min_hits,
                taxonomy_low_coverage_frac,
                orphan_analysis_enabled,
                csv_verbose: args.csv_verbose,
                mafft_missing_exon_thresh: mafft_missing_exon_threshold,
                mafft_retained_intron_thresh: mafft_retained_intron_threshold,
                features_string,
                plugins: loaded_plugins,
                rhai_runtime,
            });
            let render_max_jobs = args
                .render_max_jobs
                .or(file_cfg.render_max_jobs)
                .unwrap_or(cfg.threads.max(1))
                .max(1);
            let render_output = report_format.output_config(args.resume);
            let render_summary = run_render_pipeline(
                &cfg.out,
                &metrics,
                render_ctx,
                render_max_jobs,
                &tools,
                &checksums,
                render_output,
            )?;
            // Emit domains architecture diagnostics CSV
            if !domains_arch_dbg.is_empty() {
                let dbg_path = std::path::Path::new(&cfg.out).join("domains_arch_debug.csv");
                let mut fdbg = std::fs::File::create(dbg_path)?;
                use std::io::Write as _;
                writeln!(fdbg, "gene_id,panel_size,refs_with_domains,query_domains,core_count,accessory_count,overlap_core,overlap_accessory,recall_core,precision_acc,extras_pen,score")?;
                for (gid, ps, rwd, qd, cc, ac, oc, oa, rc, pa, ep, sc) in domains_arch_dbg {
                    writeln!(
                        fdbg,
                        "{},{},{},{},{},{},{},{},{:.4},{:.4},{:.4},{:.4}",
                        gid, ps, rwd, qd, cc, ac, oc, oa, rc, pa, ep, sc
                    )?;
                }
            }
            // Emit structvar summary (counts and simple distributions) for tuning
            if !structvar_map.is_empty() {
                let mut counts: std::collections::HashMap<&str, usize> = Default::default();
                let mut gaps: Vec<usize> = Vec::new();
                let mut covs: Vec<f64> = Vec::new();
                for sv in structvar_map.values() {
                    let k = sv.classification.as_str();
                    *counts.entry(k).or_insert(0) += 1;
                    if let Some(g) = sv.fusion_gap {
                        gaps.push(g);
                    }
                    if let Some((a, b)) = sv.fusion_cover_fracs {
                        covs.push(a);
                        covs.push(b);
                    }
                }
                covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                gaps.sort_unstable();
                let pct = |v: &[usize], p: f64| -> usize {
                    if v.is_empty() {
                        0
                    } else {
                        let i = ((p * (v.len() as f64)).clamp(0.0, (v.len() - 1) as f64)) as usize;
                        v[i]
                    }
                };
                let pctf = |v: &[f64], p: f64| -> f64 {
                    if v.is_empty() {
                        0.0
                    } else {
                        let i = ((p * (v.len() as f64)).clamp(0.0, (v.len() - 1) as f64)) as usize;
                        v[i]
                    }
                };
                let summ = serde_json::json!({
                    "counts": counts,
                    "gap_percentiles": {"p10": pct(&gaps,0.10), "p50": pct(&gaps,0.50), "p90": pct(&gaps,0.90)},
                    "cov_percentiles": {"p10": format!("{:.3}", pctf(&covs,0.10)), "p50": format!("{:.3}", pctf(&covs,0.50)), "p90": format!("{:.3}", pctf(&covs,0.90))},
                    "suggested_cutoffs": {"fusion_min_gap": pct(&gaps,0.10).max(50), "min_subject_cov": pctf(&covs,0.10)},
                });
                std::fs::write(
                    std::path::Path::new(&cfg.out).join("structvar_summary.json"),
                    serde_json::to_string_pretty(&summ)?,
                )?;
            }
            let emit_secs = step_finish("emit_outputs", t_emit, log_json);

            // Write compact run_metrics.json sidecar
            let mut totals: std::collections::HashMap<&str, serde_json::Value> = Default::default();
            let hmmer_genes = hmmsum_map.len() as u64;
            let hmmer_hits_total: u64 = hmmsum_map.values().map(|s| s.hits_count as u64).sum();
            let mafft_alignments = alignment_map.len() as u64;
            totals.insert("total_genes", serde_json::json!(metrics.len()));
            totals.insert("hmmer_genes", serde_json::json!(hmmer_genes));
            totals.insert("hmmer_hits_total", serde_json::json!(hmmer_hits_total));
            totals.insert("mafft_alignments", serde_json::json!(mafft_alignments));
            totals.insert(
                "diamond_mode",
                serde_json::json!(format!("{:?}", args.diamond_mode)),
            );
            let steps = vec![
                StepDuration {
                    name: "diamond",
                    seconds: diamond_secs,
                },
                StepDuration {
                    name: "ecs",
                    seconds: ecs_secs,
                },
                StepDuration {
                    name: "intrinsic",
                    seconds: intrinsic_secs,
                },
                StepDuration {
                    name: "consensus",
                    seconds: consensus_secs,
                },
                StepDuration {
                    name: "emit_outputs",
                    seconds: emit_secs,
                },
            ];
            // Taxonomy + backfill summaries for run_metrics
            if taxonomy_enabled_effective {
                let mut status_counts: HashMap<String, u64> = HashMap::new();
                let mut superkingdom: HashMap<String, u64> = HashMap::new();
                let mut class_counts: HashMap<String, u64> = HashMap::new();
                let mut species_counts: HashMap<String, u64> = HashMap::new();
                let mut class_scope_counts: HashMap<String, u64> = HashMap::new();
                for m in &metrics {
                    if let Some(Some(ev)) = taxsum_map.get(&m.gene_id) {
                        *status_counts.entry(format!("{}", ev.detail)).or_default() += 1;
                        if let Some(cons) = &ev.consensus {
                            // superkingdom at index 1 if available
                            if cons.lineage.len() > 1 {
                                *superkingdom.entry(cons.lineage[1].clone()).or_default() += 1;
                            }
                            // class rank label if present
                            let mut class_label: Option<String> = None;
                            for (i, tid) in cons.lineage_ids.iter().enumerate() {
                                if let Some(r) = resolver.as_ref() {
                                    if r.rank_of(*tid) == Some("class") {
                                        class_label = cons.lineage.get(i).cloned();
                                        break;
                                    }
                                }
                            }
                            if let Some(lbl) = class_label {
                                *class_counts.entry(lbl.clone()).or_default() += 1;
                                *class_scope_counts.entry(lbl).or_default() += 1;
                            }
                            // species: last lineage name if depth >=1
                            if let Some(last) = cons.lineage.last() {
                                *species_counts.entry(last.clone()).or_default() += 1;
                            }
                        }
                    }
                }
                totals.insert(
                    "taxonomy_status_counts",
                    serde_json::to_value(status_counts).unwrap(),
                );
                totals.insert(
                    "taxonomy_superkingdom_top",
                    serde_json::to_value(superkingdom).unwrap(),
                );
                totals.insert(
                    "taxonomy_class_top",
                    serde_json::to_value(class_counts).unwrap(),
                );
                totals.insert(
                    "taxonomy_species_top",
                    serde_json::to_value(species_counts).unwrap(),
                );
                // Optional: derive model proteomes from top class scopes and README mapping
                if let Some(ref cfg_rp) = file_cfg.refprot {
                    if cfg_rp.enabled.unwrap_or(false) {
                        if let (Some(readme), Some(resolver)) =
                            (cfg_rp.readme_path.as_ref(), resolver.as_ref())
                        {
                            if let Ok(entries) = refprot::parse_readme(readme) {
                                let mut v: Vec<(String, u64)> =
                                    class_scope_counts.into_iter().collect();
                                v.sort_by(|a, b| b.1.cmp(&a.1));
                                let topn = cfg_rp.max_scopes.unwrap_or(2).max(1);
                                let top_classes: Vec<String> =
                                    v.into_iter().map(|x| x.0).take(topn).collect();
                                // Map selected class names to taxids by scanning resolver lineages of proteomes
                                let mut class_taxids: Vec<u32> = Vec::new();
                                for e in &entries {
                                    let (lineage, lids, _) =
                                        resolver.reconstruct_lineage_public(e.taxid);
                                    for (i, name) in lineage.iter().enumerate() {
                                        if top_classes.iter().any(|c| c == name) {
                                            if let Some(tid) = lids.get(i) {
                                                class_taxids.push(*tid);
                                            }
                                            break;
                                        }
                                    }
                                }
                                class_taxids.sort();
                                class_taxids.dedup();
                                let selected = refprot::select_by_taxon(
                                    &entries,
                                    &class_taxids,
                                    resolver,
                                    cfg_rp.max_proteomes.unwrap_or(10),
                                );
                                let mut lines = Vec::new();
                                for e in &selected {
                                    lines.push(format!(
                                        "{}\t{}\t{}",
                                        e.proteome_id, e.taxid, e.organism
                                    ));
                                }
                                let _ = std::fs::write(
                                    std::path::Path::new(&cfg.out).join("refprot_selected.txt"),
                                    lines.join("\n"),
                                );
                                totals.insert(
                                    "refprot_selected_count",
                                    serde_json::json!(selected.len()),
                                );
                            }
                        }
                    }
                }
            }
            totals.insert(
                "cluster_backfill_used_panels",
                serde_json::json!(backfill_used),
            );
            totals.insert(
                "cluster_backfill_added_total",
                serde_json::json!(backfill_added_total),
            );
            totals.insert(
                "refprot_panels_with_hits",
                serde_json::json!(refprot_used_panels),
            );

            write_slowest_genes(&cfg.out, &alignment_secs, &hmmer_secs)?;

            let total_secs: f64 = steps.iter().map(|s| s.seconds).sum();
            let genes_per_sec = if total_secs > 0.0 {
                metrics.len() as f64 / total_secs
            } else {
                0.0
            };
            let run_summary = RunSummary {
                schema_version: OUTPUT_SCHEMA_VERSION,
                counts: RunCounts {
                    total_genes: metrics.len(),
                    rendered_genes: render_summary.total,
                    high: render_summary.high,
                    low: render_summary.low,
                    x: render_summary.x,
                    high_complete: render_summary.high_complete,
                    high_fragmented: render_summary.high_fragmented,
                    low_novel: render_summary.low_novel,
                    low_artifact: render_summary.low_artifact,
                },
                throughput: RunThroughput {
                    total_seconds: total_secs,
                    genes_per_second: genes_per_sec,
                },
                features: RunFeatures {
                    aligner: format!("{:?}", args.aligner).to_lowercase(),
                    report_format: format!("{:?}", report_format).to_lowercase(),
                    resume: args.resume,
                    taxonomy_enabled: taxonomy_enabled_effective,
                    hmmer_enabled: hmmer_requested,
                    alignment_enabled: mafft_requested,
                    genomic_enabled: args.gff.is_some(),
                    plugins: plugin_count,
                    rules: rule_count,
                    nucleotide: args.nucleotide,
                    calibration: format!("{:?}", calibration.mode).to_lowercase(),
                },
                timings: steps.clone(),
                errors: Vec::new(),
            };
            write_run_summary(&cfg.out, &run_summary)?;

            let runm = RunMetrics {
                schema_version: OUTPUT_SCHEMA_VERSION,
                steps,
                totals,
            };
            write_run_metrics(&cfg.out, &runm)?;
            Ok(())
        }
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

#[derive(Serialize, Clone)]
struct PluginManifest {
    name: String,
    path: String,
    xx64: Option<String>,
    size_bytes: Option<u64>,
    version: Option<String>,
}

#[derive(Serialize, Clone)]
struct RuleManifest {
    path: String,
    xx64: Option<String>,
}

fn collect_plugin_manifest(paths: &[String]) -> Vec<PluginManifest> {
    paths
        .iter()
        .map(|p| {
            let path = Path::new(p);
            let xx64 = if path.exists() {
                filehash_xx64(path).ok()
            } else {
                None
            };
            let size_bytes = path.metadata().map(|m| m.len()).ok();
            let mut name = path
                .file_stem()
                .or_else(|| path.file_name())
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| "plugin".to_string());
            let mut version = None;
            let sidecar = path.with_extension("json");
            if sidecar.exists() {
                if let Ok(text) = std::fs::read_to_string(&sidecar) {
                    if let Ok(val) = serde_json::from_str::<serde_json::Value>(&text) {
                        if let Some(v) = val.get("name").and_then(|v| v.as_str()) {
                            name = v.to_string();
                        }
                        if let Some(v) = val.get("version").and_then(|v| v.as_str()) {
                            version = Some(v.to_string());
                        }
                    }
                }
            }
            PluginManifest {
                name,
                path: p.clone(),
                xx64,
                size_bytes,
                version,
            }
        })
        .collect()
}

fn collect_rule_manifest(paths: &[String]) -> Vec<RuleManifest> {
    paths
        .iter()
        .map(|p| {
            let path = Path::new(p);
            let xx64 = if path.exists() {
                filehash_xx64(path).ok()
            } else {
                None
            };
            RuleManifest {
                path: p.clone(),
                xx64,
            }
        })
        .collect()
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

fn load_rendered_gene_ids(out_dir: &str) -> std::collections::HashSet<String> {
    use std::io::{BufRead, BufReader};
    let mut out = std::collections::HashSet::new();
    let json_path = Path::new(out_dir).join("qc_report.jsonl");
    if json_path.exists() {
        if let Ok(file) = File::open(&json_path) {
            let reader = BufReader::new(file);
            let mut parsed = 0usize;
            for line in reader.lines() {
                let line = match line {
                    Ok(line) => line,
                    Err(err) => {
                        log::warn!("resume: failed reading qc_report.jsonl: {err}");
                        continue;
                    }
                };
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let Ok(val) = serde_json::from_str::<serde_json::Value>(line) else {
                    log::debug!("resume: skip malformed jsonl line");
                    continue;
                };
                if val
                    .get("type")
                    .and_then(|v| v.as_str())
                    .is_some_and(|v| v == "metadata")
                {
                    continue;
                }
                if let Some(gene_id) = val.get("gene_id").and_then(|v| v.as_str()) {
                    out.insert(gene_id.to_string());
                    parsed += 1;
                }
            }
            if parsed > 0 {
                log::info!("resume: loaded {} gene ids from qc_report.jsonl", parsed);
            }
        }
        return out;
    }
    let csv_path = Path::new(out_dir).join("qc_summary.csv");
    if csv_path.exists() {
        if let Ok(file) = File::open(&csv_path) {
            let reader = BufReader::new(file);
            let mut header_skipped = false;
            for line in reader.lines() {
                let line = match line {
                    Ok(line) => line,
                    Err(err) => {
                        log::warn!("resume: failed reading qc_summary.csv: {err}");
                        continue;
                    }
                };
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                if !header_skipped {
                    header_skipped = true;
                    continue;
                }
                if let Some(gene_id) = line.split(',').next() {
                    if !gene_id.is_empty() {
                        out.insert(gene_id.to_string());
                    }
                }
            }
            if !out.is_empty() {
                log::info!("resume: loaded {} gene ids from qc_summary.csv", out.len());
            }
        }
    }
    out
}

#[derive(Serialize)]
struct ConfigSnapshot<'a> {
    fasta: &'a str,
    db: &'a str,
    out: &'a str,
    report_format: String,
    threads: usize,
    diamond_bin: &'a str,
    reference_fasta: Option<&'a str>,
    alignment_top_hits: usize,
    alignment_strategy: String,
    coverage_delta_threshold: f64,
    log_format: String,
    scoring_weights: std::collections::HashMap<String, f64>,
    calibration_mode: String,
    calibration_min_samples: usize,
    calibration_min_unique: usize,
}

fn write_run_manifest(
    cfg: &EffectiveConfig,
    tools: &preflight::ToolVersions,
    sums: &Checksums,
    snapshot: &ConfigSnapshot,
    plugins: &[PluginManifest],
    rules: &[RuleManifest],
) -> Result<(), Box<dyn std::error::Error>> {
    #[derive(Serialize)]
    struct Manifest<'a> {
        schema_version: &'a str,
        tool: &'a str,
        diamond_version: &'a str,
        mafft_version: Option<&'a str>,
        hmmscan_version: Option<&'a str>,
        fasta_xx64: Option<&'a str>,
        db_xx64: Option<&'a str>,
        config: &'a ConfigSnapshot<'a>,
        plugins: &'a [PluginManifest],
        rhai_rules: &'a [RuleManifest],
        plugin_schema_version: &'a str,
    }
    let manifest = Manifest {
        schema_version: OUTPUT_SCHEMA_VERSION,
        tool: "AnnoQC",
        diamond_version: tools.diamond.as_deref().unwrap_or_default(),
        mafft_version: tools.mafft.as_deref(),
        hmmscan_version: tools.hmmscan.as_deref(),
        fasta_xx64: sums.fasta_xx64.as_deref(),
        db_xx64: sums.db_xx64.as_deref(),
        config: snapshot,
        plugins,
        rhai_rules: rules,
        plugin_schema_version: PLUGIN_SCHEMA_VERSION,
    };
    let path = Path::new(&cfg.out).join("run.json");
    let text = serde_json::to_string_pretty(&manifest)?;
    fs::write(path, text)?;
    Ok(())
}

fn prepare_cmd(p: PrepareArgs) -> Result<(), Box<dyn std::error::Error>> {
    let diamond_bin = p.diamond_bin.unwrap_or_else(|| "diamond".to_string());
    // Step 1: ensure FASTA exists
    let fasta_path = Path::new(&p.fasta);
    if !fasta_path.exists() {
        eprintln!("prepare: FASTA {} not found. Download it per book/prepare.md or run scripts/fetch_reference_data.sh", fasta_path.display());
        return Ok(());
    }
    // Step 2: makedb (checkpointed)
    let db_out = Path::new(&p.db_out);
    let db_done_buf = format!("{}.done", p.db_out);
    let db_done = Path::new(&db_done_buf);
    let prep_json = matches!(p.log_format, LogFormat::Json);
    enforce_taxonomy_metadata(&diamond_bin, db_out, db_done, "diamond makedb");
    checkpoint::run_step(db_done, p.resume, "diamond makedb", prep_json, || {
        if db_out.exists() {
            log::info!(
                "makedb: {} exists; rebuilding due to resume=false",
                db_out.display()
            );
        }
        // Attempt to build an accession->taxid map from UniProt FASTA (uses OX= in headers)
        let map_path = format!("{}.acc_taxid.tsv", &p.db_out);
        match taxonomy::write_cache_from_fasta(&p.fasta, &map_path) {
            Ok(n) => {
                log::info!("makedb: wrote {} accessions to {}", n, map_path);
                // Convert to DIAMOND's accepted header/format: "accession.version\ttaxid"
                let ncbi_map = format!("{}.ncbi.tsv", &p.db_out);
                if let Ok(text) = std::fs::read_to_string(&map_path) {
                    use std::io::Write as _;
                    if let Ok(mut out) = std::fs::File::create(&ncbi_map) {
                        let _ = writeln!(out, "accession.version\ttaxid");
                        for line in text.lines() {
                            let mut it = line.split('\t');
                            if let (Some(acc), Some(tax)) = (it.next(), it.next()) {
                                let _ = writeln!(out, "{}\t{}", acc, tax);
                            }
                        }
                    }
                }
            }
            Err(e) => log::warn!("makedb: unable to build taxon map from FASTA: {}", e),
        }
        // Try to supply NCBI taxdump if available (downloaded by scripts/fetch_reference_data.sh)
        let taxdump_dir = std::path::Path::new("share/taxonomy/new_taxdump");
        let has_taxdump =
            taxdump_dir.join("nodes.dmp").exists() && taxdump_dir.join("names.dmp").exists();

        let mut cmd = std::process::Command::new(&diamond_bin);
        cmd.arg("makedb")
            .arg("--in")
            .arg(&p.fasta)
            .arg("--db")
            .arg(&p.db_out);
        let ncbi_map = format!("{}.ncbi.tsv", &p.db_out);
        if std::path::Path::new(&ncbi_map).exists() {
            cmd.arg("--taxonmap").arg(&ncbi_map);
        }
        if has_taxdump {
            cmd.arg("--taxonnodes")
                .arg(taxdump_dir.join("nodes.dmp"))
                .arg("--taxonnames")
                .arg(taxdump_dir.join("names.dmp"));
        }
        let status = cmd.status().map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("diamond makedb failed with status {}", status));
        }
        Ok(())
    })?;

    // Step 3: linclust -> clusters
    // Some DIAMOND builds are unstable with gzipped FASTA in linclust/cluster.
    // Decompress to a plain FASTA for these steps to improve stability.
    let fasta_for_cluster = if p.fasta.ends_with(".gz") {
        let uncompressed = Path::new(&p.fasta)
            .with_extension("")
            .to_string_lossy()
            .to_string();
        if !Path::new(&uncompressed).exists() {
            log::info!(
                "prepare: decompressing {} -> {} for linclust/cluster",
                p.fasta,
                uncompressed
            );
            let status = std::process::Command::new("sh")
                .arg("-c")
                .arg(format!("gunzip -c '{}' > '{}'", p.fasta, uncompressed))
                .status()
                .map_err(|e| e.to_string())?;
            if !status.success() {
                return Err(format!("gunzip failed for {} with status {}", p.fasta, status).into());
            }
        }
        uncompressed
    } else {
        p.fasta.clone()
    };
    let clusters = Path::new("clusters");
    let clusters_done = Path::new("clusters.done");
    checkpoint::run_step(
        clusters_done,
        p.resume,
        "diamond linclust",
        prep_json,
        || {
            diamond_linclust(
                &diamond_bin,
                &fasta_for_cluster,
                clusters,
                p.approx_id,
                p.threads,
            )
            .map_err(|e| format!("linclust failed: {}", e))
        },
    )?;

    // Step 4: cluster (sensitive) → clusters.realign
    let realign = Path::new("clusters.realign");
    let realign_done = Path::new("clusters.realign.done");
    checkpoint::run_step(realign_done, p.resume, "diamond cluster", prep_json, || {
        diamond_cluster(
            &diamond_bin,
            &fasta_for_cluster,
            realign,
            p.approx_id,
            p.threads,
        )
        .map_err(|e| format!("cluster failed: {}", e))
    })?;

    // Step 5: recluster → clusters.recluster
    let recluster = Path::new("clusters.recluster");
    let recluster_done = Path::new("clusters.recluster.done");
    checkpoint::run_step(
        recluster_done,
        p.resume,
        "diamond recluster",
        prep_json,
        || {
            diamond_recluster(
                &diamond_bin,
                &fasta_for_cluster,
                realign,
                recluster,
                p.approx_id,
                p.member_cover,
                p.threads,
            )
            .map_err(|e| format!("recluster failed: {}", e))
        },
    )?;
    if recluster.exists() {
        if let Ok(meta) = recluster.metadata() {
            log::info!(
                "prepare: recluster output {} bytes at {}",
                meta.len(),
                recluster.display()
            );
        }
    } else {
        log::warn!(
            "prepare: recluster output missing at {}",
            recluster.display()
        );
    }

    // Step 6: Reference proteomes (Aves) optional build (behind config in future; for now, auto if README exists)
    let refprot_dir = Path::new("share/uniprot/reference_proteomes");
    let aves_done = Path::new("share/refprot/aves/refprot_aves.done");
    let aves_root = Path::new("share/refprot/aves");
    let aves_db = aves_root.join("aves_refprot.dmnd");
    let aves_fasta_gz = aves_root.join("aves_refprot.fasta.gz");
    let aves_proteome_map = aves_root.join("aves_refprot.proteome_map.tsv");
    enforce_taxonomy_metadata(&diamond_bin, &aves_db, aves_done, "refprot_aves");
    ensure_refprot_integrity(&diamond_bin, &aves_db, &aves_fasta_gz, aves_done);
    checkpoint::run_step(aves_done, p.resume, "refprot_aves", prep_json, || {
        // Refresh README if older than 14 days
        let readme_path = refprot_dir.join("README");
        let mut readme_refreshed = false;
        if readme_path.exists() {
            if let Ok(meta) = std::fs::metadata(&readme_path) {
                if let Ok(modified) = meta.modified() {
                    if let Ok(age) = modified.elapsed() {
                        if age.as_secs() > 14 * 24 * 60 * 60 {
                            let url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README";
                            log::info!(
                                "refprot: README older than 14 days; refreshing from {}",
                                url
                            );
                            let dest_path = readme_path.to_string_lossy().to_string();
                            match http_download(url, &dest_path) {
                                Ok(DownloadStatus::Downloaded(_)) => {
                                    readme_refreshed = true;
                                }
                                Ok(DownloadStatus::NotModified) => {}
                                Err(e) => {
                                    log::warn!("refprot: README refresh failed: {}", e);
                                }
                            }
                        }
                    }
                }
            }
        }
        if !readme_path.exists() {
            log::info!("refprot: README not found; skipping Aves build");
            return Ok(());
        }
        std::fs::create_dir_all(aves_root).map_err(|e| e.to_string())?;
        let mut cached_count = 0usize;
        for entry in std::fs::read_dir(aves_root).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or_default();
            if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
                continue;
            }
            cached_count += 1;
        }
        let mut new_count = 0usize;
        let need_download = cached_count == 0 || readme_refreshed;
        if need_download {
            // Build resolver from taxdump for Aves filtering
            let resolver_opt = taxonomy::TaxonomyResolver::from_sources(
                None,
                None,
                Some("share/taxonomy/new_taxdump"),
            )
            .map_err(|e| e.to_string())?;
            let resolver = if let Some(r) = resolver_opt {
                r
            } else {
                return Ok(());
            };
            let entries =
                refprot::parse_readme(readme_path.to_str().unwrap()).map_err(|e| e.to_string())?;
            let aves_taxid = 8782u32; // Aves
            let aves_list =
                refprot::select_by_taxon(&entries, &[aves_taxid], &resolver, usize::MAX);
            let sel = aves_list; // download all available Aves proteomes (full set)
            log::info!(
                "refprot: selected {} Aves proteomes for download",
                sel.len()
            );
            let base = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes";
            // Limited parallel downloads (4 workers)
            let (tx, rx) = std::sync::mpsc::channel::<(std::path::PathBuf, bool)>();
            let jobs = sel
                .into_iter()
                .map(|e| (e.proteome_id, e.division, base.to_string()))
                .collect::<Vec<_>>();
            let mut handles = Vec::new();
            let workers = p.refprot_workers.max(1);
            let chunk_size = jobs.len().div_ceil(workers);
            let max_retries = p.download_retries;
            for chunk in jobs.chunks(chunk_size.max(1)) {
                let chunk = chunk.to_vec();
                let txc = tx.clone();
                let root = aves_root.to_path_buf();
                let retries = max_retries;
                handles.push(std::thread::spawn(move || {
                    for (pid, division, base_url) in chunk {
                        let div = {
                            let mut d = division.clone();
                            if d.is_empty() {
                                d = "eukaryota".into();
                            }
                            let mut ch = d.chars();
                            match ch.next() {
                                Some(c) => format!("{}{}", c.to_ascii_uppercase(), ch.as_str()),
                                None => "Eukaryota".into(),
                            }
                        };
                        let url_dir = format!("{}/{}/{}/", base_url, div, pid);
                        let html = http_get_string(&url_dir).unwrap_or_default();
                        let mut best: Option<String> = None;
                        for tok in html.split(|c: char| c == '"' || c.is_whitespace()) {
                            if tok.starts_with(&pid)
                                && tok.ends_with(".fasta.gz")
                                && !tok.contains("_DNA")
                            {
                                if !tok.contains("additional") {
                                    best = Some(tok.to_string());
                                    break;
                                }
                                if best.is_none() {
                                    best = Some(tok.to_string());
                                }
                            }
                        }
                        if let Some(fname) = best {
                            let dest = root.join(format!("{}.fasta.gz", pid));
                            let dest_string = dest.to_string_lossy().to_string();
                            let file_url = format!("{}{}", url_dir, fname);
                            let mut downloaded_now = false;
                            let mut attempts = 0usize;
                            loop {
                                match http_download(&file_url, &dest_string) {
                                    Ok(DownloadStatus::Downloaded(bytes)) => {
                                        log::info!("refprot: fetched {} ({} bytes)", pid, bytes);
                                        downloaded_now = true;
                                        break;
                                    }
                                    Ok(DownloadStatus::NotModified) => {
                                        log::debug!("refprot: {} already up to date", pid);
                                        break;
                                    }
                                    Err(err) => {
                                        if attempts >= retries {
                                            log::warn!(
                                                "refprot: failed to fetch {} after {} retries: {}",
                                                pid,
                                                retries,
                                                err
                                            );
                                            break;
                                        }
                                        let backoff =
                                            Duration::from_secs(2u64.pow(attempts.min(4) as u32));
                                        log::warn!(
                                            "refprot: retry {} for {} after {}s ({})",
                                            attempts + 1,
                                            pid,
                                            backoff.as_secs(),
                                            err
                                        );
                                        std::thread::sleep(backoff);
                                        attempts += 1;
                                    }
                                }
                            }
                            if dest.exists() {
                                let _ = txc.send((dest.clone(), downloaded_now));
                            }
                        }
                    }
                }));
            }
            drop(tx);
            for h in handles {
                let _ = h.join();
            }
            let mut downloaded: Vec<(std::path::PathBuf, bool)> = Vec::new();
            while let Ok(p) = rx.recv_timeout(std::time::Duration::from_millis(10)) {
                downloaded.push(p);
            }
            new_count = downloaded.iter().filter(|(_, fresh)| *fresh).count();
        } else {
            log::info!(
                "refprot: cached {} proteomes present; skipping download",
                cached_count
            );
        }
        let mut concat_invalid = false;
        if aves_fasta_gz.exists() {
            if let Err(e) = validate_gzip_file(&aves_fasta_gz) {
                log::warn!(
                    "refprot: combined FASTA {} failed validation ({}); rebuilding",
                    aves_fasta_gz.display(),
                    e
                );
                let _ = fs::remove_file(&aves_fasta_gz);
                concat_invalid = true;
            }
        }
        let need_concat = concat_invalid || !aves_fasta_gz.exists() || new_count > 0;
        if need_concat {
            log::info!(
                "refprot: rebuilding combined FASTA ({} new proteomes)",
                new_count
            );
            let mut files: Vec<PathBuf> = Vec::new();
            for entry in std::fs::read_dir(aves_root).map_err(|e| e.to_string())? {
                let path = entry.map_err(|e| e.to_string())?.path();
                if !path.is_file() {
                    continue;
                }
                let name = path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or_default();
                if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
                    continue;
                }
                files.push(path);
            }
            files.sort();
            if files.is_empty() {
                log::warn!(
                    "refprot: no proteome FASTAs available under {}",
                    aves_root.display()
                );
            } else {
                rebuild_combined_gzip(&files, &aves_fasta_gz)?;
            }
        }

        // Rebuild proteome map if we grabbed new files or map is missing
        if !aves_proteome_map.exists() || new_count > 0 {
            build_refprot_proteome_map(aves_root, &aves_proteome_map)
                .map_err(|e| format!("refprot: proteome map build failed: {}", e))?;
        }

        let missing_taxonomy = aves_db.exists() && !diamond_db_has_taxonomy(&diamond_bin, &aves_db);
        let db_failed_validation =
            aves_db.exists() && !diamond_db_integrity_ok(&diamond_bin, &aves_db);
        // Decide if we need to (re)build makedb: when DB is missing, lacks taxonomy, or FASTAs changed
        let need_makedb = missing_taxonomy
            || db_failed_validation
            || !aves_db.exists()
            || new_count > 0
            || !aves_fasta_gz.exists();
        if need_makedb {
            // Build accession->taxid cache from FASTA headers (OX=) and feed to makedb
            let aves_acc = aves_root.join("aves_refprot.acc_taxid.tsv");
            let aves_ncbi = aves_root.join("aves_refprot.ncbi.tsv");
            match taxonomy::write_cache_from_fasta(
                &aves_fasta_gz.to_string_lossy(),
                &aves_acc.to_string_lossy(),
            ) {
                Ok(n) => {
                    log::info!("refprot: wrote {} accessions to {}", n, aves_acc.display());
                    if let Ok(text) = std::fs::read_to_string(&aves_acc) {
                        if let Ok(mut f) = std::fs::File::create(&aves_ncbi) {
                            use std::io::Write as _;
                            let _ = writeln!(f, "accession.version\ttaxid");
                            for line in text.lines() {
                                let mut it = line.split('\t');
                                if let (Some(acc), Some(tid)) = (it.next(), it.next()) {
                                    let _ = writeln!(f, "{}\t{}", acc, tid);
                                }
                            }
                        }
                    }
                }
                Err(e) => log::warn!("refprot: unable to build taxon map: {}", e),
            }
            // Build DIAMOND DB with taxonomy (if available)
            let taxdump_dir = std::path::Path::new("share/taxonomy/new_taxdump");
            let has_taxdump =
                taxdump_dir.join("nodes.dmp").exists() && taxdump_dir.join("names.dmp").exists();
            let tmp_db = aves_root.join("aves_refprot.tmp.dmnd");
            if tmp_db.exists() {
                let _ = fs::remove_file(&tmp_db);
            }
            let mut cmd = std::process::Command::new(&diamond_bin);
            cmd.arg("makedb")
                .arg("--in")
                .arg(&aves_fasta_gz)
                .arg("--db")
                .arg(&tmp_db);
            if aves_ncbi.exists() {
                cmd.arg("--taxonmap").arg(&aves_ncbi);
            }
            if has_taxdump {
                cmd.arg("--taxonnodes")
                    .arg(taxdump_dir.join("nodes.dmp"))
                    .arg("--taxonnames")
                    .arg(taxdump_dir.join("names.dmp"));
            }
            let status = cmd.status().map_err(|e| e.to_string())?;
            if !status.success() {
                let _ = fs::remove_file(&tmp_db);
                return Err(format!(
                    "diamond makedb (refprot) failed with status {}",
                    status
                ));
            }
            let info = run_diamond_dbinfo(&diamond_bin, &tmp_db)
                .map_err(|e| format!("refprot: dbinfo failed after makedb: {}", e))?;
            let hash = dbinfo_extract_hash(&info);
            if aves_db.exists() {
                let _ = fs::remove_file(&aves_db);
            }
            fs::rename(&tmp_db, &aves_db)
                .map_err(|e| format!("refprot: rename tmp db failed: {}", e))?;
            if let Some(hash) = hash {
                let hash_path = aves_root.join("aves_refprot.hash");
                let _ = fs::write(hash_path, format!("{}\n", hash));
            }
        }
        Ok(())
    })?;
    Ok(())
}

fn build_refprot_proteome_map(
    root: &Path,
    map_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut seen: HashSet<String> = HashSet::new();
    let mut out = File::create(map_path)?;
    for entry in std::fs::read_dir(root)? {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
            continue;
        }
        let proteome_id = name.trim_end_matches(".fasta.gz").to_string();
        let mut reader = parse_fastx_file(&path)
            .map_err(|e| format!("refprot map parse {}: {}", path.display(), e))?;
        while let Some(rec) = reader.next() {
            let rec = rec.map_err(|e| format!("refprot map parse {}: {}", path.display(), e))?;
            let acc = String::from_utf8_lossy(rec.id()).to_string();
            if seen.insert(acc.clone()) {
                writeln!(out, "{}\t{}", acc, proteome_id)?;
            }
        }
    }
    Ok(())
}

fn rebuild_combined_gzip(files: &[PathBuf], out_path: &Path) -> Result<(), String> {
    use std::io::copy;
    let out_file = File::create(out_path).map_err(|e| e.to_string())?;
    let mut encoder = GzEncoder::new(out_file, Compression::default());
    for path in files {
        let file = File::open(path).map_err(|e| e.to_string())?;
        let mut decoder = MultiGzDecoder::new(file);
        copy(&mut decoder, &mut encoder)
            .map_err(|e| format!("refprot: concat {} failed: {}", path.display(), e))?;
    }
    encoder.finish().map_err(|e| e.to_string())?;
    Ok(())
}

fn validate_gzip_file(path: &Path) -> Result<(), String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut sink = io::sink();
    io::copy(&mut decoder, &mut sink).map_err(|e| e.to_string())?;
    Ok(())
}

fn diamond_db_integrity_ok(diamond_bin: &str, db_path: &Path) -> bool {
    match run_diamond_dbinfo(diamond_bin, db_path) {
        Ok(_) => true,
        Err(e) => {
            log::warn!("diamond dbinfo failed for {}: {}", db_path.display(), e);
            false
        }
    }
}

fn run_diamond_dbinfo(diamond_bin: &str, db_path: &Path) -> Result<String, String> {
    if !db_path.exists() {
        return Err("db missing".into());
    }
    let output = std::process::Command::new(diamond_bin)
        .arg("dbinfo")
        .arg("--db")
        .arg(db_path)
        .output()
        .map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!("status {}", output.status));
    }
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn dbinfo_extract_hash(text: &str) -> Option<String> {
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("Database hash") {
            return trimmed.split_whitespace().last().map(|s| s.to_string());
        }
    }
    None
}

fn ensure_refprot_integrity(
    diamond_bin: &str,
    db_path: &Path,
    fasta_path: &Path,
    done_marker: &Path,
) {
    if !done_marker.exists() {
        return;
    }
    let mut needs_rerun = false;
    if fasta_path.exists() {
        if let Err(e) = validate_gzip_file(fasta_path) {
            log::warn!(
                "refprot: detected corrupt combined FASTA {} ({}); scheduling rebuild",
                fasta_path.display(),
                e
            );
            let _ = fs::remove_file(fasta_path);
            needs_rerun = true;
        }
    }
    if db_path.exists() && !diamond_db_integrity_ok(diamond_bin, db_path) {
        log::warn!(
            "refprot: DIAMOND database {} failed validation; scheduling rebuild",
            db_path.display()
        );
        let _ = fs::remove_file(db_path);
        needs_rerun = true;
    }
    if needs_rerun {
        log::warn!(
            "refprot: marking {} stale so step will be rerun",
            done_marker.display()
        );
        let _ = fs::remove_file(done_marker);
    }
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

fn load_refprot_proteome_map(path: &Path) -> HashMap<String, String> {
    let mut map = HashMap::new();
    if let Ok(text) = fs::read_to_string(path) {
        for line in text.lines() {
            let mut parts = line.split('\t');
            if let (Some(acc), Some(pid)) = (parts.next(), parts.next()) {
                map.insert(acc.to_string(), pid.to_string());
            }
        }
    }
    map
}

fn annotate_refprot_hits(
    grouped: &mut HashMap<String, Vec<diamond::DiamondHitRow>>,
    map: &HashMap<String, String>,
) {
    for hits in grouped.values_mut() {
        for hit in hits.iter_mut() {
            let proteome_id = map
                .get(&hit.sseqid)
                .cloned()
                .unwrap_or_else(|| "refprot".to_string());
            hit.source = HitSource::RefProt(proteome_id);
        }
    }
}

fn parse_evalue_to_f64(value: &str) -> f64 {
    if value.trim().is_empty() {
        return 1.0;
    }
    value.trim().parse::<f64>().unwrap_or_else(|_| {
        match value.trim().to_ascii_lowercase().as_str() {
            "inf" => f64::INFINITY,
            _ => 1.0,
        }
    })
}

fn filter_refprot_hits(
    grouped: &mut HashMap<String, Vec<diamond::DiamondHitRow>>,
    cfg: &RefProtFallbackConfig,
) {
    let mut empty_keys = Vec::new();
    for (gene, hits) in grouped.iter_mut() {
        hits.retain(|row| {
            if row.qcov < cfg.min_qcov {
                return false;
            }
            if row.scov < cfg.min_scov {
                return false;
            }
            if row.pident < cfg.min_pident {
                return false;
            }
            let eval = parse_evalue_to_f64(&row.evalue);
            if eval.is_nan() || eval > cfg.max_evalue {
                return false;
            }
            true
        });
        hits.sort_by(|a, b| b.bitscore.total_cmp(&a.bitscore));
        if hits.len() > cfg.max_hits {
            hits.truncate(cfg.max_hits);
        }
        if hits.is_empty() {
            empty_keys.push(gene.clone());
        }
    }
    for key in empty_keys {
        grouped.remove(&key);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecs::GeneMetrics;
    use crate::taxonomy::TaxonomyEvidence;
    use tempfile::TempDir;

    fn build_ref_hit(
        gene: &str,
        subj: &str,
        qcov: f64,
        scov: f64,
        pident: f64,
        evalue: &str,
        bitscore: f64,
    ) -> diamond::DiamondHitRow {
        diamond::DiamondHitRow {
            qseqid: gene.into(),
            sseqid: subj.into(),
            bitscore,
            evalue: evalue.into(),
            length: 120,
            qcov,
            scov,
            pident,
            qstart: 1,
            qend: 100,
            sstart: 5,
            send: 105,
            qlen: 100,
            slen: 110,
            source: HitSource::RefProt("P0001".into()),
            staxid: Some(1),
            lineage: vec![],
        }
    }

    #[test]
    fn filter_refprot_hits_respects_thresholds_and_cap() {
        let mut grouped: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
        grouped.insert(
            "gene1".into(),
            vec![
                build_ref_hit("gene1", "bad_qcov", 0.2, 0.4, 40.0, "1e-20", 210.0),
                build_ref_hit("gene1", "good1", 0.8, 0.6, 45.0, "1e-20", 260.0),
                build_ref_hit("gene1", "good2", 0.9, 0.7, 50.0, "1e-15", 240.0),
                build_ref_hit("gene1", "high_eval", 0.9, 0.7, 50.0, "1e-2", 230.0),
            ],
        );
        grouped.insert(
            "gene2".into(),
            vec![build_ref_hit(
                "gene2", "fail_all", 0.1, 0.1, 10.0, "1", 200.0,
            )],
        );
        let cfg = RefProtFallbackConfig {
            min_qcov: 0.5,
            min_scov: 0.4,
            min_pident: 30.0,
            max_evalue: 1e-10,
            max_hits: 1,
            ..Default::default()
        };
        let mut cfg = cfg;
        cfg.trigger_k = 5;
        filter_refprot_hits(&mut grouped, &cfg);
        let kept = grouped.get("gene1").unwrap();
        assert_eq!(kept.len(), 1, "max_hits should truncate to 1 record");
        assert_eq!(kept[0].sseqid, "good1");
        assert!(
            !grouped.contains_key("gene2"),
            "genes with no survivors removed"
        );
    }

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
        let res = export_high_sequences(out_dir, None, &metrics, &seqs, &scores)
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
        assert!(dbinfo_text_has_taxonomy(ok));
        let bad = "Database sequences\nNo taxonomy found";
        assert!(!dbinfo_text_has_taxonomy(bad));
    }
}

fn http_get_string(url: &str) -> Result<String, String> {
    use curl::easy::Easy;
    let mut data = Vec::new();
    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    let mut transfer = easy.transfer();
    transfer
        .write_function(|new| {
            data.extend_from_slice(new);
            Ok(new.len())
        })
        .map_err(|e| e.to_string())?;
    transfer.perform().map_err(|e| e.to_string())?;
    drop(transfer);
    Ok(String::from_utf8_lossy(&data).to_string())
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct HttpCacheMeta {
    etag: Option<String>,
    last_modified: Option<String>,
}

impl HttpCacheMeta {
    fn is_empty(&self) -> bool {
        self.etag.is_none() && self.last_modified.is_none()
    }
}

enum DownloadStatus {
    NotModified,
    Downloaded(u64),
}

fn load_http_cache_meta(path: &Path) -> Option<HttpCacheMeta> {
    fs::read_to_string(path)
        .ok()
        .and_then(|txt| serde_json::from_str(&txt).ok())
}

fn save_http_cache_meta(path: &Path, meta: &HttpCacheMeta) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let text = serde_json::to_string(meta).map_err(|e| e.to_string())?;
    fs::write(path, text).map_err(|e| e.to_string())
}

fn remote_not_modified(url: &str, meta: &HttpCacheMeta) -> Result<bool, String> {
    use curl::easy::{Easy, List};
    if meta.is_empty() {
        return Ok(false);
    }
    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    easy.nobody(true).map_err(|e| e.to_string())?;
    let mut headers = List::new();
    let mut has = false;
    if let Some(etag) = &meta.etag {
        headers
            .append(&format!("If-None-Match: {}", etag))
            .map_err(|e| e.to_string())?;
        has = true;
    }
    if let Some(lm) = &meta.last_modified {
        headers
            .append(&format!("If-Modified-Since: {}", lm))
            .map_err(|e| e.to_string())?;
        has = true;
    }
    if !has {
        return Ok(false);
    }
    easy.http_headers(headers).map_err(|e| e.to_string())?;
    easy.perform().map_err(|e| e.to_string())?;
    let code = easy.response_code().map_err(|e| e.to_string())?;
    Ok(code == 304)
}

fn http_download(url: &str, dest: &str) -> Result<DownloadStatus, String> {
    use curl::easy::{Easy, List, WriteError};
    use std::io::Write as _;
    let dest_path = Path::new(dest);
    let tmp_path = PathBuf::from(format!("{}.part", dest));
    let meta_path = PathBuf::from(format!("{}.httpmeta", dest));
    let mut meta = load_http_cache_meta(&meta_path).unwrap_or_default();

    if dest_path.exists() && !tmp_path.exists() && remote_not_modified(url, &meta)? {
        return Ok(DownloadStatus::NotModified);
    }

    if tmp_path.exists() && !meta.is_empty() && !remote_not_modified(url, &meta)? {
        let _ = fs::remove_file(&tmp_path);
        meta = HttpCacheMeta::default();
    }

    let resume_from = if tmp_path.exists() {
        fs::metadata(&tmp_path).map(|m| m.len()).unwrap_or(0)
    } else {
        0
    };
    if resume_from == 0 {
        if let Some(parent) = tmp_path.parent() {
            fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        if tmp_path.exists() {
            let _ = fs::remove_file(&tmp_path);
        }
    }
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&tmp_path)
        .map_err(|e| e.to_string())?;

    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    if resume_from > 0 {
        easy.resume_from(resume_from).map_err(|e| e.to_string())?;
    }
    let mut headers = List::new();
    let mut has_headers = false;
    if resume_from == 0 {
        if let Some(etag) = &meta.etag {
            headers
                .append(&format!("If-None-Match: {}", etag))
                .map_err(|e| e.to_string())?;
            has_headers = true;
        }
        if let Some(lm) = &meta.last_modified {
            headers
                .append(&format!("If-Modified-Since: {}", lm))
                .map_err(|e| e.to_string())?;
            has_headers = true;
        }
    } else if let Some(etag) = &meta.etag {
        headers
            .append(&format!("If-Range: {}", etag))
            .map_err(|e| e.to_string())?;
        has_headers = true;
    } else if let Some(lm) = &meta.last_modified {
        headers
            .append(&format!("If-Range: {}", lm))
            .map_err(|e| e.to_string())?;
        has_headers = true;
    }
    if has_headers {
        easy.http_headers(headers).map_err(|e| e.to_string())?;
    }
    let header_meta = Arc::new(Mutex::new(HttpCacheMeta::default()));
    let header_clone = Arc::clone(&header_meta);
    easy.header_function(move |header| {
        if let Ok(text) = std::str::from_utf8(header) {
            let lower = text.to_ascii_lowercase();
            if lower.starts_with("etag:") {
                if let Ok(mut guard) = header_clone.lock() {
                    guard.etag = Some(
                        text.split_once(':')
                            .map(|(_, v)| v.trim().trim_matches('"').to_string())
                            .unwrap_or_default(),
                    );
                }
            } else if lower.starts_with("last-modified:") {
                if let Ok(mut guard) = header_clone.lock() {
                    guard.last_modified = Some(
                        text.split_once(':')
                            .map(|(_, v)| v.trim().to_string())
                            .unwrap_or_default(),
                    );
                }
            }
        }
        true
    })
    .map_err(|e| e.to_string())?;

    let mut written: u64 = 0;
    {
        let mut transfer = easy.transfer();
        transfer
            .write_function(|new| {
                file.write_all(new).map_err(|_| WriteError::Pause)?;
                written += new.len() as u64;
                Ok(new.len())
            })
            .map_err(|e| e.to_string())?;
        if let Err(e) = transfer.perform() {
            let _ = fs::remove_file(&tmp_path);
            return Err(e.to_string());
        }
    }
    let code = easy.response_code().map_err(|e| e.to_string())?;
    if code == 304 {
        let _ = fs::remove_file(&tmp_path);
        return Ok(DownloadStatus::NotModified);
    }
    if code == 416 && dest_path.exists() {
        let _ = fs::remove_file(&tmp_path);
        return Ok(DownloadStatus::NotModified);
    }
    if !(200..300).contains(&code) {
        let _ = fs::remove_file(&tmp_path);
        return Err(format!("download failed with status {}", code));
    }
    file.flush().map_err(|e| e.to_string())?;
    drop(file);
    if dest_path.exists() {
        fs::remove_file(dest_path).map_err(|e| e.to_string())?;
    }
    fs::rename(&tmp_path, dest_path).map_err(|e| e.to_string())?;
    drop(easy);
    let final_meta = Arc::try_unwrap(header_meta)
        .ok()
        .and_then(|m| m.into_inner().ok())
        .unwrap_or_default();
    if !final_meta.is_empty() {
        meta = final_meta;
    }
    save_http_cache_meta(&meta_path, &meta)?;
    let final_size = fs::metadata(dest_path).map(|m| m.len()).unwrap_or(0);
    if final_size == 0 {
        return Err("download produced empty file".into());
    }
    Ok(DownloadStatus::Downloaded(written))
}

fn diamond_db_has_taxonomy(diamond_bin: &str, db: &Path) -> bool {
    if !db.exists() {
        return false;
    }
    if let Ok(text) = run_diamond_dbinfo(diamond_bin, db) {
        if dbinfo_text_has_taxonomy(&text) {
            return true;
        }
    }
    let acc_taxid_candidates = [
        std::path::PathBuf::from(format!("{}.acc_taxid.tsv", db.display())),
        db.with_extension("acc_taxid.tsv"),
    ];
    for acc_taxid in acc_taxid_candidates {
        if let Ok(meta) = fs::metadata(&acc_taxid) {
            if meta.len() > 0 {
                return true;
            }
        }
    }
    false
}

fn dbinfo_text_has_taxonomy(text: &str) -> bool {
    let lower = text.to_lowercase();
    lower.contains("taxon") && !lower.contains("no taxonomy")
}

fn enforce_taxonomy_metadata(diamond_bin: &str, db_path: &Path, done_marker: &Path, label: &str) {
    if db_path.exists() && done_marker.exists() && !diamond_db_has_taxonomy(diamond_bin, db_path) {
        log::warn!(
            "{} missing taxonomy metadata; forcing rebuild",
            db_path.display()
        );
        let _ = fs::remove_file(done_marker);
        log::warn!("{} step will rerun to add taxonomy", label);
    }
}

fn scoring_thresholds(scoring: &Option<ScoringConfigOverride>) -> (f64, f64) {
    if let Some(cfg) = scoring {
        let h = cfg.thresholds.as_ref().and_then(|t| t.high).unwrap_or(0.8);
        let m = cfg
            .thresholds
            .as_ref()
            .and_then(|t| t.medium)
            .unwrap_or(0.5);
        (h, m)
    } else {
        (0.8, 0.5)
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

fn format_classification(base: &str, missing: &[&'static str], classify_no_data: bool) -> String {
    if !classify_no_data || missing.is_empty() {
        return base.to_string();
    }
    let list = missing.join("|");
    format!("{} (no_data:{})", base, list)
}

fn classification_base(label: &str) -> &str {
    label.split([' ', '(', '[']).next().unwrap_or(label)
}

fn print_scoring_rubric(
    scoring: &Option<ScoringConfigOverride>,
    calibration: CalibrationSettings,
    classify_no_data: bool,
) {
    let weights = scoring_weights(scoring);
    let (th_high, th_med) = scoring_thresholds(scoring);
    let sum = weights.sum().max(1e-9);
    let norm = |v: f64| v / sum;
    println!("scoring_rubric");
    println!(
        "weights_raw: homology={:.3} intrinsic={:.3} taxonomy={:.3} domains={:.3} domains_strength={:.3} length={:.3} orphan={:.3} subject_cov={:.3} termini={:.3} divergence={:.3} conserved_regions={:.3} genomic={:.3} sum={:.3}",
        weights.homology,
        weights.intrinsic,
        weights.taxonomy,
        weights.domains,
        weights.domains_strength,
        weights.length,
        weights.orphan,
        weights.subject_cov,
        weights.termini,
        weights.divergence,
        weights.conserved_regions,
        weights.genomic,
        weights.sum(),
    );
    println!(
        "weights_normalized: homology={:.3} intrinsic={:.3} taxonomy={:.3} domains={:.3} domains_strength={:.3} length={:.3} orphan={:.3} subject_cov={:.3} termini={:.3} divergence={:.3} conserved_regions={:.3} genomic={:.3}",
        norm(weights.homology),
        norm(weights.intrinsic),
        norm(weights.taxonomy),
        norm(weights.domains),
        norm(weights.domains_strength),
        norm(weights.length),
        norm(weights.orphan),
        norm(weights.subject_cov),
        norm(weights.termini),
        norm(weights.divergence),
        norm(weights.conserved_regions),
        norm(weights.genomic),
    );
    println!("thresholds: high={:.3} medium={:.3}", th_high, th_med);
    println!(
        "calibration: mode={:?} min_samples={} min_unique={}",
        calibration.mode, calibration.min_samples, calibration.min_unique
    );
    println!(
        "no_data_handling: classify_no_data={} missing_pillars_excluded_from_weights=true",
        classify_no_data
    );
}

#[derive(Clone, Copy, Debug)]
struct WeightSet {
    homology: f64,
    intrinsic: f64,
    taxonomy: f64,
    domains: f64,
    domains_strength: f64,
    length: f64,
    orphan: f64,
    subject_cov: f64,
    termini: f64,
    divergence: f64,
    conserved_regions: f64,
    genomic: f64,
}

impl WeightSet {
    fn sum(&self) -> f64 {
        self.homology
            + self.intrinsic
            + self.taxonomy
            + self.domains
            + self.domains_strength
            + self.length
            + self.orphan
            + self.subject_cov
            + self.termini
            + self.divergence
            + self.conserved_regions
            + self.genomic
    }

    fn sum_available(&self, presence: &PillarPresence) -> f64 {
        let mut sum = 0.0;
        if presence.homology {
            sum += self.homology;
        }
        if presence.intrinsic {
            sum += self.intrinsic;
        }
        if presence.taxonomy {
            sum += self.taxonomy;
        }
        if presence.domains {
            sum += self.domains;
        }
        if presence.domains_strength {
            sum += self.domains_strength;
        }
        if presence.length {
            sum += self.length;
        }
        if presence.orphan {
            sum += self.orphan;
        }
        if presence.subject_cov {
            sum += self.subject_cov;
        }
        if presence.termini {
            sum += self.termini;
        }
        if presence.divergence {
            sum += self.divergence;
        }
        if presence.conserved_regions {
            sum += self.conserved_regions;
        }
        if presence.genomic {
            sum += self.genomic;
        }
        sum
    }
}

#[derive(Clone, Debug, Default)]
struct PillarPresence {
    homology: bool,
    intrinsic: bool,
    taxonomy: bool,
    domains: bool,
    domains_strength: bool,
    length: bool,
    orphan: bool,
    subject_cov: bool,
    termini: bool,
    divergence: bool,
    conserved_regions: bool,
    genomic: bool,
}

impl PillarPresence {
    fn missing_with_weights(&self, weights: &WeightSet) -> Vec<&'static str> {
        let mut missing = Vec::new();
        if weights.homology > 0.0 && !self.homology {
            missing.push("homology");
        }
        if weights.intrinsic > 0.0 && !self.intrinsic {
            missing.push("intrinsic");
        }
        if weights.taxonomy > 0.0 && !self.taxonomy {
            missing.push("taxonomy");
        }
        if weights.domains > 0.0 && !self.domains {
            missing.push("domains");
        }
        if weights.domains_strength > 0.0 && !self.domains_strength {
            missing.push("domains_strength");
        }
        if weights.length > 0.0 && !self.length {
            missing.push("length");
        }
        if weights.orphan > 0.0 && !self.orphan {
            missing.push("orphan");
        }
        if weights.subject_cov > 0.0 && !self.subject_cov {
            missing.push("subject_cov");
        }
        if weights.termini > 0.0 && !self.termini {
            missing.push("termini");
        }
        if weights.divergence > 0.0 && !self.divergence {
            missing.push("divergence");
        }
        if weights.conserved_regions > 0.0 && !self.conserved_regions {
            missing.push("conserved_regions");
        }
        if weights.genomic > 0.0 && !self.genomic {
            missing.push("genomic");
        }
        missing
    }
}

fn scoring_weights(scoring: &Option<ScoringConfigOverride>) -> WeightSet {
    let mut ws = WeightSet {
        homology: 0.55,
        intrinsic: 0.3,
        taxonomy: 0.0,
        domains: 0.0,
        domains_strength: 0.05,
        length: 0.0,
        orphan: 0.0,
        subject_cov: 0.0,
        termini: 0.0,
        divergence: 0.0,
        conserved_regions: 0.1,
        genomic: 0.0,
    };
    if let Some(cfg) = scoring {
        if let Some(v) = cfg.weights.get("homology") {
            ws.homology = *v;
        }
        if let Some(v) = cfg.weights.get("intrinsic") {
            ws.intrinsic = *v;
        }
        if let Some(v) = cfg.weights.get("taxonomy") {
            ws.taxonomy = *v;
        }
        if let Some(v) = cfg.weights.get("domains") {
            ws.domains = *v;
        }
        if let Some(v) = cfg.weights.get("domains_strength") {
            ws.domains_strength = *v;
        }
        if let Some(v) = cfg.weights.get("length") {
            ws.length = *v;
        }
        if let Some(v) = cfg.weights.get("orphan") {
            ws.orphan = *v;
        }
        if let Some(v) = cfg.weights.get("subject_cov") {
            ws.subject_cov = *v;
        }
        if let Some(v) = cfg.weights.get("termini") {
            ws.termini = *v;
        }
        if let Some(v) = cfg.weights.get("divergence") {
            ws.divergence = *v;
        }
        if let Some(v) = cfg.weights.get("conserved_regions") {
            ws.conserved_regions = *v;
        }
        if let Some(v) = cfg.weights.get("genomic") {
            ws.genomic = *v;
        }
    }
    ws
}

#[derive(Clone, Debug, Default)]
pub(crate) struct ComponentScores {
    homology: f64,
    intrinsic: f64,
    taxonomy: Option<f64>,
    orphan: f64,
    subject_cov: f64,
    termini: Option<f64>,
    divergence: Option<f64>,
    conserved_regions: Option<f64>,
    genomic: f64,
}

#[derive(Clone, Debug)]
struct ScoreTemp {
    gene_id: String,
    raw_score: f64,
    final_score: f64,
    available_weight: f64,
    missing: Vec<&'static str>,
}

fn apply_percentile_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let len = order.len();
    for (rank, idx) in order.into_iter().enumerate() {
        let percentile = if len > 1 {
            (rank as f64 + 0.5) / len as f64
        } else {
            1.0
        };
        entries[idx].final_score = percentile.clamp(0.0, 1.0);
    }
}

fn apply_isotonic_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let n = order.len();
    #[derive(Clone)]
    struct Block {
        start: usize,
        end: usize,
        sum: f64,
        weight: usize,
    }
    let mut blocks: Vec<Block> = Vec::new();
    for (rank, _idx) in order.iter().enumerate() {
        let y = if n > 1 {
            (rank as f64 + 0.5) / n as f64
        } else {
            1.0
        };
        blocks.push(Block {
            start: rank,
            end: rank,
            sum: y,
            weight: 1,
        });
        while blocks.len() >= 2 {
            let k = blocks.len() - 1;
            let prev = &blocks[k - 1];
            let curr = &blocks[k];
            let avg_prev = prev.sum / prev.weight as f64;
            let avg_curr = curr.sum / curr.weight as f64;
            if avg_prev <= avg_curr {
                break;
            }
            let merged = Block {
                start: prev.start,
                end: curr.end,
                sum: prev.sum + curr.sum,
                weight: prev.weight + curr.weight,
            };
            blocks.pop();
            blocks.pop();
            blocks.push(merged);
        }
    }
    let mut fitted = vec![0.0; n];
    for block in blocks {
        let avg = (block.sum / block.weight as f64).clamp(0.0, 1.0);
        for value in fitted
            .iter_mut()
            .take(block.end.saturating_add(1))
            .skip(block.start)
        {
            *value = avg;
        }
    }
    for (rank, idx) in order.into_iter().enumerate() {
        entries[idx].final_score = fitted[rank];
    }
}

fn calibration_has_min_samples(
    entries: &[ScoreTemp],
    min_samples: usize,
    min_unique: usize,
) -> bool {
    if entries.len() < min_samples {
        return false;
    }
    let mut unique: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for e in entries {
        unique.insert(e.raw_score.to_bits());
        if unique.len() >= min_unique {
            return true;
        }
    }
    false
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

fn build_config_snapshot<'a>(
    cfg: &'a EffectiveConfig,
    args: &'a AnalyzeArgs,
    file_cfg: &'a FileConfig,
    calibration: &CalibrationSettings,
    report_format: ReportFormat,
) -> ConfigSnapshot<'a> {
    let mut weights = HashMap::new();
    if let Some(sc) = &file_cfg.scoring {
        weights = sc.weights.clone();
    } else {
        // Keep in sync with `scoring_weights` defaults.
        weights.insert("homology".to_string(), 0.55);
        weights.insert("intrinsic".to_string(), 0.3);
        weights.insert("taxonomy".to_string(), 0.0);
        weights.insert("domains".to_string(), 0.0);
        weights.insert("domains_strength".to_string(), 0.05);
        weights.insert("length".to_string(), 0.0);
        weights.insert("orphan".to_string(), 0.0);
        weights.insert("subject_cov".to_string(), 0.0);
        weights.insert("termini".to_string(), 0.0);
        weights.insert("divergence".to_string(), 0.0);
        weights.insert("conserved_regions".to_string(), 0.1);
        weights.insert("genomic".to_string(), 0.0);
    }
    ConfigSnapshot {
        fasta: &cfg.fasta,
        db: &cfg.db,
        out: &cfg.out,
        report_format: format!("{:?}", report_format),
        threads: cfg.threads,
        diamond_bin: &cfg.diamond_bin,
        reference_fasta: cfg.reference_fasta.as_deref(),
        alignment_top_hits: args.alignment_top_hits,
        alignment_strategy: format!("{:?}", args.alignment_strategy),
        coverage_delta_threshold: args.coverage_delta_threshold,
        log_format: format!("{:?}", args.log_format),
        scoring_weights: weights,
        calibration_mode: format!("{:?}", calibration.mode),
        calibration_min_samples: calibration.min_samples,
        calibration_min_unique: calibration.min_unique,
    }
}

fn adjust_homology_score(
    raw: f64,
    divergence_score: f64,
    divergence_present: bool,
    length_score: f64,
    length_present: bool,
) -> f64 {
    let mut score = raw;
    if divergence_present {
        let adj = (0.5 + 0.5 * divergence_score.clamp(0.0, 1.0)).clamp(0.0, 1.0);
        score *= adj;
    }
    if length_present {
        let adj = (0.5 + 0.5 * length_score.clamp(0.0, 1.0)).clamp(0.0, 1.0);
        score *= adj;
    }
    score.clamp(0.0, 1.0)
}

#[allow(clippy::too_many_arguments)]
fn build_scores_map(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
    hmmer_enabled: bool,
    hmmsum_map: Option<&HashMap<String, hmmer::HmmscanSummary>>,
    arch_map: Option<&HashMap<String, f64>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
    alignment_map: Option<&HashMap<String, mafft::AlignmentMetrics>>,
    structvar_map: Option<&HashMap<String, structvar::StructVar>>,
    genomic_map: Option<&HashMap<String, genomic::GenomicMetrics>>,
    calibration: CalibrationSettings,
    classify_no_data: bool,
) -> (HashMap<String, (f64, String)>, HashMap<String, f64>) {
    let mut weights = scoring_weights(scoring);
    if !hmmer_enabled {
        if weights.domains > 0.0 || weights.domains_strength > 0.0 || weights.orphan > 0.0 {
            log::warn!("hmmer disabled; ignoring domains/domains_strength/orphan weights");
        }
        weights.domains = 0.0;
        weights.domains_strength = 0.0;
        weights.orphan = 0.0;
    }
    let (th_high, th_med) = scoring_thresholds(scoring);
    let genomic_cfg = scoring.as_ref().and_then(|s| s.genomic.as_ref());
    let caps_cfg = scoring.as_ref().and_then(|s| s.caps.as_ref());
    let mut staging: Vec<ScoreTemp> = Vec::with_capacity(metrics.len());
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let homology_present = s.map(|s| s.count > 0).unwrap_or(false);
        let intrinsic_present = intrinsic.contains_key(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h_raw = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let (t, taxonomy_present) = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            let present = evidence
                .map(|ev| ev.considered > 0 || ev.top_hit.is_some())
                .unwrap_or(false);
            (compute_taxonomy_score(evidence), present)
        } else {
            (0.0, false)
        };
        let d = arch_map
            .and_then(|am| am.get(&m.gene_id))
            .cloned()
            .unwrap_or(0.0);
        let domains_present = arch_map
            .map(|am| am.contains_key(&m.gene_id))
            .unwrap_or(false);
        let domains_strength_score =
            compute_domains_strength_score(hmmsum_map.and_then(|hm| hm.get(&m.gene_id)));
        let domains_strength_present = hmmsum_map
            .map(|hm| hm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let l = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let length_present = len_map
            .map(|lm| lm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let subject_cov_score = compute_subject_cov_score(s);
        let subject_cov_present = homology_present;
        let o = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        let orphan_present = orphan_map
            .map(|om| om.contains_key(&m.gene_id))
            .unwrap_or(false);
        let align_entry = alignment_map.and_then(|am| am.get(&m.gene_id));
        let termini_present = align_entry
            .map(|a| a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI)
            .unwrap_or(false);
        let termini_score = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            })
            .unwrap_or(0.0);
        let divergence_present = align_entry
            .map(|a| a.sequences_aligned > 0)
            .unwrap_or(false);
        let divergence_score =
            scoring::compute_divergence_score(alignment_map.and_then(|am| am.get(&m.gene_id)));
        let conserved_regions_present = align_entry
            .map(|a| {
                a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
            })
            .unwrap_or(false);
        let conserved_regions_score =
            compute_conserved_regions_score(alignment_map.and_then(|am| am.get(&m.gene_id)));
        let genomic_score = compute_genomic_score_with_cfg(
            genomic_map.and_then(|gm| gm.get(&m.gene_id)),
            genomic_cfg.and_then(|g| g.min_canonical),
            genomic_cfg.and_then(|g| g.max_noncanonical),
            genomic_cfg.and_then(|g| g.max_weird),
        );
        let genomic_present = genomic_map
            .map(|gm| gm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let h = adjust_homology_score(
            h_raw,
            divergence_score,
            divergence_present,
            l,
            length_present,
        );
        let presence = PillarPresence {
            homology: homology_present,
            intrinsic: intrinsic_present,
            taxonomy: taxonomy_present,
            domains: domains_present,
            domains_strength: domains_strength_present,
            length: length_present,
            orphan: orphan_present,
            subject_cov: subject_cov_present,
            termini: termini_present,
            divergence: divergence_present,
            conserved_regions: conserved_regions_present,
            genomic: genomic_present,
        };
        let missing = presence.missing_with_weights(&weights);
        let available_weight = weights.sum_available(&presence);
        let mut total_weight = available_weight;
        if !presence.homology {
            total_weight += weights.homology;
        }
        if !presence.domains {
            total_weight += weights.domains;
        }
        if hmmer_enabled && !presence.domains_strength {
            total_weight += weights.domains_strength;
        }
        let mut numerator = 0.0;
        if presence.homology {
            numerator += weights.homology * h;
        }
        if presence.intrinsic {
            numerator += weights.intrinsic * i;
        }
        if presence.taxonomy {
            numerator += weights.taxonomy * t;
        }
        if presence.domains {
            numerator += weights.domains * d;
        }
        if presence.domains_strength {
            numerator += weights.domains_strength * domains_strength_score;
        }
        if presence.length {
            numerator += weights.length * l;
        }
        if presence.orphan {
            numerator += weights.orphan * o;
        }
        if presence.subject_cov {
            numerator += weights.subject_cov * subject_cov_score;
        }
        if presence.termini {
            numerator += weights.termini * termini_score;
        }
        if presence.divergence {
            numerator += weights.divergence * divergence_score;
        }
        if presence.conserved_regions {
            numerator += weights.conserved_regions * conserved_regions_score;
        }
        if presence.genomic {
            numerator += weights.genomic * genomic_score;
        }
        let score = if total_weight > 1e-6 {
            (numerator / total_weight).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let sv = structvar_map.and_then(|sv| sv.get(&m.gene_id));
        let structvar_multiplier = if homology_present {
            compute_structvar_multiplier(sv)
        } else {
            1.0
        };
        let mut score = (score * structvar_multiplier).clamp(0.0, 1.0);
        // Hard caps for catastrophic structural-variation calls: ensure they can't be averaged away.
        if let (Some(cfg), Some(sv)) = (caps_cfg, sv) {
            let cap = match sv.classification.as_str() {
                "FusionPossible" => cfg.structvar_fusion_max,
                "SplitPossible" => cfg.structvar_split_max,
                "InternalDuplicationPossible" => cfg.structvar_dup_max,
                _ => None,
            };
            if let Some(max_score) = cap {
                score = score.min(max_score.clamp(0.0, 1.0));
            }
        }
        staging.push(ScoreTemp {
            gene_id: m.gene_id.clone(),
            raw_score: score,
            final_score: score,
            available_weight,
            missing,
        });
    }
    if !matches!(calibration.mode, CalibrationMode::Off) {
        if calibration_has_min_samples(&staging, calibration.min_samples, calibration.min_unique) {
            match calibration.mode {
                CalibrationMode::Percentile => apply_percentile_calibration(&mut staging),
                CalibrationMode::Isotonic => apply_isotonic_calibration(&mut staging),
                CalibrationMode::Off => {}
            }
        } else {
            let mut unique: std::collections::HashSet<u64> = std::collections::HashSet::new();
            for e in &staging {
                unique.insert(e.raw_score.to_bits());
            }
            log::warn!(
                "calibration skipped: samples={} unique={} (min_samples={}, min_unique={})",
                staging.len(),
                unique.len(),
                calibration.min_samples,
                calibration.min_unique
            );
        }
    }
    let mut out: HashMap<String, (f64, String)> = HashMap::new();
    let mut raw_map: HashMap<String, f64> = HashMap::new();
    for entry in staging {
        let base = if entry.available_weight < 1e-6 {
            "NoData"
        } else if entry.final_score >= th_high {
            "High"
        } else if entry.final_score >= th_med {
            "Medium"
        } else {
            "Low"
        };
        let classif = format_classification(base, &entry.missing, classify_no_data);
        raw_map.insert(entry.gene_id.clone(), entry.raw_score);
        out.insert(entry.gene_id, (entry.final_score, classif));
    }
    (out, raw_map)
}

#[allow(clippy::too_many_arguments)]
fn build_component_scores(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    taxonomy_enabled: bool,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
    alignment_map: Option<&HashMap<String, mafft::AlignmentMetrics>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    genomic_map: Option<&HashMap<String, genomic::GenomicMetrics>>,
    scoring: &Option<ScoringConfigOverride>,
) -> HashMap<String, ComponentScores> {
    let genomic_cfg = scoring.as_ref().and_then(|s| s.genomic.as_ref());
    let mut out: HashMap<String, ComponentScores> = HashMap::new();
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h_raw = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let taxonomy_component = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            if evidence
                .map(|ev| ev.considered > 0 || ev.top_hit.is_some())
                .unwrap_or(false)
            {
                Some(compute_taxonomy_score(evidence))
            } else {
                None
            }
        } else {
            None
        };
        let orphan_component = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        let subject_cov_score = compute_subject_cov_score(s);
        let termini_component = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            });
        let conserved_regions_component =
            alignment_map
                .and_then(|am| am.get(&m.gene_id))
                .and_then(|a| {
                    if a.mafft_enabled
                        && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
                    {
                        Some(compute_conserved_regions_score(Some(a)))
                    } else {
                        None
                    }
                });
        let divergence_component = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .map(|a| scoring::compute_divergence_score(Some(a)));
        let divergence_present = divergence_component.is_some();
        let length_score = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let length_present = len_map
            .map(|lm| lm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let h = adjust_homology_score(
            h_raw,
            divergence_component.unwrap_or(0.0),
            divergence_present,
            length_score,
            length_present,
        );
        let genomic_component = compute_genomic_score_with_cfg(
            genomic_map.and_then(|gm| gm.get(&m.gene_id)),
            genomic_cfg.and_then(|g| g.min_canonical),
            genomic_cfg.and_then(|g| g.max_noncanonical),
            genomic_cfg.and_then(|g| g.max_weird),
        );
        out.insert(
            m.gene_id.clone(),
            ComponentScores {
                homology: h,
                intrinsic: i,
                taxonomy: taxonomy_component,
                orphan: orphan_component,
                subject_cov: subject_cov_score,
                termini: termini_component,
                divergence: divergence_component,
                conserved_regions: conserved_regions_component,
                genomic: genomic_component,
            },
        );
    }
    out
}

fn export_high_sequences(
    out_dir: &str,
    override_path: Option<&str>,
    metrics: &[GeneMetrics],
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scores_map: &HashMap<String, (f64, String)>,
) -> Result<Option<(PathBuf, usize)>, Box<dyn std::error::Error>> {
    let target_path = override_path
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(out_dir).join("high_scoring.faa"));
    let mut writer: Option<BufWriter<File>> = None;
    let mut count = 0usize;
    for m in metrics {
        let Some((_, classif)) = scores_map.get(&m.gene_id) else {
            continue;
        };
        if classification_base(classif) != "High" {
            continue;
        }
        let Some((_im, seq)) = intrinsic.get(&m.gene_id) else {
            continue;
        };
        if writer.is_none() {
            if let Some(parent) = target_path.parent() {
                if !parent.as_os_str().is_empty() {
                    fs::create_dir_all(parent)?;
                }
            }
            writer = Some(BufWriter::new(File::create(&target_path)?));
        }
        if let Some(buf) = writer.as_mut() {
            write_fasta_record(buf, &m.gene_id, seq)?;
            count += 1;
        }
    }
    if let Some(mut buf) = writer {
        buf.flush()?;
        if count == 0 {
            drop(buf);
            fs::remove_file(&target_path).ok();
            Ok(None)
        } else {
            Ok(Some((target_path, count)))
        }
    } else {
        Ok(None)
    }
}

fn write_fasta_record<W: Write>(writer: &mut W, gene_id: &str, seq: &[u8]) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(gene_id.as_bytes())?;
    writer.write_all(b"\n")?;
    for chunk in seq.chunks(60) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
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
    let genomic_metrics = ctx.genomic_map.as_ref().and_then(|map| map.get(&m.gene_id));
    let genomic_score = comp_entry
        .map(|c| c.genomic)
        .unwrap_or_else(|| compute_genomic_score(genomic_metrics));
    let (base_final_score, classif) = ctx
        .scores_map
        .get(&m.gene_id)
        .cloned()
        .unwrap_or((0.0, "Low".to_string()));

    // Apply plugin penalty
    let final_score = (base_final_score - plugin_penalty).clamp(0.0, 1.0);
    // If penalty changed score significantly, we might want to downgrade classification,
    // but without thresholds we can't reliably. We assume user checks final_score.

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
            "structvar_multiplier": structvar_multiplier,
            "termini": termini_score,
            "genomic": genomic_score
        },
        "final_score_raw": raw_final_score,
        "final_score": final_score,
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
            "orf_start_score": intrinsic.orf_start_score,
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
    let divergence_score = comp_entry
        .and_then(|c| c.divergence)
        .unwrap_or_else(|| scoring::compute_divergence_score(aln));

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
            classif.clone(),
            format!("{:.4}", homology_score),
            format!("{:.4}", intrinsic_score),
            genomic_score_str.clone(),
            taxonomy_component_str,
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
            classif.clone(),
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
        classification: classif.clone(),
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

fn transcript_root(id: &str) -> Option<String> {
    if let Some((prefix, suffix)) = id.rsplit_once('.') {
        if suffix.len() >= 2
            && suffix.starts_with('t')
            && suffix[1..].chars().all(|c| c.is_ascii_digit())
        {
            return Some(prefix.to_string());
        }
    }
    None
}

fn propagate_transcript_taxonomy(
    map: &mut HashMap<String, Option<TaxonomyEvidence>>,
    metrics: &[GeneMetrics],
) {
    let mut groups: HashMap<String, Vec<String>> = HashMap::new();
    for m in metrics {
        if let Some(root) = transcript_root(&m.gene_id) {
            groups.entry(root).or_default().push(m.gene_id.clone());
        }
    }
    for (_, members) in groups {
        let mut best: Option<TaxonomyEvidence> = None;
        for gene in &members {
            if let Some(Some(ev)) = map.get(gene) {
                if matches!(
                    ev.detail,
                    TaxonomyDetail::Consensus | TaxonomyDetail::CoarseConsensus
                ) && best
                    .as_ref()
                    .is_none_or(|b| ev.congruence_score > b.congruence_score)
                {
                    best = Some(ev.clone());
                }
            }
        }
        let Some(best_ev) = best else { continue };
        for gene in &members {
            let entry = map.entry(gene.clone()).or_insert(None);
            let should_copy = match entry {
                Some(ev) => matches!(
                    ev.detail,
                    TaxonomyDetail::NoHits | TaxonomyDetail::InsufficientHits
                ),
                None => true,
            };
            if should_copy {
                let mut clone = best_ev.clone();
                clone.detail = TaxonomyDetail::Borrowed;
                clone.support = 0;
                clone.support_fraction = 0.0;
                clone.congruence_score = (clone.congruence_score * 0.9).clamp(0.0, 1.0);
                clone.contamination_score = (1.0 - clone.congruence_score).clamp(0.0, 1.0);
                *entry = Some(clone);
            }
        }
    }
}

// ... (rest of file unchanged for brevity in plan)
