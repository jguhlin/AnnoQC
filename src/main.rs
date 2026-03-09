use std::io::Write;

use clap::{Args, Parser, Subcommand};
use mimalloc::MiMalloc;
use serde::{Deserialize, Serialize};
mod alignment_support;
mod analysis_data_support;
mod analysis_extension_support;
mod analyze;
mod app_config_support;
mod app_enums;
mod app_intrinsic_support;
mod app_support;
mod app_taxonomy_count_support;
mod debug_outputs;
use alignment_support::{build_alignment_setup, AlignmentOverrideInputs};
use analysis_support::{
    build_features_string, load_extensions, load_genomic_context, load_rnaseq_data, FeatureSummary,
};
pub(crate) use app_enums::{
    AlignerCliBackend, AlignmentStrategy, CalibrationMode, DiamondMode, LogFormat, Mode,
    ReportFormat,
};
pub(crate) use app_support::{
    collect_checksums, compute_intrinsic_for_ids, resolve_calibration_settings,
    resolve_effective_config, CalibrationSettings, Checksums, EffectiveConfig,
};
mod analysis_support;
mod checkpoint;
mod clusters;
mod config_types;
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
mod render_context_support;
mod render_job_support;
mod render_record;
mod render_support;
mod render_types;
mod reporting_metrics;
mod reporting_sidecars;
mod reporting_support;
mod rhai_rules;
mod rnaseq;
mod scoring;
mod scoring_builder;
mod scoring_calibration;
mod scoring_export;
mod scoring_support;
mod scoring_thresholds;
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
    HeavyPipelineConfig,
};
use hmmer::HmmscanSummary;
use hmmer_support::{
    build_hmmer_setup, postprocess_hmmer, DomainsArchDebugRow, HmmerOverrideInputs,
    HmmerPostProcessInputs,
};
use mafft::{AlignerBackend, AlignmentMetrics};
use orchestration::{
    step_finish, step_start, write_run_metrics, write_run_summary, write_slowest_genes,
};
use preflight::preflight;
use profiles::Profile;
use provenance::filehash_xx64;
use render_support::{build_render_context, resolve_render_max_jobs, RenderContextInputs};
pub(crate) use render_types::{RenderContext, RenderedRecord, ScoreCard};
use reporting_support::{
    augment_totals_from_taxonomy, build_base_totals, build_run_metrics as build_run_metrics_obj,
    build_run_summary as build_run_summary_obj, build_steps, write_domains_arch_debug,
    write_structvar_summary, RefprotSelectionConfig, RunSummaryInputs,
};
use scoring::{
    compute_conserved_regions_score, compute_domains_strength_score, compute_genomic_score,
    compute_genomic_score_with_cfg, compute_homology_score, compute_intrinsic_score,
    compute_rnaseq_score, compute_structvar_multiplier, compute_subject_cov_penalty,
    compute_subject_cov_score, compute_taxonomy_score,
};
use stage_config_support::{
    resolve_consensus_and_refprot_config, ConsensusOverrideInputs, RefProtOverrideInputs,
};
use structvar_support::{analyze_structvar, resolve_thresholds, StructVarOverrideInputs};
use taxonomy::TaxonomyConsensusConfig;
use taxonomy_support::{build_taxonomy_setup, TaxonomyWarnings};

pub(crate) const OUTPUT_SCHEMA_VERSION: &str = "1.1";
const PLUGIN_SCHEMA_VERSION: &str = "v1";
pub(crate) const TOOL_NAME: &str = env!("CARGO_PKG_NAME");
pub(crate) const TOOL_VERSION: &str = env!("CARGO_PKG_VERSION");

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

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
            app_support::run_taxonomy_count(args)?;
            Ok(())
        }
    }
}

#[cfg(test)]
mod main_tests;
