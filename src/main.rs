use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use clap::{Args, Parser, Subcommand, ValueEnum};
use mimalloc::MiMalloc;
use serde::{Deserialize, Serialize};
use std::time::Instant;

mod checkpoint;
mod consensus;
mod diamond;
mod ecs;
mod hmmer;
mod length;
mod mafft;
mod metrics;
mod orf;
mod preflight;
mod provenance;
mod scoring;
mod structvar;
mod taxonomy;
use diamond::{
    blastp_once, cluster as diamond_cluster, linclust as diamond_linclust, parse_tsv_stats,
    DiamondConfig,
};
use ecs::{run_alignment_pipeline, run_hmmer_pipeline, run_scheduler, EcsConfig, GeneMetrics};
use hmmer::HmmscanSummary;
use mafft::{load_sequences_by_ids, AlignmentMetrics};
use metrics::compute_intrinsic;
use preflight::preflight;
use provenance::filehash_xx64;
use scoring::{compute_homology_score, compute_intrinsic_score, compute_taxonomy_score};
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
type LengthSummary = (f64, f64, f64, String);

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
    TaxonomyCache(TaxonomyCacheArgs),
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

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct AnalyzeArgs {
    #[arg(long)]
    fasta: Option<String>,
    #[arg(long)]
    db: Option<String>,
    #[arg(long, default_value_t = 16)]
    threads: usize,
    #[arg(long, default_value_t = 50)]
    approx_id: u32,
    #[arg(long, default_value_t = 80)]
    member_cover: u32,
    #[arg(long)]
    out: Option<String>,
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
    #[arg(long, default_value_t = true)]
    classify_no_data: bool,
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
    #[arg(long)]
    diamond_max_hsps: Option<usize>,
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
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct TaxonomyCacheArgs {
    #[arg(long)]
    input: String,
    #[arg(long)]
    output: String,
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
    structvar: Option<StructVarConfigOverride>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct TaxonomyConfigOverride {
    enabled: Option<bool>,
    min_support: Option<f64>,
    top_hits: Option<usize>,
    min_consensus: Option<usize>,
    coarse_rank_index: Option<usize>,
    coarse_min_support: Option<f64>,
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
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringConfigOverride {
    #[serde(default)]
    weights: HashMap<String, f64>,
    thresholds: Option<ScoringThresholds>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringThresholds {
    high: Option<f64>,
    medium: Option<f64>,
}

#[derive(Debug, Clone)]
struct EffectiveConfig {
    fasta: String,
    db: String,
    out: String,
    threads: usize,
    diamond_bin: String,
    reference_fasta: Option<String>,
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
            let file_cfg = if let Some(path) = cli.config.as_deref() {
                let text = fs::read_to_string(path)?;
                toml::from_str::<FileConfig>(&text)?
            } else {
                FileConfig::default()
            };
            let mut cfg = resolve_effective_config(&file_cfg, &args)?;

            fs::create_dir_all(&cfg.out)?;

            // Pre-flight: tool versions
            let tools = preflight(
                &cfg.diamond_bin,
                args.mafft_bin.as_deref(),
                args.hmmscan_bin.as_deref(),
            );

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

            // ECS: compute per-gene metrics from FASTA and DIAMOND
            let t_ecs = step_start("ecs", log_json);
            let metrics: Vec<GeneMetrics> = run_scheduler(EcsConfig {
                fasta_path: cfg.fasta.clone(),
                diamond_tsv: diamond_tsv.to_string_lossy().to_string(),
                threads: cfg.threads,
                log_json: matches!(args.log_format, LogFormat::Json),
            });
            let ecs_secs = step_finish("ecs", t_ecs, log_json);

            // Intrinsic metrics per gene
            let t_intrinsic = step_start("intrinsic", log_json);
            let intrinsic_map = compute_intrinsic_for_ids(&cfg.fasta, &metrics)?;
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
            }
            let mut panel_map: HashMap<String, Vec<String>> = HashMap::new();
            let mut panel_stats_rows: Vec<(String, consensus::PanelStats)> = Vec::new();
            let mut len_map: HashMap<String, LengthSummary> = HashMap::new();
            for m in &metrics {
                if let Some(hits) = grouped.get(&m.gene_id) {
                    let selection = consensus::select_panel(hits, &cons_cfg);
                    let panel_ids = selection.ids.clone();
                    panel_stats_rows.push((m.gene_id.clone(), selection.stats.clone()));
                    if !panel_ids.is_empty() {
                        panel_map.insert(m.gene_id.clone(), panel_ids.clone());
                        // Length consistency against available subject lengths from the selected panel only
                        let id_set: HashSet<String> = panel_ids.iter().cloned().collect();
                        let slens: Vec<usize> = hits
                            .iter()
                            .filter(|h| id_set.contains(&h.sseqid))
                            .map(|h| h.slen)
                            .filter(|&x| x > 0)
                            .collect();
                        if slens.len() >= cons_cfg.min_hits {
                            if let Some(lc) = length::compute_length_consistency(m.length, &slens) {
                                len_map.insert(
                                    m.gene_id.clone(),
                                    (lc.score, lc.z, lc.ratio, lc.class_),
                                );
                            }
                        }
                    }
                } else {
                    panel_stats_rows.push((
                        m.gene_id.clone(),
                        consensus::PanelStats {
                            total_hits: 0,
                            ..Default::default()
                        },
                    ));
                }
            }
            let consensus_secs = step_finish("consensus", t_consensus, log_json);
            if !panel_stats_rows.is_empty() {
                let panel_dbg = Path::new(&cfg.out).join("panel_debug.csv");
                let mut f = File::create(panel_dbg)?;
                writeln!(f, "gene_id,total_hits,phase,selected,filtered_hits,len_ratio_window_min,len_ratio_window_max,median_len_ratio,len_ratio_min,len_ratio_max,median_pident,high_identity_dropped")?;
                for (gid, st) in panel_stats_rows {
                    writeln!(
                        f,
                        "{},{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{}",
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
                    )?;
                }
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
            if let (Some(ref_fasta), Some(mafft_bin)) =
                (cfg.reference_fasta.as_ref(), args.mafft_bin.as_ref())
            {
                let t_mafft = step_start("mafft", log_json);
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
                    .or_else(|| file_cfg.mafft_threads_per_job)
                    .unwrap_or(default_per_job)
                    .clamp(1, total_threads);
                let default_workers = (total_threads / mafft_threads_per_job).max(1);
                let requested_workers = args
                    .mafft_max_jobs
                    .or_else(|| file_cfg.mafft_max_jobs)
                    .unwrap_or(default_workers)
                    .max(1);
                let mafft_workers = requested_workers.min(default_workers).max(1);
                std::env::set_var("MAFFT_THREADS", mafft_threads_per_job.to_string());

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

                log::info!("mafft jobs queued: {}", jobs.len());
                if !jobs.is_empty() {
                    log::info!("mafft jobs queued: {}", jobs.len());
                    alignment_map = run_alignment_pipeline(
                        mafft_bin,
                        jobs,
                        query_seq_map,
                        ref_seqs,
                        mafft_workers,
                    );
                }
                let _ = step_finish("mafft", t_mafft, log_json);
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

            // Optional HMMER/Pfam domain summary (JSONL)
            let mut hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut _ref_hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut domains_arch_map: HashMap<String, f64> = HashMap::new();
            let mut domains_arch_dbg: Vec<DomainsArchDebugRow> = Vec::new();
            let mut orphan_map: HashMap<String, hmmer::OrphanAnalysis> = HashMap::new();
            let mut orphan_analysis_enabled = false;
            if let (Some(hmm), Some(pfam_db)) = (
                args.hmmscan_bin.as_ref(),
                args.pfam_db
                    .as_ref()
                    .or(file_cfg.pfam_db.as_ref())
                    .or(file_cfg.pfam_db.as_ref()),
            ) {
                let t_hmmer = step_start("hmmer", log_json);
                let items: Vec<(String, Vec<u8>)> = intrinsic_map
                    .iter()
                    .map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone()))
                    .collect();
                let top_n = args
                    .hmmer_top_n
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.top_n))
                    .unwrap_or(5);
                let hmmer_threads = args
                    .hmmer_threads
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.threads))
                    .unwrap_or(cfg.threads);
                let ievalue = args
                    .hmmer_ievalue
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.ievalue));
                if !items.is_empty() {
                    log::info!("hmmscan jobs queued: {}", items.len());
                    hmmsum_map =
                        run_hmmer_pipeline(hmm, pfam_db, items, hmmer_threads, top_n, ievalue);
                }
                let orphan_cfg = file_cfg
                    .hmmer
                    .as_ref()
                    .and_then(|h| h.orphan_analysis)
                    .unwrap_or(true);
                orphan_analysis_enabled = orphan_cfg && !args.disable_orphan_analysis;
                if orphan_analysis_enabled {
                    orphan_map = hmmsum_map
                        .iter()
                        .filter(|(_, summary)| !summary.hits.is_empty())
                        .map(|(gid, summary)| (gid.clone(), hmmer::analyze_orphan_domains(summary)))
                        .collect();
                } else {
                    orphan_map.clear();
                }
                if let Some(ref_fasta) = cfg.reference_fasta.as_ref() {
                    let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
                    let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                    if !ref_seqs.is_empty() {
                        let ref_items: Vec<(String, Vec<u8>)> = ref_seqs.into_iter().collect();
                        let ref_ievalue = args
                            .hmmer_ref_ievalue
                            .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.ref_ievalue));
                        if !ref_items.is_empty() {
                            log::info!("hmmscan reference jobs queued: {}", ref_items.len());
                            _ref_hmmsum_map = run_hmmer_pipeline(
                                hmm,
                                pfam_db,
                                ref_items,
                                hmmer_threads,
                                top_n,
                                ref_ievalue,
                            );
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
                let _ = step_finish("hmmer", t_hmmer, log_json);
            }

            // Emit outputs
            let t_emit = step_start("emit_outputs", log_json);
            let qlen_map: HashMap<String, usize> = metrics
                .iter()
                .map(|m| (m.gene_id.clone(), m.length))
                .collect();
            let stats = parse_tsv_stats(&diamond_tsv, Some(&qlen_map)).unwrap_or_default();
            let checksums = collect_checksums(&cfg)?;
            let snapshot = build_config_snapshot(&cfg, &args, &file_cfg);
            write_run_manifest(&cfg, &tools, &checksums, &snapshot)?;
            // Optional taxonomy resolution of top hits → taxid/name/lineage per gene
            let mut taxsum_map: HashMap<String, Option<TaxonomyEvidence>> = HashMap::new();
            if args.enable_taxonomy {
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
                if let Some(resolver) = resolver {
                    for m in &metrics {
                        let hit_ids: Vec<String> = grouped
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
                            if let Some(acc) =
                                stats.get(&m.gene_id).and_then(|s| s.top_sseqid.clone())
                            {
                                evidence.top_hit = resolver.lookup(&acc);
                            }
                        }
                        taxsum_map.insert(m.gene_id.clone(), Some(evidence));
                    }
                    propagate_transcript_taxonomy(&mut taxsum_map, &metrics);
                } else {
                    for m in &metrics {
                        taxsum_map.insert(m.gene_id.clone(), None);
                    }
                }
            }

            let t_scoring = step_start("scoring", log_json);
            let scores_map = build_scores_map(
                &metrics,
                &stats,
                &intrinsic_map,
                &file_cfg.scoring,
                args.enable_taxonomy,
                Some(&domains_arch_map),
                Some(&len_map),
                if orphan_analysis_enabled {
                    Some(&orphan_map)
                } else {
                    None
                },
                if args.enable_taxonomy {
                    Some(&taxsum_map)
                } else {
                    None
                },
            );
            let comp_map = build_component_scores(
                &metrics,
                &stats,
                &intrinsic_map,
                args.enable_taxonomy,
                if orphan_analysis_enabled {
                    Some(&orphan_map)
                } else {
                    None
                },
                if args.enable_taxonomy {
                    Some(&taxsum_map)
                } else {
                    None
                },
            );
            let _ = step_finish("scoring", t_scoring, log_json);
            let mafft_missing_exon_threshold = args
                .alignment_missing_exon
                .or(file_cfg.alignment_missing_exon)
                .unwrap_or(30);
            let mafft_retained_intron_threshold = args
                .alignment_retained_intron
                .or(file_cfg.alignment_retained_intron)
                .unwrap_or(30);
            write_jsonl_metrics(
                &cfg.out,
                &metrics,
                &stats,
                &intrinsic_map,
                &alignment_map,
                &hmmsum_map,
                args.coverage_delta_threshold,
                &file_cfg.scoring,
                args.enable_taxonomy,
                &taxsum_map,
                Some(&domains_arch_map),
                Some(&len_map),
                if orphan_analysis_enabled {
                    Some(&orphan_map)
                } else {
                    None
                },
                Some(&structvar_map),
                mafft_missing_exon_threshold,
                mafft_retained_intron_threshold,
            )?;
            write_csv_metrics(
                &cfg.out,
                &metrics,
                &stats,
                args.coverage_delta_threshold,
                &scores_map,
                &alignment_map,
                args.enable_taxonomy,
                &taxsum_map,
                &hmmsum_map,
                args.csv_verbose,
                Some(&comp_map),
                args.classify_no_data,
                Some(&domains_arch_map),
                Some(&len_map),
                if orphan_analysis_enabled {
                    Some(&orphan_map)
                } else {
                    None
                },
                Some(&structvar_map),
                mafft_missing_exon_threshold,
                mafft_retained_intron_threshold,
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
            let runm = RunMetrics {
                schema_version: "1.0",
                steps,
                totals,
            };
            write_run_metrics(&cfg.out, &runm)?;
            Ok(())
        }
        Commands::TaxonomyCache(t) => {
            let count = taxonomy::write_cache_from_fasta(&t.input, &t.output)
                .map_err(|e| format!("taxonomy cache failed: {}", e))?;
            eprintln!("wrote {} taxonomy entries to {}", count, t.output);
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

    Ok(EffectiveConfig {
        fasta,
        db,
        out,
        threads,
        diamond_bin,
        reference_fasta,
    })
}

#[derive(Serialize)]
struct Checksums {
    fasta_xx64: Option<String>,
    db_xx64: Option<String>,
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

#[derive(Serialize)]
struct ConfigSnapshot<'a> {
    fasta: &'a str,
    db: &'a str,
    out: &'a str,
    threads: usize,
    diamond_bin: &'a str,
    reference_fasta: Option<&'a str>,
    alignment_top_hits: usize,
    alignment_strategy: String,
    coverage_delta_threshold: f64,
    log_format: String,
    scoring_weights: std::collections::HashMap<String, f64>,
}

fn write_run_manifest(
    cfg: &EffectiveConfig,
    tools: &preflight::ToolVersions,
    sums: &Checksums,
    snapshot: &ConfigSnapshot,
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
    }
    let manifest = Manifest {
        schema_version: "1.0",
        tool: "AnnoQC",
        diamond_version: tools.diamond.as_deref().unwrap_or_default(),
        mafft_version: tools.mafft.as_deref(),
        hmmscan_version: tools.hmmscan.as_deref(),
        fasta_xx64: sums.fasta_xx64.as_deref(),
        db_xx64: sums.db_xx64.as_deref(),
        config: snapshot,
    };
    let path = Path::new(&cfg.out).join("run.json");
    let text = serde_json::to_string_pretty(&manifest)?;
    fs::write(path, text)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_jsonl_metrics(
    out_dir: &str,
    metrics: &[ecs::GeneMetrics],
    stats: &std::collections::HashMap<String, diamond::DiamondHitStats>,
    intrinsic_map: &std::collections::HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    alignment_map: &std::collections::HashMap<String, mafft::AlignmentMetrics>,
    hmmsum_map: &std::collections::HashMap<String, hmmer::HmmscanSummary>,
    cov_delta_thresh: f64,
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
    taxsum_map: &HashMap<String, Option<TaxonomyEvidence>>,
    arch_map: Option<&HashMap<String, f64>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    sv_map: Option<&HashMap<String, structvar::StructVar>>,
    mafft_missing_exon_thresh: usize,
    mafft_retained_intron_thresh: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("qc_report.jsonl");
    let mut f = File::create(path)?;
    for m in metrics {
        let warnings: Vec<String> = if m.hits == 0 {
            vec!["No DIAMOND hits".to_string()]
        } else {
            Vec::new()
        };
        let summary = stats.get(&m.gene_id);
        let intrinsic = intrinsic_map
            .get(&m.gene_id)
            .map(|t| t.0.clone())
            .unwrap_or_default();
        let aln = alignment_map.get(&m.gene_id);
        let homology_score = compute_homology_score(summary);
        let intrinsic_score = compute_intrinsic_score(&intrinsic);
        let taxonomy_evidence = taxsum_map.get(&m.gene_id).and_then(|x| x.as_ref());
        let taxonomy_score = if taxonomy_enabled {
            Some(compute_taxonomy_score(taxonomy_evidence))
        } else {
            None
        };
        let weights = scoring_weights(scoring);
        let domains_arch_score = arch_map
            .and_then(|am| am.get(&m.gene_id).cloned())
            .unwrap_or(0.0);
        let length_score = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let orphan_score = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|o| o.score)
            .unwrap_or(1.0);
        let denom = weights.sum().max(1e-6);
        let final_score = (weights.homology * homology_score
            + weights.intrinsic * intrinsic_score
            + weights.taxonomy * taxonomy_score.unwrap_or(0.0)
            + weights.domains * domains_arch_score
            + weights.length * length_score
            + weights.orphan * orphan_score)
            / denom;
        let fusion_split_flag = summary
            .map(|s| s.coverage_delta > cov_delta_thresh)
            .unwrap_or(false);
        let taxonomy = if taxonomy_enabled {
            if let Some(Some(ev)) = taxsum_map.get(&m.gene_id) {
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
        let cur_gene = &m.gene_id;
        // Build alignment/struct-var warnings
        let mut warn_extra: Vec<String> = Vec::new();
        if let Some(a) = aln {
            if a.missing_exon_run >= mafft_missing_exon_thresh {
                warn_extra.push("MissingExonPossible".into());
            }
            if a.retained_intron_run >= mafft_retained_intron_thresh {
                warn_extra.push("RetainedIntronPossible".into());
            }
        }
        let sv_obj = sv_map.and_then(|mm| mm.get(cur_gene));
        if let Some(sv) = sv_obj {
            match sv.classification.as_str() {
                "FusionPossible" => warn_extra.push("FusionPossible".into()),
                "SplitPossible" => warn_extra.push("SplitPossible".into()),
                "InternalDuplicationPossible" => {
                    warn_extra.push("InternalDuplicationPossible".into())
                }
                _ => {}
            }
            for w in &sv.warnings {
                warn_extra.push(w.clone());
            }
        }
        let domains = hmmsum_map.get(cur_gene).map(|d| {
            let score = if let Some(ev) = d.top_evalue {
                let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
                (le / 20.0).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let arch_score = arch_map
                .and_then(|am| am.get(cur_gene))
                .cloned()
                .unwrap_or(0.0);
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
        let length_block = len_map
            .and_then(|lm| lm.get(&m.gene_id))
            .map(|(s, z, r, c)| {
                serde_json::json!({
                    "length_score": s,
                    "length_z": z,
                    "length_ratio": r,
                    "length_class": c,
                })
            });
        let orphan_entry = orphan_map.and_then(|om| om.get(cur_gene));
        let record = serde_json::json!({
                "gene_id": m.gene_id,
                "taxonomy": taxonomy,
                "score_components": {"taxonomy": taxonomy_score, "homology": homology_score, "intrinsic": intrinsic_score, "domains": domains_arch_score, "length": length_score, "orphan": orphan_score},
                "final_score": final_score,
                "homology": {
                "hits_count": m.hits,
                "top_hit": summary.and_then(|s| s.top_sseqid.clone()),
                "top_bitscore": summary.map(|s| s.top_bitscore),
                "top_evalue": summary.as_ref().map(|s| s.top_evalue.clone()),
                "top_qcov": summary.map(|s| s.top_qcov),
                "top_scov": summary.map(|s| s.top_scov),
                "bitscore_density": summary.map(|s| if s.top_len>0 { s.top_bitscore / s.top_len as f64 } else { 0.0 }),
                "coverage_delta": summary.map(|s| s.coverage_delta),
                "coverage_ratio": summary.map(|s| s.coverage_ratio),
                "fusion_split_flag": fusion_split_flag,
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
            "warnings": if warn_extra.is_empty() { warnings } else { warn_extra },
        });
        writeln!(f, "{}", serde_json::to_string(&record)?)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
fn write_csv_metrics(
    out_dir: &str,
    metrics: &[ecs::GeneMetrics],
    stats: &std::collections::HashMap<String, diamond::DiamondHitStats>,
    cov_delta_thresh: f64,
    scores: &HashMap<String, (f64, String)>,
    alignment_map: &std::collections::HashMap<String, mafft::AlignmentMetrics>,
    taxonomy_enabled: bool,
    taxsum_map: &HashMap<String, Option<TaxonomyEvidence>>,
    hmmsum_map: &std::collections::HashMap<String, hmmer::HmmscanSummary>,
    csv_verbose: bool,
    comp_map: Option<&HashMap<String, ComponentScores>>,
    classify_no_data: bool,
    arch_map: Option<&HashMap<String, f64>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    sv_map: Option<&HashMap<String, structvar::StructVar>>,
    mafft_missing_exon_thresh: usize,
    mafft_retained_intron_thresh: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("qc_summary.csv");
    let mut f = File::create(path)?;
    if csv_verbose {
        writeln!(f, "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,homology_score,intrinsic_score,taxonomy_score,domains_score,domains_arch_score,orphan_domain_score,length_score,length_ratio,length_class,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,missing_exon_run,retained_intron_run,start_concordance,start_class,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_contamination,taxonomy_support,taxonomy_considered,consensus_taxon,taxonomy_status,warnings")?;
    } else {
        writeln!(f, "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,domains_score,domains_arch_score,orphan_domain_score,structvar_class,structvar_gap,structvar_left_len,structvar_right_len,structvar_cov_left,structvar_cov_right,orphan_status,taxonomy_score,taxonomy_contamination,taxonomy_support,taxonomy_considered,consensus_taxon,taxonomy_status,warnings")?;
    }
    for m in metrics {
        let mut warning_msgs: Vec<String> = Vec::new();
        if m.hits == 0 {
            warning_msgs.push("No DIAMOND hits".into());
        }
        let s = stats.get(&m.gene_id);
        let (top_hit, top_bitscore, top_evalue, top_qcov, top_scov, bsd) = if let Some(s) = s {
            let bsd = if s.top_len > 0 {
                s.top_bitscore / s.top_len as f64
            } else {
                0.0
            };
            (
                s.top_sseqid.clone().unwrap_or_default(),
                s.top_bitscore,
                s.top_evalue.clone(),
                s.top_qcov,
                s.top_scov,
                bsd,
            )
        } else {
            (String::new(), 0.0, String::new(), 0.0, 0.0, 0.0)
        };
        let cov_delta = s.map(|x| x.coverage_delta).unwrap_or(0.0);
        let cov_ratio = s.map(|x| x.coverage_ratio).unwrap_or(0.0);
        let fusion = if cov_delta > cov_delta_thresh { 1 } else { 0 };
        let (final_score, mut classif) = scores
            .get(&m.gene_id)
            .cloned()
            .unwrap_or((0.0, String::new()));
        let aln = alignment_map.get(&m.gene_id);
        let (
            mafft_enabled,
            conserved,
            pid,
            seqs_aln,
            qgap,
            gap_runs,
            max_gap,
            start_conc,
            start_class,
            missing_run,
            intron_run,
        ) = if let Some(a) = aln {
            (
                a.mafft_enabled as i32,
                a.conserved_fraction,
                a.pairwise_identity,
                a.sequences_aligned,
                a.query_gap_fraction,
                a.gap_run_count,
                a.max_gap_run,
                a.start_concordance,
                a.start_class.clone(),
                a.missing_exon_run,
                a.retained_intron_run,
            )
        } else {
            (0, 0.0, 0.0, 0, 0.0, 0, 0, 0.0, String::new(), 0, 0)
        };
        if mafft_enabled == 1 {
            if missing_run >= mafft_missing_exon_thresh {
                warning_msgs.push("MissingExonPossible".into());
            }
            if intron_run >= mafft_retained_intron_thresh {
                warning_msgs.push("RetainedIntronPossible".into());
            }
        }
        let (
            structvar_class,
            structvar_gap,
            structvar_left_len,
            structvar_right_len,
            structvar_cov_left,
            structvar_cov_right,
        ) = if let Some(sv) = sv_map.and_then(|sm| sm.get(&m.gene_id)) {
            match sv.classification.as_str() {
                "FusionPossible" => warning_msgs.push("FusionPossible".into()),
                "SplitPossible" => warning_msgs.push("SplitPossible".into()),
                "InternalDuplicationPossible" => {
                    warning_msgs.push("InternalDuplicationPossible".into())
                }
                _ => {}
            }
            for w in &sv.warnings {
                warning_msgs.push(w.clone());
            }
            (
                sv.classification.clone(),
                sv.fusion_gap.map(|v| v.to_string()).unwrap_or_default(),
                sv.fusion_left_len
                    .map(|v| v.to_string())
                    .unwrap_or_default(),
                sv.fusion_right_len
                    .map(|v| v.to_string())
                    .unwrap_or_default(),
                sv.fusion_cover_fracs
                    .map(|t| format!("{:.3}", t.0))
                    .unwrap_or_default(),
                sv.fusion_cover_fracs
                    .map(|t| format!("{:.3}", t.1))
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
        let (
            taxonomy_score_value,
            taxonomy_contamination_field,
            taxonomy_support_field,
            taxonomy_considered_field,
            taxonomy_status_str,
            consensus_label,
        ) = if taxonomy_enabled {
            if let Some(Some(ev)) = taxsum_map.get(&m.gene_id) {
                let consensus_label = ev
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
                    compute_taxonomy_score(Some(ev)),
                    format!("{:.4}", ev.contamination_score),
                    ev.support.to_string(),
                    ev.considered.to_string(),
                    ev.detail.to_string(),
                    consensus_label,
                )
            } else {
                (
                    0.0,
                    String::new(),
                    String::new(),
                    String::new(),
                    "NoResolver".to_string(),
                    String::new(),
                )
            }
        } else {
            (
                0.0,
                String::new(),
                String::new(),
                String::new(),
                "disabled".to_string(),
                String::new(),
            )
        };
        let domains_score_field = if let Some(d) = hmmsum_map.get(&m.gene_id) {
            let s = if let Some(ev) = d.top_evalue {
                let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
                (le / 20.0).clamp(0.0, 1.0)
            } else {
                0.0
            };
            format!("{:.4}", s)
        } else {
            String::new()
        };
        let domains_arch_field = arch_map
            .and_then(|am| am.get(&m.gene_id))
            .map(|v| format!("{:.4}", v))
            .unwrap_or_default();
        let (orphan_status_str, orphan_score_field) = if let Some(map) = orphan_map {
            if let Some(oa) = map.get(&m.gene_id) {
                (oa.status.as_str().to_string(), format!("{:.4}", oa.score))
            } else {
                (String::new(), String::new())
            }
        } else {
            (String::new(), String::new())
        };
        // Optional NoData classification override
        if classify_no_data && m.hits == 0 {
            let nonzero_domains = hmmsum_map
                .get(&m.gene_id)
                .map(|d| d.hits_count > 0)
                .unwrap_or(false);
            if !nonzero_domains {
                classif = "NoData".to_string();
            }
        }
        let warnings_field = if warning_msgs.is_empty() {
            String::new()
        } else {
            warning_msgs.join(";")
        };
        if csv_verbose {
            let (hs, is, ts_opt, orphan_component) = comp_map
                .and_then(|cm| cm.get(&m.gene_id))
                .map(|c| (c.homology, c.intrinsic, c.taxonomy, c.orphan))
                .unwrap_or((0.0, 0.0, None, 1.0));
            let ts_str = if taxonomy_enabled {
                format!("{:.4}", ts_opt.unwrap_or(0.0))
            } else {
                String::new()
            };
            let (len_s, _len_z, len_r, len_c) = len_map
                .and_then(|lm| lm.get(&m.gene_id))
                .cloned()
                .unwrap_or((0.0, 0.0, 0.0, String::new()));
            let len_ratio_str = if len_r > 0.0 {
                format!("{:.3}", len_r)
            } else {
                String::new()
            };
            let orphan_component_str = if orphan_map.is_some() {
                if orphan_score_field.is_empty() {
                    format!("{:.4}", orphan_component)
                } else {
                    orphan_score_field.clone()
                }
            } else {
                String::new()
            };
            let row = vec![
                m.gene_id.clone(),
                m.hits.to_string(),
                top_hit.clone(),
                format!("{:.3}", top_bitscore),
                top_evalue.clone(),
                format!("{:.3}", top_qcov),
                format!("{:.3}", top_scov),
                format!("{:.3}", bsd),
                format!("{:.3}", cov_delta),
                format!("{:.3}", cov_ratio),
                fusion.to_string(),
                format!("{:.3}", final_score),
                classif.clone(),
                format!("{:.4}", hs),
                format!("{:.4}", is),
                ts_str,
                domains_score_field.clone(),
                domains_arch_field.clone(),
                orphan_component_str,
                format!("{:.4}", len_s),
                len_ratio_str,
                len_c,
                mafft_enabled.to_string(),
                format!("{:.3}", conserved),
                format!("{:.3}", pid),
                seqs_aln.to_string(),
                format!("{:.3}", qgap),
                gap_runs.to_string(),
                max_gap.to_string(),
                missing_run.to_string(),
                intron_run.to_string(),
                format!("{:.3}", start_conc),
                start_class,
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
                consensus_label.clone(),
                taxonomy_status_str.clone(),
                warnings_field.clone(),
            ]
            .join(",");
            writeln!(f, "{}", row)?;
            continue;
        }
        let row = vec![
            m.gene_id.clone(),
            m.hits.to_string(),
            top_hit,
            format!("{:.3}", top_bitscore),
            top_evalue,
            format!("{:.3}", top_qcov),
            format!("{:.3}", top_scov),
            format!("{:.3}", bsd),
            format!("{:.3}", cov_delta),
            format!("{:.3}", cov_ratio),
            fusion.to_string(),
            format!("{:.3}", final_score),
            classif,
            mafft_enabled.to_string(),
            format!("{:.3}", conserved),
            format!("{:.3}", pid),
            seqs_aln.to_string(),
            format!("{:.3}", qgap),
            gap_runs.to_string(),
            max_gap.to_string(),
            domains_score_field,
            domains_arch_field,
            orphan_score_field,
            structvar_class,
            structvar_gap,
            structvar_left_len,
            structvar_right_len,
            structvar_cov_left,
            structvar_cov_right,
            orphan_status_str,
            format!("{:.4}", taxonomy_score_value),
            taxonomy_contamination_field,
            taxonomy_support_field,
            taxonomy_considered_field,
            consensus_label,
            taxonomy_status_str,
            warnings_field,
        ]
        .join(",");
        writeln!(f, "{}", row)?;
    }
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
    checkpoint::run_step(db_done, p.resume, "diamond makedb", prep_json, || {
        if db_out.exists() {
            log::info!(
                "makedb: {} exists; rebuilding due to resume=false",
                db_out.display()
            );
        }
        let status = std::process::Command::new(&diamond_bin)
            .arg("makedb")
            .arg("--in")
            .arg(&p.fasta)
            .arg("--db")
            .arg(&p.db_out)
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("diamond makedb failed with status {}", status));
        }
        Ok(())
    })?;

    // Step 3: linclust -> clusters
    let clusters = Path::new("clusters");
    let clusters_done = Path::new("clusters.done");
    checkpoint::run_step(
        clusters_done,
        p.resume,
        "diamond linclust",
        prep_json,
        || {
            diamond_linclust(&diamond_bin, &p.fasta, clusters, p.approx_id, p.threads)
                .map_err(|e| format!("linclust failed: {}", e))
        },
    )?;

    // Step 4: cluster (sensitive) → clusters.realign
    let realign = Path::new("clusters.realign");
    let realign_done = Path::new("clusters.realign.done");
    checkpoint::run_step(realign_done, p.resume, "diamond cluster", prep_json, || {
        diamond_cluster(&diamond_bin, &p.fasta, realign, p.approx_id, p.threads)
            .map_err(|e| format!("cluster failed: {}", e))
    })?;

    // Step 5: recluster (placeholder) → clusters.recluster
    let recluster = Path::new("clusters.recluster");
    let recluster_done = Path::new("clusters.recluster.done");
    checkpoint::run_step(
        recluster_done,
        p.resume,
        "diamond recluster (placeholder)",
        prep_json,
        || {
            if recluster.exists() {
                std::fs::remove_file(recluster).ok();
            }
            std::fs::copy(realign, recluster)
                .map(|_| ())
                .map_err(|e| e.to_string())
        },
    )?;
    Ok(())
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

#[derive(Clone, Copy, Debug)]
struct WeightSet {
    homology: f64,
    intrinsic: f64,
    taxonomy: f64,
    domains: f64,
    length: f64,
    orphan: f64,
}

impl WeightSet {
    fn sum(&self) -> f64 {
        self.homology + self.intrinsic + self.taxonomy + self.domains + self.length + self.orphan
    }
}

fn scoring_weights(scoring: &Option<ScoringConfigOverride>) -> WeightSet {
    let mut ws = WeightSet {
        homology: 0.6,
        intrinsic: 0.4,
        taxonomy: 0.0,
        domains: 0.0,
        length: 0.0,
        orphan: 0.0,
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
        if let Some(v) = cfg.weights.get("length") {
            ws.length = *v;
        }
        if let Some(v) = cfg.weights.get("orphan") {
            ws.orphan = *v;
        }
    }
    ws
}

#[derive(Clone, Debug, Default)]
struct ComponentScores {
    homology: f64,
    intrinsic: f64,
    taxonomy: Option<f64>,
    orphan: f64,
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
) -> ConfigSnapshot<'a> {
    let mut weights = HashMap::new();
    if let Some(sc) = &file_cfg.scoring {
        weights = sc.weights.clone();
    } else {
        weights.insert("homology".to_string(), 0.6);
        weights.insert("intrinsic".to_string(), 0.4);
    }
    ConfigSnapshot {
        fasta: &cfg.fasta,
        db: &cfg.db,
        out: &cfg.out,
        threads: cfg.threads,
        diamond_bin: &cfg.diamond_bin,
        reference_fasta: cfg.reference_fasta.as_deref(),
        alignment_top_hits: args.alignment_top_hits,
        alignment_strategy: format!("{:?}", args.alignment_strategy),
        coverage_delta_threshold: args.coverage_delta_threshold,
        log_format: format!("{:?}", args.log_format),
        scoring_weights: weights,
    }
}

#[allow(clippy::too_many_arguments)]
fn build_scores_map(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
    arch_map: Option<&HashMap<String, f64>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
) -> HashMap<String, (f64, String)> {
    let weights = scoring_weights(scoring);
    let (th_high, th_med) = scoring_thresholds(scoring);
    let mut out: HashMap<String, (f64, String)> = HashMap::new();
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let t = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            compute_taxonomy_score(evidence)
        } else {
            0.0
        };
        let d = arch_map
            .and_then(|am| am.get(&m.gene_id))
            .cloned()
            .unwrap_or(0.0);
        let l = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let o = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        let denom = weights.sum().max(1e-6);
        let score = (weights.homology * h
            + weights.intrinsic * i
            + weights.taxonomy * t
            + weights.domains * d
            + weights.length * l
            + weights.orphan * o)
            / denom;
        let classif = if score >= th_high {
            "High"
        } else if score >= th_med {
            "Medium"
        } else {
            "Low"
        };
        out.insert(m.gene_id.clone(), (score, classif.to_string()));
    }
    out
}

fn build_component_scores(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    taxonomy_enabled: bool,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
) -> HashMap<String, ComponentScores> {
    let mut out: HashMap<String, ComponentScores> = HashMap::new();
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let taxonomy_component = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            Some(compute_taxonomy_score(evidence))
        } else {
            None
        };
        let orphan_component = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        out.insert(
            m.gene_id.clone(),
            ComponentScores {
                homology: h,
                intrinsic: i,
                taxonomy: taxonomy_component,
                orphan: orphan_component,
            },
        );
    }
    out
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
