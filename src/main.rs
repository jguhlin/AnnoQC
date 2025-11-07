use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use clap::{Args, Parser, Subcommand, ValueEnum};
use mimalloc::MiMalloc;
use serde::{Deserialize, Serialize};
use std::time::Instant;

mod checkpoint;
mod diamond;
mod ecs;
mod consensus;
mod hmmer;
mod mafft;
mod metrics;
mod preflight;
mod provenance;
mod scoring;
mod taxonomy;
use diamond::{
    blastp_once, cluster as diamond_cluster, linclust as diamond_linclust, parse_tsv_stats,
    DiamondConfig,
};
use ecs::{run_scheduler, EcsConfig, GeneMetrics};
use hmmer::HmmscanSummary;
use mafft::{load_sequences_by_ids, run_mafft, AlignmentMetrics};
use metrics::compute_intrinsic;
use preflight::preflight;
use provenance::filehash_xx64;
use scoring::{
    combine_scores3, compute_homology_score, compute_intrinsic_score, compute_taxonomy_score,
};
use taxonomy::{TaxonSummary, TaxonomyResolver};

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
    reference_fasta: Option<String>,
    #[arg(long, default_value_t = 5)]
    alignment_top_hits: usize,
    #[arg(long, value_enum, default_value_t = AlignmentStrategy::Auto)]
    alignment_strategy: AlignmentStrategy,
    #[arg(long, default_value_t = 0.8)]
    conserved_identity_min: f64,
    #[arg(long, default_value_t = 8)]
    batch_size: usize,
    #[arg(long)]
    enable_taxonomy: bool,
    #[arg(long)]
    taxonomy_min_support: Option<f64>,
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
    #[arg(long, value_enum, default_value_t = DiamondMode::Auto)]
    diamond_mode: DiamondMode,
    #[arg(long, value_enum, default_value_t = LogFormat::Text)]
    log_format: LogFormat,
    #[arg(long, default_value_t = 0.25)]
    coverage_delta_threshold: f64,
    #[arg(long)]
    diamond_auto_threshold: Option<usize>,
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
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct TaxonomyConfigOverride {
    enabled: Option<bool>,
    min_support: Option<f64>,
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
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct DiamondConfigOverride {
    auto_threshold: Option<usize>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ConsensusConfigOverride {
    min_hits: Option<usize>,
    max_panel: Option<usize>,
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

fn write_run_metrics(out_dir: &str, metrics: &RunMetrics) -> Result<(), Box<dyn std::error::Error>> {
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
            let cfg = resolve_effective_config(&file_cfg, &args)?;

            fs::create_dir_all(&cfg.out)?;

            // Pre-flight: tool versions
            let tools = preflight(
                &cfg.diamond_bin,
                args.mafft_bin.as_deref(),
                args.hmmscan_bin.as_deref(),
            );

            // DIAMOND pre-run
            let dia_cfg = DiamondConfig {
                bin: cfg.diamond_bin.clone(),
                db: cfg.db.clone(),
                query_fasta: cfg.fasta.clone(),
                threads: cfg.threads,
                out_dir: cfg.out.clone(),
                out_name: "diamond.blastp.tsv".to_string(),
                retries: 2,
            };
            let log_json = matches!(args.log_format, LogFormat::Json);
            let t_diamond = step_start("diamond", log_json);
            let diamond_tsv = match args.diamond_mode {
                DiamondMode::Single => blastp_once(&dia_cfg)?,
                DiamondMode::Batch => diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?,
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
            let qlen_map: HashMap<String, usize> = metrics.iter().map(|m| (m.gene_id.clone(), m.length)).collect();
            let grouped = diamond::parse_tsv_grouped(&diamond_tsv, Some(&qlen_map), Some(500)).unwrap_or_default();
            let cons_cfg = consensus::ConsensusConfig {
                min_hits: file_cfg
                    .consensus
                    .as_ref()
                    .and_then(|c| c.min_hits)
                    .unwrap_or(5),
                max_panel: file_cfg
                    .consensus
                    .as_ref()
                    .and_then(|c| c.max_panel)
                    .unwrap_or(10),
                ..Default::default()
            };
            let mut panel_map: HashMap<String, Vec<String>> = HashMap::new();
            for m in &metrics {
                if let Some(hits) = grouped.get(&m.gene_id) {
                    let ids = consensus::select_panel(hits, &cons_cfg);
                    if !ids.is_empty() {
                        panel_map.insert(m.gene_id.clone(), ids);
                    }
                }
            }
            let consensus_secs = step_finish("consensus", t_consensus, log_json);

            // Optional MAFFT alignment metrics
            let mut alignment_map: HashMap<String, AlignmentMetrics> = HashMap::new();
            if let (Some(ref_fasta), Some(mafft_bin)) =
                (cfg.reference_fasta.as_ref(), args.mafft_bin.as_ref())
            {
                let t_mafft = step_start("mafft", log_json);
                let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
                let seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                for g in &metrics {
                    let ids = panel_map.get(&g.gene_id).cloned().unwrap_or_default();
                    if ids.is_empty() {
                        continue;
                    }
                    let mut hits: HashMap<String, Vec<u8>> = HashMap::new();
                    for id in ids {
                        let key = taxonomy::canonical_accession(&id);
                        if let Some(s) = seqs.get(&key) {
                            hits.insert(key.clone(), s.clone());
                        }
                    }
                    if let Some((_im, qseq)) = intrinsic_map.get(&g.gene_id) {
                        if !hits.is_empty() {
                            if let Ok(m) = run_mafft(mafft_bin, &g.gene_id, qseq, &hits) {
                                alignment_map.insert(g.gene_id.clone(), m);
                            }
                        }
                    }
                }
                let _ = step_finish("mafft", t_mafft, log_json);
            }

            // Optional HMMER/Pfam domain summary (JSONL)
            let mut hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut _ref_hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
            let mut domains_arch_map: HashMap<String, f64> = HashMap::new();
            if let (Some(hmm), Some(pfam_db)) = (args.hmmscan_bin.as_ref(), args.pfam_db.as_ref().or(file_cfg.pfam_db.as_ref())) {
                let t_hmmer = step_start("hmmer", log_json);
                let items: Vec<(String, Vec<u8>)> = intrinsic_map.iter().map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone())).collect();
                let top_n = args
                    .hmmer_top_n
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.top_n))
                    .unwrap_or(5);
                let hmmer_threads = args
                    .hmmer_threads
                    .or_else(|| file_cfg.hmmer.as_ref().and_then(|h| h.threads))
                    .unwrap_or(cfg.threads);
                if let Ok(m) = hmmer::run_hmmscan_batch(hmm, pfam_db, items, hmmer_threads, top_n) {
                    hmmsum_map = m;
                }
                if let Some(ref_fasta) = cfg.reference_fasta.as_ref() {
                    let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
                    let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                    if !ref_seqs.is_empty() {
                        let ref_items: Vec<(String, Vec<u8>)> = ref_seqs.into_iter().collect();
                        if let Ok(m) = hmmer::run_hmmscan_batch(hmm, pfam_db, ref_items, hmmer_threads, top_n) {
                            _ref_hmmsum_map = m;
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
                        let qsum_c = if let Some(ref clans) = clan_map { hmmer::collapse_by_clan(qsum, clans) } else { qsum.clone() };
                        let ref_ids = panel_map.get(&g.gene_id).cloned().unwrap_or_default();
                        let mut ref_map_c: HashMap<String, HmmscanSummary> = HashMap::new();
                        for rid in &ref_ids {
                            let key = taxonomy::canonical_accession(rid);
                            if let Some(s) = _ref_hmmsum_map.get(&key) {
                                let val = if let Some(ref clans) = clan_map { hmmer::collapse_by_clan(s, clans) } else { s.clone() };
                                ref_map_c.insert(key.clone(), val);
                            }
                        }
                        let score = hmmer::domains_architecture_score(&qsum_c, &ref_ids, &ref_map_c);
                        domains_arch_map.insert(g.gene_id.clone(), score);
                    }
                }
                let _ = step_finish("hmmer", t_hmmer, log_json);
            }

            // Emit outputs
            let t_emit = step_start("emit_outputs", log_json);
            let qlen_map: HashMap<String, usize> = metrics.iter().map(|m| (m.gene_id.clone(), m.length)).collect();
            let stats = parse_tsv_stats(&diamond_tsv, Some(&qlen_map)).unwrap_or_default();
            let checksums = collect_checksums(&cfg)?;
            let snapshot = build_config_snapshot(&cfg, &args, &file_cfg);
            write_run_manifest(&cfg, &tools, &checksums, &snapshot)?;
            // Optional taxonomy resolution of top hits → taxid/name/lineage per gene
            let mut taxsum_map: HashMap<String, Option<TaxonSummary>> = HashMap::new();
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
                let resolver = TaxonomyResolver::from_sources(
                    cache_path,
                    cfg.reference_fasta.as_deref(),
                    taxdump_dir,
                )
                .map_err(|e| format!("taxonomy setup failed: {}", e))?;
                if let Some(resolver) = resolver {
                    for m in &metrics {
                        let acc = stats.get(&m.gene_id).and_then(|s| s.top_sseqid.clone());
                        if let Some(acc) = acc {
                            taxsum_map.insert(m.gene_id.clone(), resolver.lookup(&acc));
                        } else {
                            taxsum_map.insert(m.gene_id.clone(), None);
                        }
                    }
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
            );
            let comp_map = build_component_scores(
                &metrics,
                &stats,
                &intrinsic_map,
                &file_cfg.scoring,
                args.enable_taxonomy,
            );
            let _ = step_finish("scoring", t_scoring, log_json);
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
            )?;
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
            totals.insert("diamond_mode", serde_json::json!(format!("{:?}", args.diamond_mode)));
            let steps = vec![
                StepDuration { name: "diamond", seconds: diamond_secs },
                StepDuration { name: "ecs", seconds: ecs_secs },
                StepDuration { name: "intrinsic", seconds: intrinsic_secs },
                StepDuration { name: "consensus", seconds: consensus_secs },
                StepDuration { name: "emit_outputs", seconds: emit_secs },
            ];
            let runm = RunMetrics { schema_version: "1.0", steps, totals };
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
    taxsum_map: &HashMap<String, Option<TaxonSummary>>,
    arch_map: Option<&HashMap<String, f64>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("qc_report.jsonl");
    let mut f = File::create(path)?;
    for m in metrics {
        let warnings = if m.hits == 0 {
            vec!["No DIAMOND hits"]
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
        let taxonomy_score = if taxonomy_enabled {
            Some(compute_taxonomy_score(
                taxsum_map
                    .get(&m.gene_id)
                    .and_then(|x| x.as_ref())
                    .is_some(),
            ))
        } else {
            None
        };
        let (w_h, w_i, w_t) = scoring_weights3(scoring);
        let w_d = scoring
            .as_ref()
            .and_then(|s| s.weights.get("domains").cloned())
            .unwrap_or(0.0);
        let domains_arch_score = arch_map.and_then(|am| am.get(&m.gene_id).cloned()).unwrap_or(0.0);
        let denom = (w_h + w_i + w_t + w_d).max(1e-6);
        let final_score = (w_h * homology_score
            + w_i * intrinsic_score
            + w_t * taxonomy_score.unwrap_or(0.0)
            + w_d * domains_arch_score)
            / denom;
        let fusion_split_flag = summary
            .map(|s| s.coverage_delta > cov_delta_thresh)
            .unwrap_or(false);
        let taxonomy = if taxonomy_enabled {
            if let Some(Some(ts)) = taxsum_map.get(&m.gene_id) {
                serde_json::json!({
                    "status": "enabled",
                    "taxid": ts.taxid,
                    "name": ts.name,
                    "lineage": ts.lineage,
                })
            } else {
                serde_json::json!({"status":"enabled"})
            }
        } else {
            serde_json::json!({"status":"disabled"})
        };
        let cur_gene = &m.gene_id;
        let domains = hmmsum_map.get(cur_gene).map(|d| {
            let score = if let Some(ev) = d.top_evalue { let le = if ev>0.0 { -ev.log10() } else { 100.0 }; (le/20.0).clamp(0.0, 1.0) } else { 0.0 };
            let arch_score = arch_map.and_then(|am| am.get(cur_gene)).cloned().unwrap_or(0.0);
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
            let record = serde_json::json!({
                "gene_id": m.gene_id,
                "taxonomy": taxonomy,
                "score_components": {"taxonomy": taxonomy_score, "homology": homology_score, "intrinsic": intrinsic_score, "domains": domains_arch_score},
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
            })),
            "domains": domains,
            "warnings": warnings,
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
    taxsum_map: &HashMap<String, Option<TaxonSummary>>,
    hmmsum_map: &std::collections::HashMap<String, hmmer::HmmscanSummary>,
    csv_verbose: bool,
    comp_map: Option<&HashMap<String, (f64, f64, Option<f64>, f64)>>,
    classify_no_data: bool,
    arch_map: Option<&HashMap<String, f64>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("qc_summary.csv");
    let mut f = File::create(path)?;
    if csv_verbose {
        writeln!(f, "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,homology_score,intrinsic_score,taxonomy_score,domains_score,domains_arch_score,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,start_concordance,start_class,taxonomy_status,warnings")?;
    } else {
        writeln!(f, "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,coverage_delta,coverage_ratio,fusion_split,final_score,classification,mafft_enabled,conserved_fraction,pairwise_identity,sequences_aligned,query_gap_fraction,gap_run_count,max_gap_run,domains_score,domains_arch_score,taxonomy_score,taxonomy_status,warnings")?;
    }
    for m in metrics {
        let warnings = if m.hits == 0 { "No DIAMOND hits" } else { "" };
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
        let (mafft_enabled, conserved, pid, seqs_aln, qgap, gap_runs, max_gap, start_conc, start_class) =
            if let Some(a) = aln {
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
                )
            } else {
                (0, 0.0, 0.0, 0, 0.0, 0, 0, 0.0, String::new())
            };
        let taxonomy_score = if taxonomy_enabled {
            compute_taxonomy_score(
                taxsum_map
                    .get(&m.gene_id)
                    .and_then(|x| x.as_ref())
                    .is_some(),
            )
        } else {
            0.0
        };
        let taxonomy_status = if taxonomy_enabled {
            "enabled"
        } else {
            "disabled"
        };
        let domains_score_field = if let Some(d) = hmmsum_map.get(&m.gene_id) {
            let s = if let Some(ev) = d.top_evalue { let le = if ev>0.0 { -ev.log10() } else { 100.0 }; (le/20.0).clamp(0.0, 1.0) } else { 0.0 };
            format!("{:.4}", s)
        } else {
            String::new()
        };
        let domains_arch_field = arch_map.and_then(|am| am.get(&m.gene_id)).map(|v| format!("{:.4}", v)).unwrap_or_else(|| String::new());
        // Optional NoData classification override
        if classify_no_data && m.hits == 0 {
            let nonzero_domains = hmmsum_map.get(&m.gene_id).map(|d| d.hits_count > 0).unwrap_or(false);
            if !nonzero_domains {
                classif = "NoData".to_string();
            }
        }
        if csv_verbose {
            let (hs, is, ts_opt, _fs) = comp_map.and_then(|cm| cm.get(&m.gene_id).cloned()).unwrap_or((0.0,0.0,None,final_score));
            let ts_str = if taxonomy_enabled { format!("{:.4}", ts_opt.unwrap_or(0.0)) } else { String::new() };
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
                domains_score_field,
                domains_arch_field.clone(),
                mafft_enabled.to_string(),
                format!("{:.3}", conserved),
                pid.to_string(),
                seqs_aln.to_string(),
                format!("{:.3}", qgap),
                gap_runs.to_string(),
                max_gap.to_string(),
                format!("{:.3}", start_conc),
                start_class,
                taxonomy_status.to_string(),
                warnings.to_string(),
            ].join(",");
            writeln!(f, "{}", row)?;
            continue;
        }
        writeln!(
            f,
            "{},{},{},{:.3},{},{:.3},{:.3},{:.3},{:.3},{:.3},{},{:.3},{},{},{:.3},{:.3},{},{:.3},{},{},{},{},{:.4},{},{}",
            m.gene_id,
            m.hits,
            top_hit,
            top_bitscore,
            top_evalue,
            top_qcov,
            top_scov,
            bsd,
            cov_delta,
            cov_ratio,
            fusion,
            final_score,
            classif,
            mafft_enabled,
            conserved,
            pid,
            seqs_aln,
            qgap,
            gap_runs,
            max_gap,
            domains_arch_field,
            domains_score_field,
            taxonomy_score,
            taxonomy_status,
            warnings
        )?;
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

fn scoring_weights3(scoring: &Option<ScoringConfigOverride>) -> (f64, f64, f64) {
    if let Some(cfg) = scoring {
        let wh = *cfg.weights.get("homology").unwrap_or(&0.6);
        let wi = *cfg.weights.get("intrinsic").unwrap_or(&0.4);
        let wt = *cfg.weights.get("taxonomy").unwrap_or(&0.0);
        (wh, wi, wt)
    } else {
        (0.6, 0.4, 0.0)
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

fn top_ids_for_gene(
    gene_id: &str,
    stats: &HashMap<String, diamond::DiamondHitStats>,
    top_n: usize,
) -> Vec<String> {
    stats
        .get(gene_id)
        .and_then(|s| s.top_sseqid.clone())
        .map(|t| vec![t])
        .unwrap_or_default()
        .into_iter()
        .take(top_n)
        .collect()
}

fn stats_top_ids_for(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    top_n: usize,
) -> Vec<String> {
    let mut ids = Vec::new();
    for m in metrics {
        ids.extend(top_ids_for_gene(&m.gene_id, stats, top_n));
    }
    ids.sort();
    ids.dedup();
    ids
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

fn build_scores_map(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
    arch_map: Option<&HashMap<String, f64>>,
) -> HashMap<String, (f64, String)> {
    let (w_h, w_i, w_t) = scoring_weights3(scoring);
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
            compute_taxonomy_score(s.is_some())
        } else {
            0.0
        };
        let w_d = scoring.as_ref().and_then(|sc| sc.weights.get("domains")).cloned().unwrap_or(0.0);
        let d = arch_map.and_then(|am| am.get(&m.gene_id)).cloned().unwrap_or(0.0);
        let denom = (w_h + w_i + w_t + w_d).max(1e-6);
        let score = (w_h*h + w_i*i + w_t*t + w_d*d) / denom;
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
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
) -> HashMap<String, (f64, f64, Option<f64>, f64)> {
    let (w_h, w_i, w_t) = scoring_weights3(scoring);
    let w_d = scoring.as_ref().and_then(|s| s.weights.get("domains").cloned()).unwrap_or(0.0);
    let mut out: HashMap<String, (f64, f64, Option<f64>, f64)> = HashMap::new();
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic.get(&m.gene_id).map(|t| &t.0).unwrap_or(&default_im);
        let h = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let t_opt = if taxonomy_enabled { Some(compute_taxonomy_score(s.is_some())) } else { None };
        let t_val = t_opt.unwrap_or(0.0);
        let denom = (w_h + w_i + w_t + w_d).max(1e-6);
        let final_score = (w_h*h + w_i*i + w_t*t_val) / denom; // domains added elsewhere for JSONL; CSV uses this for classification unless w_d>0
        out.insert(m.gene_id.clone(), (h, i, t_opt, final_score));
    }
    out
}

// ... (rest of file unchanged for brevity in plan)
