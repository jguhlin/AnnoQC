use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use clap::{Args, Parser, Subcommand, ValueEnum};
use mimalloc::MiMalloc;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};

mod checkpoint;
mod diamond;
mod ecs;
mod mafft;
mod metrics;
mod preflight;
mod provenance;
mod taxonomy;
use diamond::{
    blastp_once, cluster as diamond_cluster, linclust as diamond_linclust, parse_tsv_stats,
    version as diamond_version, DiamondConfig,
};
use ecs::{run_scheduler, EcsConfig, GeneMetrics};
use mafft::{load_sequences_by_ids, run_mafft, AlignmentMetrics};
use metrics::compute_intrinsic;
use preflight::preflight;
use provenance::filehash_xx64;
use taxonomy::TaxonomyResolver;

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
    Analyze(AnalyzeArgs),
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
    #[arg(long, value_enum, default_value_t = DiamondMode::Auto)]
    diamond_mode: DiamondMode,
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
}

#[derive(Args, Debug, Clone, Serialize, Deserialize)]
struct TaxonomyCacheArgs {
    #[arg(long)]
    input: String,
    #[arg(long)]
    output: String,
}

#[derive(Debug, Clone, Default, Deserialize)]
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
    diamond_mode: Option<DiamondMode>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct TaxonomyConfigOverride {
    enabled: Option<bool>,
    min_support: Option<f64>,
    profile_db: Option<String>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
struct ScoringConfigOverride {
    #[serde(default)]
    weights: HashMap<String, f64>,
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Prepare(p) => {
            prepare_cmd(p)?;
            Ok(())
        }
        Commands::Analyze(args) => {
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
            };
            let diamond_tsv = blastp_once(&dia_cfg)?;

            // ECS: compute per-gene metrics from FASTA and DIAMOND
            let metrics: Vec<GeneMetrics> = run_scheduler(EcsConfig {
                fasta_path: cfg.fasta.clone(),
                diamond_tsv: diamond_tsv.to_string_lossy().to_string(),
                threads: cfg.threads,
            });

            // Intrinsic metrics per gene
            let intrinsic_map = compute_intrinsic_for_ids(&cfg.fasta, &metrics)?;

            // Optional MAFFT alignment metrics
            let mut alignment_map: HashMap<String, AlignmentMetrics> = HashMap::new();
            if let (Some(ref_fasta), Some(mafft_bin)) =
                (cfg.reference_fasta.as_ref(), args.mafft_bin.as_ref())
            {
                let stats_all = parse_tsv_stats(&diamond_tsv).unwrap_or_default();
                let all_ids = stats_top_ids_for(&metrics, &stats_all, args.alignment_top_hits);
                let seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
                for g in &metrics {
                    let top_ids = top_ids_for_gene(&g.gene_id, &stats_all, args.alignment_top_hits);
                    if top_ids.is_empty() {
                        continue;
                    }
                    let mut hits: HashMap<String, Vec<u8>> = HashMap::new();
                    for id in top_ids {
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
            }

            // Emit outputs
            let stats = parse_tsv_stats(&diamond_tsv).unwrap_or_default();
            let checksums = collect_checksums(&cfg)?;
            write_run_manifest(&cfg, &tools, &checksums)?;
            write_jsonl_metrics(&cfg.out, &metrics, &stats, &intrinsic_map, &alignment_map)?;
            write_csv_metrics(&cfg.out, &metrics, &stats)?;
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

fn write_run_manifest(
    cfg: &EffectiveConfig,
    tools: &preflight::ToolVersions,
    sums: &Checksums,
) -> Result<(), Box<dyn std::error::Error>> {
    #[derive(Serialize)]
    struct Manifest<'a> {
        tool: &'a str,
        diamond_version: &'a str,
        mafft_version: Option<&'a str>,
        hmmscan_version: Option<&'a str>,
        fasta: &'a str,
        db: &'a str,
        threads: usize,
        out: &'a str,
        fasta_xx64: Option<&'a str>,
        db_xx64: Option<&'a str>,
    }
    let manifest = Manifest {
        tool: "AnnoQC",
        diamond_version: tools.diamond.as_deref().unwrap_or_default(),
        mafft_version: tools.mafft.as_deref(),
        hmmscan_version: tools.hmmscan.as_deref(),
        fasta: &cfg.fasta,
        db: &cfg.db,
        threads: cfg.threads,
        out: &cfg.out,
        fasta_xx64: sums.fasta_xx64.as_deref(),
        db_xx64: sums.db_xx64.as_deref(),
    };
    let path = Path::new(&cfg.out).join("run.json");
    let text = serde_json::to_string_pretty(&manifest)?;
    fs::write(path, text)?;
    Ok(())
}

fn write_jsonl_metrics(
    out_dir: &str,
    metrics: &[ecs::GeneMetrics],
    stats: &std::collections::HashMap<String, diamond::DiamondHitStats>,
    intrinsic_map: &std::collections::HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    alignment_map: &std::collections::HashMap<String, mafft::AlignmentMetrics>,
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
        let record = serde_json::json!({
            "gene_id": m.gene_id,
            "taxonomy": {"status": "disabled"},
            "score_components": {"taxonomy": serde_json::Value::Null},
            "homology": {
                "hits_count": m.hits,
                "top_hit": summary.and_then(|s| s.top_sseqid.clone()),
                "top_bitscore": summary.map(|s| s.top_bitscore),
                "top_evalue": summary.as_ref().map(|s| s.top_evalue.clone()),
                "top_qcov": summary.map(|s| s.top_qcov),
                "top_scov": summary.map(|s| s.top_scov),
                "bitscore_density": summary.map(|s| if s.top_len>0 { s.top_bitscore / s.top_len as f64 } else { 0.0 }),
            },
            "intrinsic": {
                "ambiguous_fraction": intrinsic.ambiguous_fraction,
                "max_homopolymer": intrinsic.max_homopolymer,
                "low_complexity_fraction": intrinsic.low_complexity_fraction,
                "low_complexity_windows": intrinsic.low_complexity_windows,
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
            "warnings": warnings,
        });
        writeln!(f, "{}", serde_json::to_string(&record)?)?;
    }
    Ok(())
}

fn write_csv_metrics(
    out_dir: &str,
    metrics: &[ecs::GeneMetrics],
    stats: &std::collections::HashMap<String, diamond::DiamondHitStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("qc_summary.csv");
    let mut f = File::create(path)?;
    writeln!(f, "gene_id,hits_count,top_hit,top_bitscore,top_evalue,top_qcov,top_scov,bitscore_density,taxonomy_score,taxonomy_status,warnings")?;
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
        writeln!(
            f,
            "{},{},{},{:.3},{},{:.3},{:.3},{:.3},{},{},{}",
            m.gene_id,
            m.hits,
            top_hit,
            top_bitscore,
            top_evalue,
            top_qcov,
            top_scov,
            bsd,
            "",
            "disabled",
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
    checkpoint::run_step(db_done, p.resume, "diamond makedb", || {
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
    checkpoint::run_step(clusters_done, p.resume, "diamond linclust", || {
        diamond_linclust(&diamond_bin, &p.fasta, clusters, p.approx_id, p.threads)
            .map_err(|e| format!("linclust failed: {}", e))
    })?;

    // Step 4: cluster (sensitive) → clusters.realign
    let realign = Path::new("clusters.realign");
    let realign_done = Path::new("clusters.realign.done");
    checkpoint::run_step(realign_done, p.resume, "diamond cluster", || {
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

// ... (rest of file unchanged for brevity in plan)
