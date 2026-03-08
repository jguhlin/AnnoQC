use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

#[derive(Debug, Clone, Serialize)]
pub struct StepDuration<'a> {
    pub name: &'a str,
    pub seconds: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunMetrics<'a> {
    pub schema_version: &'a str,
    pub steps: Vec<StepDuration<'a>>,
    pub totals: HashMap<&'a str, serde_json::Value>,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunCounts {
    pub total_genes: usize,
    pub rendered_genes: usize,
    pub high: usize,
    pub low: usize,
    pub x: usize,
    pub high_complete: usize,
    pub high_fragmented: usize,
    pub low_novel: usize,
    pub low_artifact: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunThroughput {
    pub total_seconds: f64,
    pub genes_per_second: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunFeatures {
    pub aligner: String,
    pub report_format: String,
    pub resume: bool,
    pub taxonomy_enabled: bool,
    pub hmmer_enabled: bool,
    pub alignment_enabled: bool,
    pub genomic_enabled: bool,
    pub plugins: usize,
    pub rules: usize,
    pub nucleotide: bool,
    pub calibration: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct RunSummary<'a> {
    pub schema_version: &'a str,
    pub counts: RunCounts,
    pub throughput: RunThroughput,
    pub features: RunFeatures,
    pub timings: Vec<StepDuration<'a>>,
    pub errors: Vec<String>,
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

pub fn step_start(name: &str, json_logs: bool) -> Instant {
    if json_logs {
        log::info!("{}", serde_json::json!({"event":"step_start","name":name}));
    } else {
        log::info!("step_start: {}", name);
    }
    Instant::now()
}

pub fn step_finish(name: &str, start: Instant, json_logs: bool) -> f64 {
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

pub fn write_run_metrics(
    out_dir: &str,
    metrics: &RunMetrics,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("run_metrics.json");
    let text = serde_json::to_string_pretty(metrics)?;
    std::fs::write(path, text)?;
    Ok(())
}

pub fn write_run_summary(
    out_dir: &str,
    summary: &RunSummary,
) -> Result<(), Box<dyn std::error::Error>> {
    let path = Path::new(out_dir).join("run_summary.json");
    let text = serde_json::to_string_pretty(summary)?;
    std::fs::write(path, text)?;
    Ok(())
}

pub fn write_slowest_genes(
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
    let path = Path::new(out_dir).join("slowest_genes.json");
    let text = serde_json::to_string_pretty(&payload)?;
    std::fs::write(path, text)?;
    Ok(())
}
