use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::Serialize;

use crate::{
    filehash_xx64, preflight, AnalyzeArgs, CalibrationSettings, Checksums, EffectiveConfig,
    FileConfig, ReportFormat, OUTPUT_SCHEMA_VERSION, PLUGIN_SCHEMA_VERSION,
};

#[derive(Serialize, Clone)]
pub(crate) struct PluginManifest {
    name: String,
    path: String,
    xx64: Option<String>,
    size_bytes: Option<u64>,
    version: Option<String>,
}

#[derive(Serialize, Clone)]
pub(crate) struct RuleManifest {
    path: String,
    xx64: Option<String>,
}

pub(crate) fn collect_plugin_manifest(paths: &[String]) -> Vec<PluginManifest> {
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
                if let Ok(text) = fs::read_to_string(&sidecar) {
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

pub(crate) fn collect_rule_manifest(paths: &[String]) -> Vec<RuleManifest> {
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

pub(crate) fn load_rendered_gene_ids(out_dir: &str) -> std::collections::HashSet<String> {
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
pub(crate) struct ConfigSnapshot<'a> {
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
    scoring_weights: HashMap<String, f64>,
    calibration_mode: String,
    calibration_min_samples: usize,
    calibration_min_unique: usize,
}

pub(crate) fn write_run_manifest(
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

pub(crate) fn build_config_snapshot<'a>(
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
