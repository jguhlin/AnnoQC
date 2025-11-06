use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone)]
pub struct DiamondConfig {
    pub bin: String,
    pub db: String,
    pub query_fasta: String,
    pub threads: usize,
    pub out_dir: String,
    pub out_name: String,
}

impl DiamondConfig {
    pub fn out_path(&self) -> PathBuf {
        Path::new(&self.out_dir).join(&self.out_name)
    }
}

pub fn version(bin: &str) -> Result<String, String> {
    let output = Command::new(bin)
        .arg("--version")
        .output()
        .map_err(|e| format!("failed to execute '{} --version': {}", bin, e))?;
    if !output.status.success() {
        return Err(format!(
            "'{} --version' exited with status {}",
            bin, output.status
        ));
    }
    let v = String::from_utf8_lossy(&output.stdout).trim().to_string();
    Ok(v)
}

/// Run DIAMOND blastp and write tabular output. Skips if output exists and is non-empty.
pub fn blastp_once(cfg: &DiamondConfig) -> Result<PathBuf, String> {
    let out_path = cfg.out_path();
    if out_path.exists() {
        if let Ok(meta) = fs::metadata(&out_path) {
            if meta.len() > 0 {
                return Ok(out_path);
            }
        }
    }

    fs::create_dir_all(&cfg.out_dir)
        .map_err(|e| format!("create out_dir {}: {}", cfg.out_dir, e))?;

    let outfmt = "6 qseqid sseqid bitscore evalue length qcovhsp scovhsp pident";
    let status = Command::new(&cfg.bin)
        .arg("blastp")
        .arg("--db")
        .arg(&cfg.db)
        .arg("--query")
        .arg(&cfg.query_fasta)
        .arg("--outfmt")
        .arg(outfmt)
        .arg("--threads")
        .arg(cfg.threads.to_string())
        .arg("--max-target-seqs")
        .arg("25")
        .arg("--quiet")
        .arg("--out")
        .arg(&out_path)
        .status()
        .map_err(|e| format!("failed to run diamond blastp: {}", e))?;
    if !status.success() {
        return Err(format!("diamond blastp exited with status {}", status));
    }
    Ok(out_path)
}

/// Run DIAMOND linclust to quickly cluster a reference FASTA. Idempotent via `.done` file.
pub fn linclust(
    bin: &str,
    reference_fasta: &str,
    out_path: &Path,
    approx_id: u32,
    threads: usize,
) -> Result<(), String> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let status = Command::new(bin)
        .arg("linclust")
        .arg("-d")
        .arg(reference_fasta)
        .arg("-o")
        .arg(out_path)
        .arg("--approx-id")
        .arg(approx_id.to_string())
        .arg("--threads")
        .arg(threads.to_string())
        .status()
        .map_err(|e| format!("failed to run diamond linclust: {}", e))?;
    if !status.success() {
        return Err(format!("diamond linclust exited with status {}", status));
    }
    Ok(())
}

/// Run DIAMOND cluster (more sensitive) to refine clusters. Idempotent via `.done` file.
pub fn cluster(
    bin: &str,
    reference_fasta: &str,
    out_path: &Path,
    approx_id: u32,
    threads: usize,
) -> Result<(), String> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let status = Command::new(bin)
        .arg("cluster")
        .arg("-d")
        .arg(reference_fasta)
        .arg("-o")
        .arg(out_path)
        .arg("--approx-id")
        .arg(approx_id.to_string())
        .arg("--threads")
        .arg(threads.to_string())
        .status()
        .map_err(|e| format!("failed to run diamond cluster: {}", e))?;
    if !status.success() {
        return Err(format!("diamond cluster exited with status {}", status));
    }
    Ok(())
}

#[derive(Debug, Clone, Default)]
pub struct DiamondHitStats {
    pub count: usize,
    pub top_sseqid: Option<String>,
    pub top_bitscore: f64,
    pub top_evalue: String,
    pub top_qcov: f64,
    pub top_scov: f64,
    pub top_len: usize,
    pub top_pident: f64,
}

/// Parse diamond tsv produced by `blastp_once` and compute per-query top-hit stats.
pub fn parse_tsv_stats(
    tsv: &Path,
) -> Result<std::collections::HashMap<String, DiamondHitStats>, String> {
    let mut map: std::collections::HashMap<String, DiamondHitStats> = Default::default();
    if !tsv.exists() {
        return Ok(map);
    }
    let file = std::fs::File::open(tsv).map_err(|e| e.to_string())?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        let q = cols[0];
        let sseqid = cols[1].to_string();
        let bitscore = cols[2].parse::<f64>().unwrap_or(0.0);
        let evalue = cols[3].to_string();
        let alen = cols[4].parse::<usize>().unwrap_or(0);
        let qcov = cols[5].parse::<f64>().unwrap_or(0.0);
        let scov = cols[6].parse::<f64>().unwrap_or(0.0);
        let pident = cols[7].parse::<f64>().unwrap_or(0.0);
        let entry = map.entry(q.to_string()).or_default();
        entry.count += 1;
        if bitscore > entry.top_bitscore {
            entry.top_bitscore = bitscore;
            entry.top_sseqid = Some(sseqid);
            entry.top_evalue = evalue;
            entry.top_qcov = qcov;
            entry.top_scov = scov;
            entry.top_len = alen;
            entry.top_pident = pident;
        }
    }
    Ok(map)
}
