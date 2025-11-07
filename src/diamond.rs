use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct DiamondConfig {
    pub bin: String,
    pub db: String,
    pub query_fasta: String,
    pub threads: usize,
    pub out_dir: String,
    pub out_name: String,
    pub retries: usize,
}

impl DiamondConfig {
    pub fn out_path(&self) -> PathBuf {
        Path::new(&self.out_dir).join(&self.out_name)
    }
}

// version helper moved to preflight module

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

    // Request explicit outfmt 6 columns by passing tokens separately.
    let outfmt_tokens = [
        "6",
        "qseqid",
        "sseqid",
        "bitscore",
        "evalue",
        "length",
        "qcovhsp",
        "scovhsp",
        "pident",
    ];
    let mut attempts = 0usize;
    loop {
        attempts += 1;
        let mut cmd = Command::new(&cfg.bin);
        let output = cmd
            .arg("blastp")
            .arg("--db")
            .arg(&cfg.db)
            .arg("--query")
            .arg(&cfg.query_fasta)
            .arg("--outfmt")
            .args(outfmt_tokens)
            .arg("--threads")
            .arg(cfg.threads.to_string())
            .arg("--max-target-seqs")
            .arg("25")
            .arg("--quiet")
            .arg("--out")
            .arg(&out_path)
            .output()
            .map_err(|e| format!("failed to run diamond blastp: {}", e))?;
        if output.status.success() {
            break;
        }
        if attempts > cfg.retries.max(1) {
            let mut ctx = String::from_utf8_lossy(&output.stderr).to_string();
            if ctx.len() > 400 { ctx.truncate(400); }
            return Err(format!("diamond blastp failed (attempt {}): status={} stderr='{}'", attempts, output.status, ctx));
        }
        std::thread::sleep(std::time::Duration::from_millis(500 * attempts as u64));
    }
    Ok(out_path)
}

/// Chunked mode: split queries into chunks of `chunk_size` records and append outputs.
/// If `log_json` is true, emits per-chunk JSON progress events.
pub fn blastp_chunked(cfg: &DiamondConfig, chunk_size: usize, log_json: bool) -> Result<PathBuf, String> {
    use needletail::parse_fastx_file;
    let out_path = cfg.out_path();
    if out_path.exists() {
        std::fs::remove_file(&out_path).ok();
    }
    fs::create_dir_all(&cfg.out_dir).map_err(|e| e.to_string())?;
    let outfmt_tokens = [
        "6",
        "qseqid",
        "sseqid",
        "bitscore",
        "evalue",
        "length",
        "qcovhsp",
        "scovhsp",
        "pident",
    ];
    let mut reader = parse_fastx_file(&cfg.query_fasta).map_err(|e| e.to_string())?;
    let mut batch: Vec<(String, Vec<u8>)> = Vec::new();
    let mut tmp_idx = 0usize;
    let start = Instant::now();
    let mut processed = 0usize;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| e.to_string())?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        batch.push((id, rec.seq().to_vec()));
        if batch.len() >= chunk_size {
            run_chunk(&batch, &mut tmp_idx, &out_path, &outfmt_tokens, cfg)?;
            processed += batch.len();
            if log_json {
                let secs = start.elapsed().as_secs_f64();
                let rate = if secs > 0.0 { processed as f64 / secs } else { 0.0 };
                log::info!("{}", serde_json::json!({
                    "event":"diamond_chunk","chunk_index": tmp_idx-1,
                    "queries": batch.len(),"processed": processed,
                    "seconds": format!("{:.2}", secs),
                    "rate": format!("{:.2}", rate)
                }));
            }
            batch.clear();
        }
    }
    if !batch.is_empty() {
        run_chunk(&batch, &mut tmp_idx, &out_path, &outfmt_tokens, cfg)?;
        processed += batch.len();
        if log_json {
            let secs = start.elapsed().as_secs_f64();
            let rate = if secs > 0.0 { processed as f64 / secs } else { 0.0 };
            log::info!("{}", serde_json::json!({
                "event":"diamond_chunk","chunk_index": tmp_idx-1,
                "queries": batch.len(),"processed": processed,
                "seconds": format!("{:.2}", secs),
                "rate": format!("{:.2}", rate)
            }));
        }
    }
    Ok(out_path)
}

fn run_chunk(
    batch: &[(String, Vec<u8>)],
    idx: &mut usize,
    out_path: &Path,
    outfmt_tokens: &[&str],
    cfg: &DiamondConfig,
) -> Result<(), String> {
    use std::io::Write;
    let tmpfasta = out_path.with_extension(format!("chunk{}.fa", *idx));
    *idx += 1;
    {
        let mut f = std::fs::File::create(&tmpfasta).map_err(|e| e.to_string())?;
        for (id, seq) in batch {
            writeln!(f, ">{}", id).map_err(|e| e.to_string())?;
            writeln!(f, "{}", String::from_utf8_lossy(seq)).map_err(|e| e.to_string())?;
        }
    }
    let tmpout = out_path.with_extension(format!("chunk{}.tsv", *idx));
    let mut attempts = 0usize;
    loop {
        attempts += 1;
        let mut cmd = Command::new(&cfg.bin);
        let output = cmd
            .arg("blastp")
            .arg("--db")
            .arg(&cfg.db)
            .arg("--query")
            .arg(&tmpfasta)
            .arg("--outfmt")
            .args(outfmt_tokens)
            .arg("--threads")
            .arg("1")
            .arg("--max-target-seqs")
            .arg("25")
            .arg("--quiet")
            .arg("--out")
            .arg(&tmpout)
            .output()
            .map_err(|e| e.to_string())?;
        if output.status.success() {
            break;
        }
        if attempts > cfg.retries.max(1) {
            let mut ctx = String::from_utf8_lossy(&output.stderr).to_string();
            if ctx.len() > 400 { ctx.truncate(400); }
            return Err(format!("diamond chunk blastp failed after {} attempts: status={} stderr='{}'", attempts, output.status, ctx));
        }
        std::thread::sleep(std::time::Duration::from_millis(300 * attempts as u64));
    }
    // Append
    let mut out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out_path)
        .map_err(|e| e.to_string())?;
    let chunk = std::fs::read(&tmpout).map_err(|e| e.to_string())?;
    out.write_all(&chunk).map_err(|e| e.to_string())?;
    std::fs::remove_file(tmpfasta).ok();
    std::fs::remove_file(tmpout).ok();
    Ok(())
}

/// Count FASTA records for auto mode selection.
pub fn estimate_query_count(fasta: &str) -> Result<usize, String> {
    use needletail::parse_fastx_file;
    let mut reader = parse_fastx_file(fasta).map_err(|e| e.to_string())?;
    let mut n = 0usize;
    while let Some(r) = reader.next() { r.map_err(|e| e.to_string())?; n += 1; }
    Ok(n)
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
    pub coverage_delta: f64,
    pub coverage_ratio: f64,
}

#[derive(Debug, Clone, Default)]
pub struct DiamondHitRow {
    pub qseqid: String,
    pub sseqid: String,
    pub bitscore: f64,
    pub evalue: String,
    pub length: usize,
    pub qcov: f64,
    pub scov: f64,
    pub pident: f64,
}

/// Parse diamond tsv produced by `blastp_once` and compute per-query top-hit stats.
pub fn parse_tsv_stats(
    tsv: &Path,
    qlen_map: Option<&std::collections::HashMap<String, usize>>,
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
        let q = cols[0];
        let sseqid = cols[1].to_string();
        let (bitscore, evalue, alen, qcov, scov, pident) = if cols.len() >= 8 && cols.len() < 12 {
            // our requested outfmt: 6 qseqid sseqid bitscore evalue length qcovhsp scovhsp pident
            (
                cols[2].parse::<f64>().unwrap_or(0.0),
                cols[3].to_string(),
                cols[4].parse::<usize>().unwrap_or(0),
                cols[5].parse::<f64>().unwrap_or(0.0) / 100.0,
                cols[6].parse::<f64>().unwrap_or(0.0) / 100.0,
                cols[7].parse::<f64>().unwrap_or(0.0),
            )
        } else if cols.len() >= 12 {
            // default BLAST 6 order
            let pident = cols[2].parse::<f64>().unwrap_or(0.0);
            let alen = cols[3].parse::<usize>().unwrap_or(0);
            let evalue = cols[10].to_string();
            let bitscore = cols[11].parse::<f64>().unwrap_or(0.0);
            // compute qcov from qstart/qend if we know query length
            let qcov = if let Some(map) = qlen_map {
                if let Some(qlen) = map.get(q) {
                    let qstart = cols[6].parse::<f64>().unwrap_or(0.0);
                    let qend = cols[7].parse::<f64>().unwrap_or(0.0);
                    let span = (qend - qstart).abs() + 1.0;
                    if *qlen > 0 { (span / (*qlen as f64)).clamp(0.0, 1.0) } else { 0.0 }
                } else { 0.0 }
            } else { 0.0 };
            (bitscore, evalue, alen, qcov, 0.0, pident)
        } else {
            continue;
        };
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
            entry.coverage_delta = (qcov - scov).abs();
            entry.coverage_ratio = if scov > 0.0 { qcov / scov } else { 0.0 };
        }
    }
    Ok(map)
}

/// Parse diamond tsv into grouped hit rows per query. Best-effort parsing of either our
/// explicit outfmt or default BLAST 6, using optional qlen_map to compute qcov.
pub fn parse_tsv_grouped(
    tsv: &Path,
    qlen_map: Option<&std::collections::HashMap<String, usize>>,
    max_per_query: Option<usize>,
) -> Result<std::collections::HashMap<String, Vec<DiamondHitRow>>, String> {
    let mut map: std::collections::HashMap<String, Vec<DiamondHitRow>> = Default::default();
    if !tsv.exists() { return Ok(map); }
    let file = std::fs::File::open(tsv).map_err(|e| e.to_string())?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 { continue; }
        let q = cols[0].to_string();
        let s = cols[1].to_string();
        let (bitscore, evalue, alen, qcov, scov, pident) = if cols.len() >= 8 && cols.len() < 12 {
            (
                cols[2].parse::<f64>().unwrap_or(0.0),
                cols[3].to_string(),
                cols[4].parse::<usize>().unwrap_or(0),
                cols[5].parse::<f64>().unwrap_or(0.0) / 100.0,
                cols[6].parse::<f64>().unwrap_or(0.0) / 100.0,
                cols[7].parse::<f64>().unwrap_or(0.0),
            )
        } else if cols.len() >= 12 {
            let pident = cols[2].parse::<f64>().unwrap_or(0.0);
            let alen = cols[3].parse::<usize>().unwrap_or(0);
            let evalue = cols[10].to_string();
            let bitscore = cols[11].parse::<f64>().unwrap_or(0.0);
            let qcov = if let Some(map) = qlen_map {
                if let Some(qlen) = map.get(&q) {
                    let qstart = cols[6].parse::<f64>().unwrap_or(0.0);
                    let qend = cols[7].parse::<f64>().unwrap_or(0.0);
                    let span = (qend - qstart).abs() + 1.0;
                    if *qlen > 0 { (span / (*qlen as f64)).clamp(0.0, 1.0) } else { 0.0 }
                } else { 0.0 }
            } else { 0.0 };
            (bitscore, evalue, alen, qcov, 0.0, pident)
        } else { continue; };
        let row = DiamondHitRow { qseqid: q.clone(), sseqid: s, bitscore, evalue, length: alen, qcov, scov, pident };
        let entry = map.entry(q).or_default();
        entry.push(row);
        if let Some(k) = max_per_query { if entry.len() >= k { continue; } }
    }
    // Sort each vector by decreasing bitscore
    for v in map.values_mut() {
        v.sort_by(|a,b| b.bitscore.partial_cmp(&a.bitscore).unwrap_or(std::cmp::Ordering::Equal));
    }
    Ok(map)
}
