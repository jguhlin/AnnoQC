use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

fn format_exit_status(status: &std::process::ExitStatus) -> String {
    if let Some(code) = status.code() {
        return format!("exit_code={}", code);
    }
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        if let Some(sig) = status.signal() {
            return format!("signal={}", sig);
        }
    }
    "unknown_exit_status".to_string()
}

fn log_diamond_command(bin: &str, args: &[String]) {
    log::info!("diamond cmd: {} {}", bin, args.join(" "));
}

fn read_trimmed(path: &str) -> Option<String> {
    std::fs::read_to_string(path)
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
}

fn log_resource_snapshot(context: &str) {
    let pid = std::process::id();
    let mem_max = read_trimmed("/sys/fs/cgroup/memory.max");
    let mem_current = read_trimmed("/sys/fs/cgroup/memory.current");
    let mem_high = read_trimmed("/sys/fs/cgroup/memory.high");
    let cpu_max = read_trimmed("/sys/fs/cgroup/cpu.max");
    let oom_score = read_trimmed("/proc/self/oom_score");
    let oom_adj = read_trimmed("/proc/self/oom_score_adj");
    log::info!(
        "resource snapshot {}: pid={} memory.max={:?} memory.current={:?} memory.high={:?} cpu.max={:?} oom_score={:?} oom_score_adj={:?}",
        context,
        pid,
        mem_max,
        mem_current,
        mem_high,
        cpu_max,
        oom_score,
        oom_adj
    );
}

fn read_proc_status_kb(pid: u32, key: &str) -> Option<u64> {
    let path = format!("/proc/{}/status", pid);
    let text = std::fs::read_to_string(path).ok()?;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix(key) {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if let Some(val) = parts.first() {
                if let Ok(kb) = val.parse::<u64>() {
                    return Some(kb);
                }
            }
        }
    }
    None
}

fn write_running_status(
    status_path: &Path,
    stderr_path: &Path,
    out_path: &Path,
    pid: u32,
    extra: serde_json::Value,
) {
    let out_bytes = std::fs::metadata(out_path).map(|m| m.len()).unwrap_or(0);
    let stderr_bytes = std::fs::metadata(stderr_path).map(|m| m.len()).unwrap_or(0);
    let mem_current = read_trimmed("/sys/fs/cgroup/memory.current");
    let vm_rss_kb = read_proc_status_kb(pid, "VmRSS:");
    let vm_size_kb = read_proc_status_kb(pid, "VmSize:");
    write_status(
        status_path,
        serde_json::json!({
            "state": "running",
            "heartbeat_at": now_unix_seconds(),
            "pid": pid,
            "out_path": out_path.to_string_lossy(),
            "stderr_path": stderr_path.to_string_lossy(),
            "out_bytes": out_bytes,
            "stderr_bytes": stderr_bytes,
            "memory_current": mem_current,
            "vm_rss_kb": vm_rss_kb,
            "vm_size_kb": vm_size_kb,
            "extra": extra,
        }),
    );
}

fn now_unix_seconds() -> f64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs_f64())
        .unwrap_or(0.0)
}

fn status_path_for(out_path: &Path) -> PathBuf {
    out_path.with_extension("status.json")
}

fn stderr_path_for(out_path: &Path) -> PathBuf {
    out_path.with_extension("stderr.log")
}

fn write_status(path: &Path, payload: serde_json::Value) {
    if let Err(e) = std::fs::write(path, payload.to_string()) {
        log::warn!(
            "diamond status write failed: file={} err={}",
            path.display(),
            e
        );
    }
}

fn append_stderr(path: &Path, stderr: &[u8]) {
    if stderr.is_empty() {
        return;
    }
    let mut f = match std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    {
        Ok(f) => f,
        Err(e) => {
            log::warn!(
                "diamond stderr open failed: file={} err={}",
                path.display(),
                e
            );
            return;
        }
    };
    use std::io::Write;
    let _ = writeln!(f, "---- stderr @ {:.3} ----", now_unix_seconds());
    let _ = f.write_all(stderr);
    let _ = writeln!(f);
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub enum HitSource {
    #[default]
    SwissProt,
    RefProt(String),
    #[allow(dead_code)]
    Cluster,
}

#[derive(Debug, Clone)]
pub struct DiamondConfig {
    pub bin: String,
    pub db: String,
    pub query_fasta: String,
    pub threads: usize,
    pub out_dir: String,
    pub out_name: String,
    pub retries: usize,
    pub max_hsps: usize,
}

impl DiamondConfig {
    pub fn out_path(&self) -> PathBuf {
        Path::new(&self.out_dir).join(&self.out_name)
    }
}

// version helper moved to preflight module

/// Run DIAMOND blastp and write tabular output. Skips if output exists and is non-empty.
/// Retries failures with a short backoff, blocking until success or retries exhausted.
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

    let status_path = status_path_for(&out_path);
    let stderr_path = stderr_path_for(&out_path);

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
        "qstart",
        "qend",
        "sstart",
        "send",
        "qlen",
        "slen",
        "staxids",
        "slineages",
    ];
    let mut args: Vec<String> = vec![
        "blastp".to_string(),
        "--db".to_string(),
        cfg.db.clone(),
        "--query".to_string(),
        cfg.query_fasta.clone(),
        "--threads".to_string(),
        cfg.threads.to_string(),
        "--max-target-seqs".to_string(),
        "50".to_string(),
        "--max-hsps".to_string(),
        cfg.max_hsps.to_string(),
        "--sensitive".to_string(),
        "--motif-masking".to_string(),
        "0".to_string(),
        "--quiet".to_string(),
        "--out".to_string(),
        out_path.to_string_lossy().to_string(),
    ];
    args.push("--outfmt".to_string());
    args.extend(outfmt_tokens.iter().map(|s| s.to_string()));
    log::info!("diamond stderr: {}", stderr_path.display());
    write_status(
        &status_path,
        serde_json::json!({
            "state": "start",
            "started_at": now_unix_seconds(),
            "out_path": out_path.to_string_lossy(),
            "stderr_path": stderr_path.to_string_lossy(),
            "db": cfg.db,
            "query_fasta": cfg.query_fasta,
            "threads": cfg.threads,
            "args": args.clone(),
        }),
    );
    let mut attempts = 0usize;
    loop {
        attempts += 1;
        if attempts == 1 {
            log_diamond_command(&cfg.bin, &args);
            log_resource_snapshot("diamond_blastp_start");
        }
        let mut child = Command::new(&cfg.bin)
            .args(&args)
            // Avoid buffering long-running stderr into memory; let it stream to the caller.
            .stdout(Stdio::null())
            .stderr(
                std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&stderr_path)
                    .map(Stdio::from)
                    .map_err(|e| {
                        format!(
                            "failed to open diamond stderr log {}: {}",
                            stderr_path.display(),
                            e
                        )
                    })?,
            )
            .spawn()
            .map_err(|e| format!("failed to spawn diamond blastp: {}", e))?;
        let child_pid = child.id();
        write_running_status(
            &status_path,
            &stderr_path,
            &out_path,
            child_pid,
            serde_json::json!({ "attempt": attempts }),
        );
        log::info!("diamond pid: {}", child_pid);
        let stop = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let stop_clone = stop.clone();
        let status_path_clone = status_path.clone();
        let stderr_path_clone = stderr_path.clone();
        let out_path_clone = out_path.clone();
        let heartbeat = std::thread::spawn(move || {
            while !stop_clone.load(std::sync::atomic::Ordering::Relaxed) {
                write_running_status(
                    &status_path_clone,
                    &stderr_path_clone,
                    &out_path_clone,
                    child_pid,
                    serde_json::json!({}),
                );
                std::thread::sleep(std::time::Duration::from_secs(5));
            }
        });
        let status = child
            .wait()
            .map_err(|e| format!("failed to wait on diamond blastp: {}", e))?;
        stop.store(true, std::sync::atomic::Ordering::Relaxed);
        let _ = heartbeat.join();
        if status.success() {
            break;
        }
        if attempts > cfg.retries.max(1) {
            write_status(
                &status_path,
                serde_json::json!({
                    "state": "failed",
                    "finished_at": now_unix_seconds(),
                    "out_path": out_path.to_string_lossy(),
                    "stderr_path": stderr_path.to_string_lossy(),
                    "exit": format_exit_status(&status),
                }),
            );
            return Err(format!(
                "diamond blastp failed (attempt {}): {} (stderr: {})",
                attempts,
                format_exit_status(&status),
                stderr_path.display()
            ));
        }
        std::thread::sleep(std::time::Duration::from_millis(500 * attempts as u64));
    }
    if let Ok(meta) = fs::metadata(&out_path) {
        if meta.len() == 0 {
            log::warn!(
                "diamond blastp produced an empty output at {}",
                out_path.display()
            );
        }
    }
    let out_bytes = fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0);
    write_status(
        &status_path,
        serde_json::json!({
            "state": "done",
            "finished_at": now_unix_seconds(),
            "out_path": out_path.to_string_lossy(),
            "stderr_path": stderr_path.to_string_lossy(),
            "out_bytes": out_bytes,
        }),
    );
    Ok(out_path)
}

/// Chunked mode: split queries into chunks of `chunk_size` records and append outputs.
/// If `log_json` is true, emits per-chunk JSON progress events.
pub fn blastp_chunked(
    cfg: &DiamondConfig,
    chunk_size: usize,
    log_json: bool,
) -> Result<PathBuf, String> {
    use needletail::parse_fastx_file;
    let out_path = cfg.out_path();
    let status_path = status_path_for(&out_path);
    let stderr_path = stderr_path_for(&out_path);
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
        "qstart",
        "qend",
        "sstart",
        "send",
        "qlen",
        "slen",
        "staxids",
        "slineages",
    ];
    let mut reader = parse_fastx_file(&cfg.query_fasta).map_err(|e| e.to_string())?;
    let mut batch: Vec<(String, Vec<u8>)> = Vec::new();
    let mut tmp_idx = 0usize;
    let start = Instant::now();
    let mut processed = 0usize;
    log::info!("diamond stderr: {}", stderr_path.display());
    write_status(
        &status_path,
        serde_json::json!({
            "state": "start",
            "started_at": now_unix_seconds(),
            "out_path": out_path.to_string_lossy(),
            "stderr_path": stderr_path.to_string_lossy(),
            "db": cfg.db,
            "query_fasta": cfg.query_fasta,
            "threads": cfg.threads,
            "chunk_size": chunk_size,
        }),
    );
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| e.to_string())?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        batch.push((id, rec.seq().to_vec()));
        if batch.len() >= chunk_size {
            run_chunk(
                &batch,
                &mut tmp_idx,
                &out_path,
                &outfmt_tokens,
                &stderr_path,
                cfg,
            )?;
            processed += batch.len();
            write_status(
                &status_path,
                serde_json::json!({
                    "state": "running",
                    "heartbeat_at": now_unix_seconds(),
                    "out_path": out_path.to_string_lossy(),
                    "stderr_path": stderr_path.to_string_lossy(),
                    "out_bytes": std::fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0),
                    "processed": processed,
                    "chunk_index": tmp_idx.saturating_sub(1),
                }),
            );
            if log_json {
                let secs = start.elapsed().as_secs_f64();
                let rate = if secs > 0.0 {
                    processed as f64 / secs
                } else {
                    0.0
                };
                log::info!(
                    "{}",
                    serde_json::json!({
                        "event":"diamond_chunk","chunk_index": tmp_idx-1,
                        "queries": batch.len(),"processed": processed,
                        "seconds": format!("{:.2}", secs),
                        "rate": format!("{:.2}", rate)
                    })
                );
            }
            batch.clear();
        }
    }
    if !batch.is_empty() {
        run_chunk(
            &batch,
            &mut tmp_idx,
            &out_path,
            &outfmt_tokens,
            &stderr_path,
            cfg,
        )?;
        processed += batch.len();
        write_status(
            &status_path,
            serde_json::json!({
                "state": "running",
                "heartbeat_at": now_unix_seconds(),
                "out_path": out_path.to_string_lossy(),
                "stderr_path": stderr_path.to_string_lossy(),
                "out_bytes": std::fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0),
                "processed": processed,
                "chunk_index": tmp_idx.saturating_sub(1),
            }),
        );
        if log_json {
            let secs = start.elapsed().as_secs_f64();
            let rate = if secs > 0.0 {
                processed as f64 / secs
            } else {
                0.0
            };
            log::info!(
                "{}",
                serde_json::json!({
                    "event":"diamond_chunk","chunk_index": tmp_idx-1,
                    "queries": batch.len(),"processed": processed,
                    "seconds": format!("{:.2}", secs),
                    "rate": format!("{:.2}", rate)
                })
            );
        }
    }
    let out_bytes = fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0);
    write_status(
        &status_path,
        serde_json::json!({
            "state": "done",
            "finished_at": now_unix_seconds(),
            "out_path": out_path.to_string_lossy(),
            "stderr_path": stderr_path.to_string_lossy(),
            "out_bytes": out_bytes,
        }),
    );
    Ok(out_path)
}

fn run_chunk(
    batch: &[(String, Vec<u8>)],
    idx: &mut usize,
    out_path: &Path,
    outfmt_tokens: &[&str],
    stderr_path: &Path,
    cfg: &DiamondConfig,
) -> Result<(), String> {
    use std::io::{BufReader, Write};
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
    log::info!(
        "diamond chunk: index={} fasta={} out={}",
        idx.saturating_sub(1),
        tmpfasta.display(),
        tmpout.display()
    );
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
            .arg("50")
            .arg("--max-hsps")
            .arg(cfg.max_hsps.to_string())
            .arg("--sensitive")
            .arg("--motif-masking")
            .arg("0")
            .arg("--quiet")
            .arg("--out")
            .arg(&tmpout)
            .output()
            .map_err(|e| e.to_string())?;
        append_stderr(stderr_path, &output.stderr);
        if output.status.success() {
            break;
        }
        if attempts > cfg.retries.max(1) {
            let mut ctx = String::from_utf8_lossy(&output.stderr).to_string();
            if ctx.len() > 400 {
                ctx.truncate(400);
            }
            return Err(format!(
                "diamond chunk blastp failed after {} attempts: {} stderr='{}' (stderr: {})",
                attempts,
                format_exit_status(&output.status),
                ctx,
                stderr_path.display()
            ));
        }
        std::thread::sleep(std::time::Duration::from_millis(300 * attempts as u64));
    }
    // Append
    let mut out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out_path)
        .map_err(|e| e.to_string())?;
    let chunk = std::fs::File::open(&tmpout).map_err(|e| e.to_string())?;
    let mut chunk_reader = BufReader::new(chunk);
    std::io::copy(&mut chunk_reader, &mut out).map_err(|e| e.to_string())?;
    std::fs::remove_file(tmpfasta).ok();
    std::fs::remove_file(tmpout).ok();
    Ok(())
}

/// Count FASTA records for auto mode selection.
pub fn estimate_query_count(fasta: &str) -> Result<usize, String> {
    use needletail::parse_fastx_file;
    let mut reader = parse_fastx_file(fasta).map_err(|e| e.to_string())?;
    let mut n = 0usize;
    while let Some(r) = reader.next() {
        r.map_err(|e| e.to_string())?;
        n += 1;
    }
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

/// Run DIAMOND recluster to fix clustering errors using a prior clusters file as input.
pub fn recluster(
    bin: &str,
    reference_fasta: &str,
    clusters_in: &Path,
    out_path: &Path,
    approx_id: u32,
    member_cover: u32,
    threads: usize,
) -> Result<(), String> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let status = Command::new(bin)
        .arg("recluster")
        .arg("-d")
        .arg(reference_fasta)
        .arg("--clusters")
        .arg(clusters_in)
        .arg("-o")
        .arg(out_path)
        .arg("--approx-id")
        .arg(approx_id.to_string())
        .arg("--member-cover")
        .arg(member_cover.to_string())
        .arg("--threads")
        .arg(threads.to_string())
        .status()
        .map_err(|e| format!("failed to run diamond recluster: {}", e))?;
    if !status.success() {
        return Err(format!("diamond recluster exited with status {}", status));
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
    #[allow(dead_code)]
    pub qcov: f64,
    #[allow(dead_code)]
    pub scov: f64,
    pub pident: f64,
    pub qstart: usize,
    pub qend: usize,
    pub sstart: usize,
    pub send: usize,
    pub qlen: usize,
    pub slen: usize,
    pub source: HitSource,
    pub staxid: Option<u32>,
    pub lineage: Vec<String>,
}

#[derive(Default)]
struct TsvParseDiag {
    lines: usize,
    skipped_empty: usize,
    skipped_short: usize,
    malformed: usize,
    parse_errors: usize,
}

fn parse_f64(diag: &mut TsvParseDiag, value: &str) -> f64 {
    match value.parse::<f64>() {
        Ok(v) => v,
        Err(_) => {
            diag.parse_errors += 1;
            0.0
        }
    }
}

fn parse_usize(diag: &mut TsvParseDiag, value: &str) -> usize {
    match value.parse::<usize>() {
        Ok(v) => v,
        Err(_) => {
            diag.parse_errors += 1;
            0
        }
    }
}

fn log_tsv_diag(tsv: &Path, diag: &TsvParseDiag) {
    if diag.malformed > 0 || diag.parse_errors > 0 || diag.skipped_short > 0 {
        log::warn!(
            "diamond tsv parse: file={} lines={} malformed={} short={} parse_errors={}",
            tsv.display(),
            diag.lines,
            diag.malformed,
            diag.skipped_short,
            diag.parse_errors
        );
    }
}

/// Parse diamond tsv produced by `blastp_once` and compute per-query top-hit stats.
pub fn parse_tsv_stats(
    tsv: &Path,
    qlen_map: Option<&std::collections::HashMap<String, usize>>,
) -> Result<std::collections::HashMap<String, DiamondHitStats>, String> {
    let mut map: std::collections::HashMap<String, DiamondHitStats> = Default::default();
    let mut diag = TsvParseDiag::default();
    if !tsv.exists() {
        return Ok(map);
    }
    let file = std::fs::File::open(tsv).map_err(|e| e.to_string())?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        diag.lines += 1;
        if line.trim().is_empty() {
            diag.skipped_empty += 1;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            diag.skipped_short += 1;
            continue;
        }
        let q = cols[0];
        let sseqid = cols[1].to_string();
        let (bitscore, evalue, alen, qcov, scov, pident) = if cols.len() >= 16 {
            // requested outfmt columns include qcovhsp/scovhsp/pident and extra metadata fields
            (
                parse_f64(&mut diag, cols[2]),
                cols[3].to_string(),
                parse_usize(&mut diag, cols[4]),
                parse_f64(&mut diag, cols[5]) / 100.0,
                parse_f64(&mut diag, cols[6]) / 100.0,
                parse_f64(&mut diag, cols[7]),
            )
        } else if cols.len() == 8 {
            // minimal custom outfmt: qseqid sseqid bitscore evalue length qcovhsp scovhsp pident
            (
                parse_f64(&mut diag, cols[2]),
                cols[3].to_string(),
                parse_usize(&mut diag, cols[4]),
                parse_f64(&mut diag, cols[5]) / 100.0,
                parse_f64(&mut diag, cols[6]) / 100.0,
                parse_f64(&mut diag, cols[7]),
            )
        } else if cols.len() >= 12 {
            // default BLAST 6 order
            let pident = parse_f64(&mut diag, cols[2]);
            let alen = parse_usize(&mut diag, cols[3]);
            let evalue = cols[10].to_string();
            let bitscore = parse_f64(&mut diag, cols[11]);
            // compute qcov from qstart/qend if we know query length
            let qcov = if let Some(map) = qlen_map {
                if let Some(qlen) = map.get(q) {
                    let qstart = parse_f64(&mut diag, cols[6]);
                    let qend = parse_f64(&mut diag, cols[7]);
                    let span = (qend - qstart).abs() + 1.0;
                    if *qlen > 0 {
                        (span / (*qlen as f64)).clamp(0.0, 1.0)
                    } else {
                        0.0
                    }
                } else {
                    0.0
                }
            } else {
                0.0
            };
            (bitscore, evalue, alen, qcov, 0.0, pident)
        } else {
            diag.malformed += 1;
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
            entry.coverage_delta = if qcov.is_finite() && scov.is_finite() {
                (qcov - scov).abs()
            } else {
                0.0
            };
            entry.coverage_ratio = if scov.is_finite() && scov > 0.0 && qcov.is_finite() {
                qcov / scov
            } else {
                0.0
            };
        }
    }
    log_tsv_diag(tsv, &diag);
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
    let mut diag = TsvParseDiag::default();
    if !tsv.exists() {
        return Ok(map);
    }
    let file = std::fs::File::open(tsv).map_err(|e| e.to_string())?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        diag.lines += 1;
        if line.trim().is_empty() {
            diag.skipped_empty += 1;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            diag.skipped_short += 1;
            continue;
        }
        let q = cols[0].to_string();
        let s = cols[1].to_string();
        let mut staxid: Option<u32> = None;
        let mut lineage: Vec<String> = Vec::new();
        let (bitscore, evalue, alen, qcov, scov, pident, qstart, qend, sstart, send, qlen, slen) =
            if cols.len() >= 16 {
                staxid = parse_taxid(cols[14]);
                lineage = parse_lineage(cols[15]);
                (
                    parse_f64(&mut diag, cols[2]),
                    cols[3].to_string(),
                    parse_usize(&mut diag, cols[4]),
                    parse_f64(&mut diag, cols[5]) / 100.0,
                    parse_f64(&mut diag, cols[6]) / 100.0,
                    parse_f64(&mut diag, cols[7]),
                    parse_usize(&mut diag, cols[8]),
                    parse_usize(&mut diag, cols[9]),
                    parse_usize(&mut diag, cols[10]),
                    parse_usize(&mut diag, cols[11]),
                    parse_usize(&mut diag, cols[12]),
                    parse_usize(&mut diag, cols[13]),
                )
            } else if cols.len() >= 14 {
                (
                    parse_f64(&mut diag, cols[2]),
                    cols[3].to_string(),
                    parse_usize(&mut diag, cols[4]),
                    parse_f64(&mut diag, cols[5]) / 100.0,
                    parse_f64(&mut diag, cols[6]) / 100.0,
                    parse_f64(&mut diag, cols[7]),
                    parse_usize(&mut diag, cols[8]),
                    parse_usize(&mut diag, cols[9]),
                    parse_usize(&mut diag, cols[10]),
                    parse_usize(&mut diag, cols[11]),
                    parse_usize(&mut diag, cols[12]),
                    parse_usize(&mut diag, cols[13]),
                )
            } else if cols.len() >= 12 {
                let pident = parse_f64(&mut diag, cols[2]);
                let alen = parse_usize(&mut diag, cols[3]);
                let evalue = cols[10].to_string();
                let bitscore = parse_f64(&mut diag, cols[11]);
                let qcov = if let Some(map) = qlen_map {
                    if let Some(qlen) = map.get(&q) {
                        let qstart = parse_f64(&mut diag, cols[6]);
                        let qend = parse_f64(&mut diag, cols[7]);
                        let span = (qend - qstart).abs() + 1.0;
                        if *qlen > 0 {
                            (span / (*qlen as f64)).clamp(0.0, 1.0)
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                (
                    bitscore,
                    evalue,
                    alen,
                    qcov,
                    0.0,
                    pident,
                    parse_usize(&mut diag, cols[6]),
                    parse_usize(&mut diag, cols[7]),
                    parse_usize(&mut diag, cols[8]),
                    parse_usize(&mut diag, cols[9]),
                    qlen_map.and_then(|m| m.get(&q).cloned()).unwrap_or(0),
                    0,
                )
            } else {
                diag.malformed += 1;
                continue;
            };
        let row = DiamondHitRow {
            qseqid: q.clone(),
            sseqid: s,
            bitscore,
            evalue,
            length: alen,
            qcov,
            scov,
            pident,
            qstart,
            qend,
            sstart,
            send,
            qlen,
            slen,
            source: HitSource::SwissProt,
            staxid,
            lineage,
        };
        let entry = map.entry(q).or_default();
        entry.push(row);
        if let Some(k) = max_per_query {
            if entry.len() >= k {
                continue;
            }
        }
    }
    // Sort each vector by decreasing bitscore
    for v in map.values_mut() {
        v.sort_by(|a, b| {
            b.bitscore
                .partial_cmp(&a.bitscore)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    log_tsv_diag(tsv, &diag);
    Ok(map)
}

fn parse_taxid(raw: &str) -> Option<u32> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let parsed = trimmed
        .split([';', ',', '|', ' '])
        .find_map(|tok| tok.trim().parse::<u32>().ok());
    if parsed.is_none() {
        log::debug!("diamond parse_taxid: no valid taxid in '{}'", trimmed);
    }
    parsed
}

fn parse_lineage(raw: &str) -> Vec<String> {
    if raw.trim().is_empty() {
        return Vec::new();
    }
    raw.split([';', '|', ','])
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}
