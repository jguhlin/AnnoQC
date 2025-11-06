use std::io::BufRead;
use std::path::Path;
use std::process::{Command, Stdio};

#[derive(Debug, Clone, Default)]
pub struct HmmscanHit {
    pub target_name: String,
    pub accession: String,
    pub evalue: f64,
    pub score: f64,
    pub bias: f64,
}

#[derive(Debug, Clone, Default)]
pub struct HmmscanSummary {
    pub hits_count: usize,
    pub top_accession: Option<String>,
    pub top_evalue: Option<f64>,
}

pub fn run_hmmscan(bin: &str, db_path: &str, id: &str, seq: &[u8]) -> Result<HmmscanSummary, String> {
    let mut child = Command::new(bin)
        .arg("--domtblout").arg("/dev/stdout")
        .arg(db_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("failed to start hmmscan: {}", e))?;
    {
        use std::io::Write;
        let stdin = child.stdin.as_mut().ok_or("hmmscan stdin unavailable")?;
        writeln!(stdin, ">{}", id).map_err(|e| e.to_string())?;
        stdin.write_all(seq).map_err(|e| e.to_string())?;
        writeln!(stdin).map_err(|e| e.to_string())?;
    }
    let out = child.wait_with_output().map_err(|e| e.to_string())?;
    if !out.status.success() { return Err(format!("hmmscan exited with status {}", out.status)); }
    parse_domtblout(&out.stdout[..])
}

pub fn parse_domtblout(bytes: &[u8]) -> Result<HmmscanSummary, String> {
    let mut hits: Vec<HmmscanHit> = Vec::new();
    for line in std::io::BufReader::new(bytes).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim_start().starts_with('#') || line.trim().is_empty() { continue; }
        // domtblout columns: target name, accession, tlen, query name, accession, qlen, ... , i-Evalue, score, bias, ...
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 13 { continue; }
        let target_name = cols[0].to_string();
        let accession = cols[1].to_string();
        let i_eval = cols[12].parse::<f64>().unwrap_or(1.0);
        let score = cols.get(13).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        let bias = cols.get(14).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
        hits.push(HmmscanHit { target_name, accession, evalue: i_eval, score, bias });
    }
    hits.sort_by(|a,b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let top = hits.first();
    Ok(HmmscanSummary {
        hits_count: hits.len(),
        top_accession: top.map(|h| h.accession.clone()),
        top_evalue: top.map(|h| h.evalue),
    })
}

