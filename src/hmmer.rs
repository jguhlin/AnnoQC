use std::io::BufRead;
use std::sync::{Arc, Mutex};
use std::thread;
use std::process::{Command, Stdio};

#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
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
    pub hits: Vec<HmmscanHit>,
}

pub fn run_hmmscan(
    bin: &str,
    db_path: &str,
    id: &str,
    seq: &[u8],
) -> Result<HmmscanSummary, String> {
    let mut child = Command::new(bin)
        .arg("--domtblout")
        .arg("/dev/stdout")
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
    if !out.status.success() {
        return Err(format!("hmmscan exited with status {}", out.status));
    }
    parse_domtblout(&out.stdout[..])
}

pub fn parse_domtblout(bytes: &[u8]) -> Result<HmmscanSummary, String> {
    let mut hits: Vec<HmmscanHit> = Vec::new();
    for line in std::io::BufReader::new(bytes).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim_start().starts_with('#') || line.trim().is_empty() {
            continue;
        }
        // domtblout columns: target name, accession, tlen, query name, accession, qlen, ... , i-Evalue, score, bias, ...
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 13 {
            continue;
        }
        let target_name = cols[0].to_string();
        let accession = cols[1].to_string();
        let i_eval = cols[12].parse::<f64>().unwrap_or(1.0);
        let score = cols
            .get(13)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let bias = cols
            .get(14)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        hits.push(HmmscanHit {
            target_name,
            accession,
            evalue: i_eval,
            score,
            bias,
        });
    }
    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let top = hits.first();
    Ok(HmmscanSummary {
        hits_count: hits.len(),
        top_accession: top.map(|h| h.accession.clone()),
        top_evalue: top.map(|h| h.evalue),
        hits,
    })
}

#[cfg(test)]
mod tests {
    use super::parse_domtblout;

    #[test]
    fn parse_simple_domtblout() {
        // Minimal domtblout snippet with one hit line
        let data = b"# target name  acc  tlen  query name  acc  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to  from  to  acc  description\n\
PF00001.1  PF00001.1  250  Q12345  -  120  1e-20  100.0  0.1  1  1  2e-20  1e-20  99.0  0.0  1  100  1  100  1  100  0.95  Some description\n";
        let sum = parse_domtblout(data).expect("parse");
        assert_eq!(sum.hits_count, 1);
        assert_eq!(sum.top_accession.as_deref(), Some("PF00001.1"));
        assert!(sum.top_evalue.unwrap() <= 1e-20 * 1.0001);
        assert_eq!(sum.hits.len(), 1);
        let h = &sum.hits[0];
        assert_eq!(h.accession, "PF00001.1");
    }
}

pub fn run_hmmscan_batch(
    bin: &str,
    db_path: &str,
    items: Vec<(String, Vec<u8>)>,
    threads: usize,
    top_n: usize,
) -> Result<std::collections::HashMap<String, HmmscanSummary>, String> {
    let nthreads = threads.max(1);
    let queue = Arc::new(Mutex::new(items.into_iter()));
    let results: Arc<Mutex<std::collections::HashMap<String, HmmscanSummary>>> =
        Arc::new(Mutex::new(Default::default()));
    let mut handles = Vec::new();
    for _ in 0..nthreads {
        let q = Arc::clone(&queue);
        let r = Arc::clone(&results);
        let bin_s = bin.to_string();
        let db_s = db_path.to_string();
        let handle = thread::spawn(move || {
            loop {
                let next = {
                    let mut guard = q.lock().unwrap();
                    guard.next()
                };
                let Some((gid, seq)) = next else { break; };
                let sum = run_hmmscan(&bin_s, &db_s, &gid, &seq)
                    .unwrap_or_default();
                let mut trimmed = sum.clone();
                if trimmed.hits.len() > top_n {
                    trimmed.hits.truncate(top_n);
                }
                let mut out = r.lock().unwrap();
                out.insert(gid, trimmed);
            }
        });
        handles.push(handle);
    }
    for h in handles { h.join().map_err(|_| "hmmscan thread panicked".to_string())?; }
    let map = Arc::try_unwrap(results).map_err(|_| "results arc busy".to_string())
        .and_then(|m| m.into_inner().map_err(|_| "results poisoned".to_string()))?;
    Ok(map)
}

#[allow(dead_code)]
pub fn domains_architecture_score(
    query: &HmmscanSummary,
    ref_ids: &[String],
    ref_map: &std::collections::HashMap<String, HmmscanSummary>,
) -> f64 {
    if ref_ids.is_empty() { return 0.0; }
    use std::collections::HashMap;
    let mut freq: HashMap<String, usize> = HashMap::new();
    let mut denom = 0usize;
    for rid in ref_ids {
        if let Some(s) = ref_map.get(rid) {
            let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
            for h in &s.hits {
                let acc = h.accession.as_str();
                if seen.insert(acc) {
                    *freq.entry(acc.to_string()).or_insert(0) += 1;
                }
            }
            denom += 1;
        }
    }
    if denom == 0 { return 0.0; }
    let denom_f = denom as f64;
    let core_thresh = 0.7;
    let acc_thresh = 0.3;
    let mut core: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut acc: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for (k, v) in &freq {
        let f = (*v as f64) / denom_f;
        if f >= core_thresh { core.insert(k.as_str()); }
        else if f >= acc_thresh { acc.insert(k.as_str()); }
    }
    let qset: std::collections::HashSet<&str> = query.hits.iter().map(|h| h.accession.as_str()).collect();
    let core_count = core.len() as f64;
    let recall_core = if core_count > 0.0 {
        let have = qset.iter().filter(|d| core.contains(**d)).count() as f64;
        have / core_count
    } else { 1.0 };
    let dq_minus_core: Vec<&str> = qset.iter().copied().filter(|d| !core.contains(*d)).collect();
    let denom_acc = dq_minus_core.len() as f64;
    let precision_acc = if denom_acc > 0.0 {
        let good = dq_minus_core.iter().filter(|d| acc.contains(**d)).count() as f64;
        good / denom_acc
    } else { 1.0 };
    let extras = dq_minus_core.iter().filter(|d| !acc.contains(**d)).count() as f64;
    let extras_pen = if denom_acc > 0.0 { extras / denom_acc } else { 0.0 };
    let w_core = 0.6; let w_acc = 0.3; let w_extra = 0.1; let w_ord = 0.0;
    let order_pen = 0.0;
    (w_core*recall_core + w_acc*precision_acc - w_extra*extras_pen - w_ord*order_pen).clamp(0.0, 1.0)
}

/// Load Pfam clans mapping from a TSV with columns: Pfam_Acc\tClan_Acc
pub fn load_pfam_clans(path: &str) -> Result<std::collections::HashMap<String, String>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut map = std::collections::HashMap::new();
    for line in text.lines() {
        if line.trim().is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 { continue; }
        map.insert(cols[0].to_string(), cols[1].to_string());
    }
    Ok(map)
}

/// Collapse a domain summary by clan: replace accessions with clan ids when available, and
/// deduplicate hits per clan keeping the lowest e-value/ highest score.
pub fn collapse_by_clan(
    sum: &HmmscanSummary,
    clan_map: &std::collections::HashMap<String, String>,
) -> HmmscanSummary {
    use std::collections::HashMap;
    let mut best: HashMap<String, HmmscanHit> = HashMap::new();
    for h in &sum.hits {
        let key = clan_map.get(&h.accession).cloned().unwrap_or_else(|| h.accession.clone());
        let entry = best.entry(key).or_insert_with(|| HmmscanHit {
            target_name: h.target_name.clone(),
            accession: h.accession.clone(),
            evalue: h.evalue,
            score: h.score,
            bias: h.bias,
        });
        // keep better (lower evalue, then higher score)
        if h.evalue < entry.evalue || (h.evalue == entry.evalue && h.score > entry.score) {
            *entry = HmmscanHit {
                target_name: h.target_name.clone(),
                accession: h.accession.clone(),
                evalue: h.evalue,
                score: h.score,
                bias: h.bias,
            };
        }
    }
    let mut hits: Vec<HmmscanHit> = best.into_values().collect();
    hits.sort_by(|a,b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    let top = hits.first();
    HmmscanSummary {
        hits_count: hits.len(),
        top_accession: top.map(|h| h.accession.clone()),
        top_evalue: top.map(|h| h.evalue),
        hits,
    }
}
