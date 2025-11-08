use std::io::BufRead;
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;

#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct HmmscanHit {
    pub target_name: String,
    pub accession: String,
    pub evalue: f64,
    pub score: f64,
    pub bias: f64,
    pub hmm_from: usize,
    pub hmm_to: usize,
    pub hmm_len: usize,
    pub ali_from: usize,
    pub ali_to: usize,
    pub env_from: usize,
    pub env_to: usize,
    pub query_len: usize,
}

#[derive(Debug, Clone, Default)]
pub struct HmmscanSummary {
    pub hits_count: usize,
    pub top_accession: Option<String>,
    pub top_evalue: Option<f64>,
    pub hits: Vec<HmmscanHit>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OrphanStatus {
    #[default]
    None,
    NTerminal,
    CTerminal,
    Both,
}

impl OrphanStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            OrphanStatus::None => "None",
            OrphanStatus::NTerminal => "NTerminalOrphan",
            OrphanStatus::CTerminal => "CTerminalOrphan",
            OrphanStatus::Both => "BothOrphans",
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct OrphanDomainDetail {
    pub accession: String,
    pub domain_index: usize,
    pub total_domains: usize,
    pub completeness: f64,
    pub hmm_from: usize,
    pub hmm_to: usize,
    pub hmm_len: usize,
}

#[derive(Debug, Clone)]
pub struct OrphanAnalysis {
    pub status: OrphanStatus,
    pub score: f64,
    pub details: Vec<OrphanDomainDetail>,
}

impl Default for OrphanAnalysis {
    fn default() -> Self {
        Self {
            status: OrphanStatus::None,
            score: 1.0,
            details: Vec::new(),
        }
    }
}

const ORPHAN_MARGIN: usize = 20;
const ORPHAN_MIN_MODEL: usize = 60;

pub fn analyze_orphan_domains(summary: &HmmscanSummary) -> OrphanAnalysis {
    if summary.hits.is_empty() {
        return OrphanAnalysis::default();
    }
    let mut hits = summary.hits.clone();
    hits.sort_by(|a, b| a.ali_from.cmp(&b.ali_from));
    let total = hits.len();
    let mut flagged = Vec::new();
    let mut has_n = false;
    let mut has_c = false;

    let first = &hits[0];
    if first.hmm_len >= ORPHAN_MIN_MODEL && first.hmm_from > ORPHAN_MARGIN {
        has_n = true;
        flagged.push(build_orphan_detail(first, 0, total));
    }
    let last = &hits[total - 1];
    let missing_c = last.hmm_len.saturating_sub(last.hmm_to);
    if last.hmm_len >= ORPHAN_MIN_MODEL && missing_c > ORPHAN_MARGIN {
        has_c = true;
        if total == 1 {
            // avoid duplicating the same detail twice for single-domain proteins
            if !has_n {
                flagged.push(build_orphan_detail(last, total - 1, total));
            }
        } else {
            flagged.push(build_orphan_detail(last, total - 1, total));
        }
    }

    let status = match (has_n, has_c) {
        (true, true) => OrphanStatus::Both,
        (true, false) => OrphanStatus::NTerminal,
        (false, true) => OrphanStatus::CTerminal,
        (false, false) => OrphanStatus::None,
    };
    if status == OrphanStatus::None {
        return OrphanAnalysis::default();
    }
    let score = (1.0 - 0.5 * (has_n as u8 + has_c as u8) as f64).clamp(0.0, 1.0);
    OrphanAnalysis {
        status,
        score,
        details: flagged,
    }
}

fn build_orphan_detail(hit: &HmmscanHit, index: usize, total: usize) -> OrphanDomainDetail {
    let covered = if hit.hmm_to >= hit.hmm_from {
        hit.hmm_to - hit.hmm_from + 1
    } else {
        0
    };
    let completeness = if hit.hmm_len > 0 {
        covered as f64 / hit.hmm_len as f64
    } else {
        0.0
    };
    OrphanDomainDetail {
        accession: hit.accession.clone(),
        domain_index: index + 1,
        total_domains: total,
        completeness,
        hmm_from: hit.hmm_from,
        hmm_to: hit.hmm_to,
        hmm_len: hit.hmm_len,
    }
}

pub fn run_hmmscan(
    bin: &str,
    db_path: &str,
    id: &str,
    seq: &[u8],
) -> Result<HmmscanSummary, String> {
    let tmpdir = std::env::temp_dir();
    let domtbl = tmpdir.join(format!("hmmscan_{}_{}.domtblout", std::process::id(), id));
    let mut child = Command::new(bin)
        .arg("-o")
        .arg("/dev/null")
        .arg("--noali")
        .arg("--domtblout")
        .arg(&domtbl)
        .arg(db_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
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
    let status = child.wait().map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!("hmmscan exited with status {}", status));
    }
    let bytes = std::fs::read(&domtbl).map_err(|e| format!("read domtblout: {}", e))?;
    let res = parse_domtblout(&bytes[..]);
    let _ = std::fs::remove_file(&domtbl);
    res
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
        if cols.len() < 22 {
            continue;
        }
        let target_name = cols[0].to_string();
        let accession = cols[1].to_string();
        let hmm_len = cols[2].parse::<usize>().unwrap_or(0);
        let query_len = cols[5].parse::<usize>().unwrap_or(0);
        let i_eval = cols[12].parse::<f64>().unwrap_or(1.0);
        let score = cols
            .get(13)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let bias = cols
            .get(14)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        let hmm_from = cols
            .get(15)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let hmm_to = cols
            .get(16)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let ali_from = cols
            .get(17)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let ali_to = cols
            .get(18)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let env_from = cols
            .get(19)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let env_to = cols
            .get(20)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        hits.push(HmmscanHit {
            target_name,
            accession,
            evalue: i_eval,
            score,
            bias,
            hmm_from,
            hmm_to,
            hmm_len,
            ali_from,
            ali_to,
            env_from,
            env_to,
            query_len,
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
    use super::{
        analyze_orphan_domains, parse_domtblout, HmmscanHit, HmmscanSummary, OrphanStatus,
    };

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
        assert_eq!(h.hmm_len, 250);
        assert_eq!(h.query_len, 120);
        assert_eq!(h.hmm_from, 1);
        assert_eq!(h.hmm_to, 100);
        assert_eq!(h.ali_from, 1);
        assert_eq!(h.ali_to, 100);
    }

    #[test]
    fn orphan_detection_flags_n_and_c() {
        let mut summary = HmmscanSummary::default();
        summary.hits = vec![
            HmmscanHit {
                accession: "PF00001".into(),
                hmm_from: 35,
                hmm_to: 120,
                hmm_len: 200,
                ali_from: 5,
                ali_to: 140,
                ..Default::default()
            },
            HmmscanHit {
                accession: "PF00002".into(),
                hmm_from: 10,
                hmm_to: 215,
                hmm_len: 220,
                ali_from: 200,
                ali_to: 360,
                ..Default::default()
            },
        ];
        summary.hits_count = summary.hits.len();
        let analysis = analyze_orphan_domains(&summary);
        assert_eq!(analysis.status, OrphanStatus::NTerminal);
        assert!(analysis.score < 1.0);
    }

    #[test]
    fn orphan_detection_handles_both_ends() {
        let mut summary = HmmscanSummary::default();
        summary.hits = vec![HmmscanHit {
            accession: "PF00003".into(),
            hmm_from: 40,
            hmm_to: 160,
            hmm_len: 300,
            ali_from: 5,
            ali_to: 270,
            ..Default::default()
        }];
        summary.hits_count = 1;
        let analysis = analyze_orphan_domains(&summary);
        assert_eq!(analysis.status, OrphanStatus::Both);
        assert_eq!(analysis.details.len(), 1);
        assert!(analysis.score <= 0.0 + f64::EPSILON);
    }
}

/// Batch hmmscan with optional i-Evalue filtering before truncation.
pub fn run_hmmscan_batch_opts(
    bin: &str,
    db_path: &str,
    items: Vec<(String, Vec<u8>)>,
    threads: usize,
    top_n: usize,
    max_ievalue: Option<f64>,
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
        let max_ev = max_ievalue;
        let handle = thread::spawn(move || loop {
            let next = {
                let mut guard = q.lock().unwrap();
                guard.next()
            };
            let Some((gid, seq)) = next else {
                break;
            };
            let mut sum = run_hmmscan(&bin_s, &db_s, &gid, &seq).unwrap_or_default();
            if let Some(th) = max_ev {
                sum.hits.retain(|h| h.evalue <= th);
            }
            let mut trimmed = sum.clone();
            if trimmed.hits.len() > top_n {
                trimmed.hits.truncate(top_n);
            }
            let mut out = r.lock().unwrap();
            out.insert(gid, trimmed);
        });
        handles.push(handle);
    }
    for h in handles {
        h.join()
            .map_err(|_| "hmmscan thread panicked".to_string())?;
    }
    let map = Arc::try_unwrap(results)
        .map_err(|_| "results arc busy".to_string())
        .and_then(|m| m.into_inner().map_err(|_| "results poisoned".to_string()))?;
    Ok(map)
}

#[allow(dead_code)]
pub fn domains_architecture_score(
    query: &HmmscanSummary,
    ref_ids: &[String],
    ref_map: &std::collections::HashMap<String, HmmscanSummary>,
) -> f64 {
    if ref_ids.is_empty() {
        return 0.0;
    }
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
    if denom == 0 {
        return 0.0;
    }
    let denom_f = denom as f64;
    let core_thresh = 0.7;
    let acc_thresh = 0.3;
    let mut core: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut acc: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for (k, v) in &freq {
        let f = (*v as f64) / denom_f;
        if f >= core_thresh {
            core.insert(k.as_str());
        } else if f >= acc_thresh {
            acc.insert(k.as_str());
        }
    }
    let qset: std::collections::HashSet<&str> =
        query.hits.iter().map(|h| h.accession.as_str()).collect();
    let core_count = core.len() as f64;
    let recall_core = if core_count > 0.0 {
        let have = qset.iter().filter(|d| core.contains(**d)).count() as f64;
        have / core_count
    } else {
        1.0
    };
    let dq_minus_core: Vec<&str> = qset
        .iter()
        .copied()
        .filter(|d| !core.contains(*d))
        .collect();
    let denom_acc = dq_minus_core.len() as f64;
    let precision_acc = if denom_acc > 0.0 {
        let good = dq_minus_core.iter().filter(|d| acc.contains(**d)).count() as f64;
        good / denom_acc
    } else {
        1.0
    };
    let extras = dq_minus_core.iter().filter(|d| !acc.contains(**d)).count() as f64;
    let extras_pen = if denom_acc > 0.0 {
        extras / denom_acc
    } else {
        0.0
    };
    let w_core = 0.6;
    let w_acc = 0.3;
    let w_extra = 0.1;
    let w_ord = 0.0;
    let order_pen = 0.0;
    (w_core * recall_core + w_acc * precision_acc - w_extra * extras_pen - w_ord * order_pen)
        .clamp(0.0, 1.0)
}

#[derive(Debug, Clone, Default)]
pub struct DomainsArchDiagnostics {
    pub panel_size: usize,
    pub refs_with_domains: usize,
    pub query_domains: usize,
    pub core_count: usize,
    pub accessory_count: usize,
    pub overlap_core: usize,
    pub overlap_accessory: usize,
    pub extras_count: usize,
    pub recall_core: f64,
    pub precision_acc: f64,
    pub extras_pen: f64,
    pub score: f64,
}

/// Same logic as `domains_architecture_score` but returns detailed diagnostics to aid debugging.
/// Assumes inputs are already clan-collapsed if desired by caller.
pub fn domains_architecture_diagnostics(
    query: &HmmscanSummary,
    ref_ids: &[String],
    ref_map: &std::collections::HashMap<String, HmmscanSummary>,
) -> DomainsArchDiagnostics {
    use std::collections::{HashMap, HashSet};
    // Always record query-side domain count (useful even if refs have none)
    let qset: HashSet<&str> = query.hits.iter().map(|h| h.accession.as_str()).collect();
    let mut diag = DomainsArchDiagnostics {
        panel_size: ref_ids.len(),
        query_domains: qset.len(),
        ..Default::default()
    };
    if ref_ids.is_empty() {
        return diag;
    }
    let mut freq: HashMap<String, usize> = HashMap::new();
    let mut denom = 0usize;
    for rid in ref_ids {
        if let Some(s) = ref_map.get(rid) {
            let mut seen: HashSet<&str> = HashSet::new();
            let mut had = false;
            for h in &s.hits {
                let acc = h.accession.as_str();
                if seen.insert(acc) {
                    *freq.entry(acc.to_string()).or_insert(0) += 1;
                    had = true;
                }
            }
            if had {
                diag.refs_with_domains += 1;
            }
            denom += 1;
        }
    }
    if denom == 0 {
        return diag;
    }
    let denom_f = denom as f64;
    let core_thresh = 0.7;
    let acc_thresh = 0.3;
    let mut core: HashSet<&str> = HashSet::new();
    let mut acc: HashSet<&str> = HashSet::new();
    for (k, v) in &freq {
        let f = (*v as f64) / denom_f;
        if f >= core_thresh {
            core.insert(k.as_str());
        } else if f >= acc_thresh {
            acc.insert(k.as_str());
        }
    }
    diag.core_count = core.len();
    diag.accessory_count = acc.len();
    // qset/diag.query_domains already computed above
    let core_count_f = diag.core_count as f64;
    diag.overlap_core = qset.iter().filter(|d| core.contains(**d)).count();
    diag.recall_core = if core_count_f > 0.0 {
        diag.overlap_core as f64 / core_count_f
    } else {
        1.0
    };
    let dq_minus_core: Vec<&str> = qset
        .iter()
        .copied()
        .filter(|d| !core.contains(*d))
        .collect();
    let denom_acc = dq_minus_core.len() as f64;
    diag.overlap_accessory = dq_minus_core.iter().filter(|d| acc.contains(**d)).count();
    diag.precision_acc = if denom_acc > 0.0 {
        diag.overlap_accessory as f64 / denom_acc
    } else {
        1.0
    };
    diag.extras_count = dq_minus_core.iter().filter(|d| !acc.contains(**d)).count();
    diag.extras_pen = if denom_acc > 0.0 {
        diag.extras_count as f64 / denom_acc
    } else {
        0.0
    };
    let w_core = 0.6;
    let w_acc = 0.3;
    let w_extra = 0.1;
    let w_ord = 0.0;
    let order_pen = 0.0;
    diag.score = (w_core * diag.recall_core + w_acc * diag.precision_acc
        - w_extra * diag.extras_pen
        - w_ord * order_pen)
        .clamp(0.0, 1.0);
    diag
}

/// Load Pfam clans mapping from a TSV with columns: Pfam_Acc\tClan_Acc
pub fn load_pfam_clans(path: &str) -> Result<std::collections::HashMap<String, String>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut map = std::collections::HashMap::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
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
        // Normalize PFAM accession by dropping version suffix (PF00476.27 -> PF00476)
        let base_acc = h
            .accession
            .split('.')
            .next()
            .unwrap_or(&h.accession)
            .to_string();
        // Use clan id as the grouping key; if no clan, fall back to base accession
        let key = clan_map.get(&base_acc).cloned().unwrap_or(base_acc);
        // Store the collapsed hit using the clan id (or base accession) in accession field
        let entry = best.entry(key.clone()).or_insert_with(|| {
            let mut clone = h.clone();
            clone.accession = key.clone();
            clone
        });
        // keep better (lower evalue, then higher score)
        if h.evalue < entry.evalue || (h.evalue == entry.evalue && h.score > entry.score) {
            let mut clone = h.clone();
            clone.accession = key.clone();
            *entry = clone;
        }
    }
    let mut hits: Vec<HmmscanHit> = best.into_values().collect();
    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let top = hits.first();
    HmmscanSummary {
        hits_count: hits.len(),
        top_accession: top.map(|h| h.accession.clone()),
        top_evalue: top.map(|h| h.evalue),
        hits,
    }
}
