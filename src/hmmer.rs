use std::io::BufRead;
use std::process::{Command, Stdio};
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

/// Analyze terminal orphan domains for a summary.
///
/// # Parameters
/// - `summary`: Parsed hmmscan summary for a single query.
///
/// # Returns
/// Orphan analysis with a score where 1.0 indicates no orphan penalty.
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
        stdin.flush().map_err(|e| e.to_string())?;
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
    for (line_idx, line) in std::io::BufReader::new(bytes).lines().enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim_start().starts_with('#') || line.trim().is_empty() {
            continue;
        }
        // domtblout columns: target name, accession, tlen, query name, accession, qlen, ... , i-Evalue, score, bias, ...
        let mut target_name = None;
        let mut accession = None;
        let mut hmm_len = None;
        let mut query_len = None;
        let mut i_eval = None;
        let mut score = None;
        let mut bias = None;
        let mut hmm_from = None;
        let mut hmm_to = None;
        let mut ali_from = None;
        let mut ali_to = None;
        let mut env_from = None;
        let mut env_to = None;
        let mut col_count = 0usize;
        for (idx, col) in line.split_whitespace().enumerate() {
            col_count = idx + 1;
            match idx {
                0 => target_name = Some(col),
                1 => accession = Some(col),
                2 => hmm_len = Some(col),
                5 => query_len = Some(col),
                12 => i_eval = Some(col),
                13 => score = Some(col),
                14 => bias = Some(col),
                15 => hmm_from = Some(col),
                16 => hmm_to = Some(col),
                17 => ali_from = Some(col),
                18 => ali_to = Some(col),
                19 => env_from = Some(col),
                20 => env_to = Some(col),
                _ => {}
            }
        }
        if col_count < 22 {
            continue;
        }
        let line_no = line_idx + 1;
        let target_name = target_name
            .ok_or_else(|| format!("domtblout line {} missing target name", line_no))?
            .to_string();
        let accession = accession
            .ok_or_else(|| format!("domtblout line {} missing accession", line_no))?
            .to_string();
        let hmm_len = parse_domtblout_usize(hmm_len, "hmm_len", line_no)?;
        let query_len = parse_domtblout_usize(query_len, "query_len", line_no)?;
        let i_eval = parse_domtblout_f64(i_eval, "i_eval", line_no)?;
        let score = parse_domtblout_f64(score, "score", line_no)?;
        let bias = parse_domtblout_f64(bias, "bias", line_no)?;
        let hmm_from = parse_domtblout_usize(hmm_from, "hmm_from", line_no)?;
        let hmm_to = parse_domtblout_usize(hmm_to, "hmm_to", line_no)?;
        let ali_from = parse_domtblout_usize(ali_from, "ali_from", line_no)?;
        let ali_to = parse_domtblout_usize(ali_to, "ali_to", line_no)?;
        let env_from = parse_domtblout_usize(env_from, "env_from", line_no)?;
        let env_to = parse_domtblout_usize(env_to, "env_to", line_no)?;
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

fn parse_domtblout_usize(field: Option<&str>, name: &str, line_no: usize) -> Result<usize, String> {
    let value = field.ok_or_else(|| format!("domtblout line {} missing {}", line_no, name))?;
    value
        .parse::<usize>()
        .map_err(|e| format!("domtblout line {} invalid {}: {}", line_no, name, e))
}

fn parse_domtblout_f64(field: Option<&str>, name: &str, line_no: usize) -> Result<f64, String> {
    let value = field.ok_or_else(|| format!("domtblout line {} missing {}", line_no, name))?;
    value
        .parse::<f64>()
        .map_err(|e| format!("domtblout line {} invalid {}: {}", line_no, name, e))
}

#[cfg(test)]
mod tests {
    use super::{
        analyze_orphan_domains, collapse_by_clan, domains_architecture_diagnostics,
        domains_architecture_score, parse_domtblout, HmmscanHit, HmmscanSummary, OrphanStatus,
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

    #[test]
    fn orphan_detection_flags_c_only() {
        let mut summary = HmmscanSummary::default();
        summary.hits = vec![HmmscanHit {
            accession: "PF99999".into(),
            hmm_from: 5,
            hmm_to: 170,
            hmm_len: 220,
            ali_from: 10,
            ali_to: 200,
            ..Default::default()
        }];
        summary.hits_count = 1;
        let analysis = analyze_orphan_domains(&summary);
        assert_eq!(analysis.status, OrphanStatus::CTerminal);
        assert!(analysis.score < 1.0);
    }

    #[test]
    fn parse_domtblout_handles_multiple_hits() {
        let data = b"# comment\nPF00001.1 PF00001.1 120 Q1 - 90 1e-10 80 0.0 1 1 1e-10 1e-10 80 0.0 1 80 5 85 3 82 0.95 desc\nPF00002.2 PF00002.2 200 Q1 - 90 5e-5 50 0.0 1 1 5e-5 5e-5 50 0.0 10 160 15 170 12 175 0.80 other\n";
        let sum = parse_domtblout(data).expect("parse domtblout");
        assert_eq!(sum.hits_count, 2);
        assert_eq!(sum.hits.len(), 2);
        assert_eq!(sum.top_accession.as_deref(), Some("PF00001.1"));
        assert!(sum.top_evalue.unwrap() < 1e-9);
    }

    #[test]
    fn parse_domtblout_skips_truncated_lines() {
        let data = b"# comment\nPF00001.1 PF00001.1 120 Q1 - 90 1e-10 80 0.0 1 1 1e-10 1e-10 80 0.0 1 80 5 85 3 82 0.95 desc\nPF00002.2 PF00002.2 200 Q1 - 90 5e-5\n";
        let sum = parse_domtblout(data).expect("parse domtblout");
        assert_eq!(sum.hits_count, 1);
        assert_eq!(sum.hits.len(), 1);
        assert_eq!(sum.top_accession.as_deref(), Some("PF00001.1"));
    }

    #[test]
    fn collapse_by_clan_coalesces_accessions() {
        let mut summary = HmmscanSummary::default();
        summary.hits = vec![
            HmmscanHit {
                accession: "PF00001.27".into(),
                evalue: 1e-30,
                ..Default::default()
            },
            HmmscanHit {
                accession: "PF00002.1".into(),
                evalue: 1e-20,
                ..Default::default()
            },
            HmmscanHit {
                accession: "PF00001.27".into(),
                evalue: 1e-25,
                ..Default::default()
            },
        ];
        let mut clan_map = std::collections::HashMap::new();
        clan_map.insert("PF00001".into(), "CL0001".into());
        clan_map.insert("PF00002".into(), "CL0002".into());
        let collapsed = collapse_by_clan(&summary, &clan_map);
        assert_eq!(collapsed.hits_count, 2);
        assert!(collapsed
            .hits
            .iter()
            .any(|h| h.accession == "CL0001" && h.evalue == 1e-30));
        assert!(collapsed.hits.iter().any(|h| h.accession == "CL0002"));
    }

    fn make_summary(domains: &[(&str, f64)]) -> HmmscanSummary {
        let hits: Vec<HmmscanHit> = domains
            .iter()
            .map(|(acc, eval)| HmmscanHit {
                accession: acc.to_string(),
                evalue: *eval,
                ..Default::default()
            })
            .collect();
        HmmscanSummary {
            hits_count: hits.len(),
            top_accession: hits.first().map(|h| h.accession.clone()),
            top_evalue: hits.first().map(|h| h.evalue),
            hits,
        }
    }

    #[test]
    fn architecture_diagnostics_scores_overlap() {
        let query = make_summary(&[("CL0001", 1e-30), ("CL0002", 1e-20), ("CL9999", 1e-5)]);
        let mut refs = std::collections::HashMap::new();
        refs.insert(
            "ref1".into(),
            make_summary(&[("CL0001", 1e-10), ("CL0002", 1e-9)]),
        );
        refs.insert("ref2".into(), make_summary(&[("CL0001", 1e-8)]));
        let ref_ids = vec!["ref1".to_string(), "ref2".to_string()];
        let diag = domains_architecture_diagnostics(&query, &ref_ids, &refs, 0.0);
        assert_eq!(diag.panel_size, ref_ids.len());
        assert_eq!(diag.refs_with_domains, 2);
        assert_eq!(diag.core_count, 1);
        assert_eq!(diag.accessory_count, 1);
        assert!(diag.recall_core <= 1.0);
        assert!(diag.precision_acc <= 1.0);
        assert!(diag.score >= 0.5);

        let score = domains_architecture_score(&query, &ref_ids, &refs);
        assert!((score - diag.score).abs() < 1e-9);
    }
}

/// Batch hmmscan with optional i-Evalue filtering before truncation.
#[allow(dead_code)]
pub fn run_hmmscan_batch_opts(
    bin: &str,
    db_path: &str,
    items: Vec<(String, Vec<u8>)>,
    threads: usize,
    top_n: usize,
    max_ievalue: Option<f64>,
) -> Result<std::collections::HashMap<String, HmmscanSummary>, String> {
    let nthreads = threads.max(1);
    let mut buckets = vec![Vec::new(); nthreads];
    for (idx, item) in items.into_iter().enumerate() {
        buckets[idx % nthreads].push(item);
    }
    let mut handles = Vec::with_capacity(nthreads);
    for bucket in buckets {
        let bin_s = bin.to_string();
        let db_s = db_path.to_string();
        let max_ev = max_ievalue;
        let handle = thread::spawn(move || -> Result<_, String> {
            let mut local = std::collections::HashMap::new();
            for (gid, seq) in bucket {
                let mut sum = run_hmmscan(&bin_s, &db_s, &gid, &seq)
                    .map_err(|e| format!("hmmscan failed for {}: {}", gid, e))?;
                if let Some(th) = max_ev {
                    sum.hits.retain(|h| h.evalue <= th);
                }
                if sum.hits.len() > top_n {
                    sum.hits.truncate(top_n);
                }
                let trimmed = refresh_summary(sum);
                local.insert(gid, trimmed);
            }
            Ok(local)
        });
        handles.push(handle);
    }
    let mut results: std::collections::HashMap<String, HmmscanSummary> = Default::default();
    for h in handles {
        let local = h
            .join()
            .map_err(|_| "hmmscan thread panicked".to_string())??;
        results.extend(local);
    }
    Ok(results)
}

#[allow(dead_code)]
/// Score domain architecture overlap against a reference panel.
///
/// # Parameters
/// - `query`: Query hmmscan summary (clan-collapsed if desired).
/// - `ref_ids`: Reference ids to consider.
/// - `ref_map`: Reference summaries keyed by id (clan-collapsed if desired).
///
/// # Returns
/// Score in [0.0, 1.0]; returns 0.0 when no valid references are provided.
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
    let (core_thresh, acc_thresh) = normalize_domain_thresholds(0.7, 0.3);
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

/// Calculate domain order penalty by comparing query domain order against reference panel.
/// Uses pairwise inversion counting: compares each pair of adjacent domains in query vs panel.
/// Returns a penalty in [0, 1], where 0 means perfect order match and 1 means completely inverted.
fn calculate_domain_order_penalty(
    query_hits: &[HmmscanHit],
    ref_ids: &[String],
    ref_map: &std::collections::HashMap<String, HmmscanSummary>,
) -> f64 {
    // Extract query domain order (unique accessions in order of occurrence)
    let query_domains: Vec<&str> = query_hits
        .iter()
        .filter(|h| h.evalue < 1e-5) // Only consider significant hits
        .map(|h| h.accession.as_str())
        .collect();

    if query_domains.is_empty() || ref_ids.is_empty() {
        return 0.0;
    }

    // Collect panel domain orders
    let mut panel_orders: Vec<Vec<&str>> = Vec::new();
    for rid in ref_ids {
        if let Some(summary) = ref_map.get(rid) {
            let domains: Vec<&str> = summary
                .hits
                .iter()
                .filter(|h| h.evalue < 1e-5)
                .map(|h| h.accession.as_str())
                .collect();
            if !domains.is_empty() {
                panel_orders.push(domains);
            }
        }
    }

    if panel_orders.is_empty() {
        return 0.0;
    }

    // Calculate average order penalty across all panel sequences
    let mut total_penalty = 0.0;
    let mut comparisons = 0;

    for panel_domains in &panel_orders {
        let penalty = calculate_order_penalty_between(&query_domains, panel_domains);
        total_penalty += penalty;
        comparisons += 1;
    }

    if comparisons == 0 {
        return 0.0;
    }

    (total_penalty / comparisons as f64).min(1.0)
}

/// Calculate order penalty between two domain orderings using pairwise inversion counting.
/// For each pair of domains, check if their relative order matches between query and panel.
fn calculate_order_penalty_between(query: &[&str], panel: &[&str]) -> f64 {
    if query.len() < 2 || panel.len() < 2 {
        return 0.0;
    }

    // Create position maps for quick lookup
    let mut query_pos: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (i, &domain) in query.iter().enumerate() {
        query_pos.insert(domain, i);
    }

    let mut panel_pos: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (i, &domain) in panel.iter().enumerate() {
        panel_pos.insert(domain, i);
    }

    // Count inversions: for each pair of domains in query, check if order matches panel
    let mut inversions = 0;
    let mut total_comparisons = 0;

    // Only compare domains that exist in both query and panel
    let common_domains: Vec<&str> = query_pos
        .keys()
        .cloned()
        .filter(|d| panel_pos.contains_key(*d))
        .collect();

    for i in 0..common_domains.len() {
        for j in (i + 1)..common_domains.len() {
            let d1 = common_domains[i];
            let d2 = common_domains[j];

            let q_order = query_pos.get(d1).unwrap_or(&0) < query_pos.get(d2).unwrap_or(&0);
            let p_order = panel_pos.get(d1).unwrap_or(&0) < panel_pos.get(d2).unwrap_or(&0);

            if q_order != p_order {
                inversions += 1;
            }
            total_comparisons += 1;
        }
    }

    if total_comparisons == 0 {
        return 0.0;
    }

    (inversions as f64 / total_comparisons as f64).min(1.0)
}

/// Same logic as `domains_architecture_score` but returns detailed diagnostics to aid debugging.
/// Assumes inputs are already clan-collapsed if desired by caller.
pub fn domains_architecture_diagnostics(
    query: &HmmscanSummary,
    ref_ids: &[String],
    ref_map: &std::collections::HashMap<String, HmmscanSummary>,
    order_weight: f64,
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
    let (core_thresh, acc_thresh) = normalize_domain_thresholds(0.7, 0.3);
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

    // Calculate domain order penalty
    // Compare query domain order against reference panel domain order
    let order_pen = calculate_domain_order_penalty(&query.hits, ref_ids, ref_map);

    let w_core = 0.6;
    let w_acc = 0.3;
    let w_extra = 0.1;
    let w_ord = order_weight; // Use configurable weight
    diag.score = (w_core * diag.recall_core + w_acc * diag.precision_acc
        - w_extra * diag.extras_pen
        - w_ord * order_pen)
        .clamp(0.0, 1.0);
    diag
}

fn normalize_domain_thresholds(core: f64, accessory: f64) -> (f64, f64) {
    let mut core = core.clamp(0.0, 1.0);
    let mut accessory = accessory.clamp(0.0, 1.0);
    if accessory > core {
        std::mem::swap(&mut core, &mut accessory);
    }
    (core, accessory)
}

fn refresh_summary(mut sum: HmmscanSummary) -> HmmscanSummary {
    sum.hits_count = sum.hits.len();
    sum.top_accession = sum.hits.first().map(|h| h.accession.clone());
    sum.top_evalue = sum.hits.first().map(|h| h.evalue);
    sum
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
