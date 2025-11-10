use std::collections::HashMap;
use std::io::Write;
use std::process::{Command, Stdio};

use needletail::parse_fastx_file;

use crate::taxonomy::canonical_accession;

#[derive(Debug, Clone, Default)]
pub struct AlignmentMetrics {
    pub mafft_enabled: bool,
    pub strategy_used: String,
    pub conserved_fraction: f64,
    pub pairwise_identity: f64,
    pub sequences_aligned: usize,
    pub query_gap_fraction: f64,
    pub gap_run_count: usize,
    pub max_gap_run: usize,
    pub motif_mismatch_fraction: f64,
    pub start_concordance: f64,
    pub start_class: String,
    pub missing_exon_run: usize,
    pub retained_intron_run: usize,
}

pub fn load_sequences_by_ids(
    reference_fasta: &str,
    ids: &[String],
) -> Result<HashMap<String, Vec<u8>>, String> {
    let want: std::collections::HashSet<String> =
        ids.iter().map(|s| canonical_accession(s)).collect();
    let mut found = HashMap::new();
    let mut reader = parse_fastx_file(reference_fasta).map_err(|e| e.to_string())?;
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| e.to_string())?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let key = canonical_accession(&id);
        if want.contains(&key) {
            found.insert(key, rec.seq().to_vec());
        }
    }
    Ok(found)
}

pub fn run_mafft(
    mafft_bin: &str,
    query_id: &str,
    query_seq: &[u8],
    hit_seqs: &HashMap<String, Vec<u8>>,
) -> Result<AlignmentMetrics, String> {
    if hit_seqs.is_empty() {
        return Ok(AlignmentMetrics::default());
    }
    // Write a temporary FASTA with query first and then hits
    let mut fasta_data = Vec::new();
    writeln!(&mut fasta_data, ">{}", query_id).unwrap();
    writeln!(&mut fasta_data, "{}", String::from_utf8_lossy(query_seq)).unwrap();
    for (id, seq) in hit_seqs.iter() {
        writeln!(&mut fasta_data, ">{}", id).unwrap();
        writeln!(&mut fasta_data, "{}", String::from_utf8_lossy(seq)).unwrap();
    }

    let threads_env = std::env::var("MAFFT_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1);
    let mut child = Command::new(mafft_bin)
        .arg("--auto")
        .arg("--quiet")
        .arg("--thread")
        .arg(threads_env.to_string())
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("failed to start mafft: {}", e))?;
    {
        let stdin = child.stdin.as_mut().ok_or("mafft stdin unavailable")?;
        stdin.write_all(&fasta_data).map_err(|e| e.to_string())?;
    }
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!("mafft exited with status {}", output.status));
    }
    let aln = String::from_utf8_lossy(&output.stdout);
    let seqs = parse_fasta_sequences(&aln);
    Ok(compute_alignment_metrics(&seqs))
}

fn parse_fasta_sequences(s: &str) -> Vec<Vec<u8>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    for line in s.lines() {
        if line.starts_with('>') {
            if !cur.is_empty() {
                seqs.push(cur.clone());
                cur.clear();
            }
        } else {
            cur.extend_from_slice(line.as_bytes());
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    seqs
}

fn compute_alignment_metrics(seqs: &[Vec<u8>]) -> AlignmentMetrics {
    if seqs.is_empty() {
        return AlignmentMetrics::default();
    }
    let cols = seqs[0].len();
    let n = seqs.len();
    let query = &seqs[0];
    if n == 0 || cols == 0 {
        return AlignmentMetrics::default();
    }
    let mut conserved = 0usize;
    // formerly tracked pid_sum and valid_cols; not used in current metrics
    for (c, _) in seqs[0].iter().enumerate() {
        let mut base = None;
        let mut all_same = true;
        let mut non_gap = 0usize;
        for seq in seqs.iter() {
            let ch = seq[c];
            if ch != b'-' {
                non_gap += 1;
                if let Some(b) = base {
                    if b != ch {
                        all_same = false;
                    }
                } else {
                    base = Some(ch);
                }
            }
        }
        if all_same && non_gap == n {
            conserved += 1;
        }
    }
    // pairwise identity: between first and others averaged
    let mut matches = 0usize;
    let mut compared = 0usize;
    for seq in seqs.iter().skip(1) {
        for (&a, &b) in query.iter().zip(seq.iter()) {
            if a == b && a != b'-' {
                matches += 1;
            }
            if a != b'-' && b != b'-' {
                compared += 1;
            }
        }
    }
    let pairwise_identity = if compared > 0 {
        matches as f64 / compared as f64
    } else {
        0.0
    };

    // gap metrics on query
    let mut gap_count = 0usize;
    let mut gap_run = 0usize;
    let mut max_gap_run = 0usize;
    let mut gap_runs = 0usize;
    for &ch in query {
        if ch == b'-' {
            gap_count += 1;
            gap_run += 1;
        } else {
            if gap_run > 0 {
                gap_runs += 1;
            }
            max_gap_run = max_gap_run.max(gap_run);
            gap_run = 0;
        }
    }
    if gap_run > 0 {
        gap_runs += 1;
        max_gap_run = max_gap_run.max(gap_run);
    }

    // Start-concordance: compute first non-gap column per sequence
    let starts: Vec<usize> = seqs
        .iter()
        .map(|s| s.iter().position(|&c| c != b'-').unwrap_or(0))
        .collect();
    let query_start = *starts.first().unwrap_or(&0);
    let mut others: Vec<usize> = starts.iter().cloned().skip(1).collect();
    others.sort_unstable();
    let modal = if others.is_empty() {
        query_start
    } else {
        // median as a robust proxy for consensus start
        let mid = others.len() / 2;
        others[mid]
    };
    let diff = if query_start > modal {
        (query_start - modal) as i64
    } else {
        -((modal - query_start) as i64)
    };
    let start_concordance = (1.0 - (diff.unsigned_abs() as f64 / 30.0)).clamp(0.0, 1.0);
    let start_class = if diff.abs() <= 3 {
        "LikelyComplete"
    } else if diff > 3 {
        // query starts later -> likely N-truncated
        "LikelyNTruncated"
    } else {
        // query starts earlier -> likely N-extended
        "LikelyNExtended"
    };

    AlignmentMetrics {
        mafft_enabled: true,
        strategy_used: "auto".to_string(),
        conserved_fraction: if cols > 0 {
            conserved as f64 / cols as f64
        } else {
            0.0
        },
        pairwise_identity,
        sequences_aligned: n,
        query_gap_fraction: if cols > 0 {
            gap_count as f64 / cols as f64
        } else {
            0.0
        },
        gap_run_count: gap_runs,
        max_gap_run,
        motif_mismatch_fraction: 0.0,
        start_concordance,
        start_class: start_class.to_string(),
        missing_exon_run: compute_consensus_gap_run(seqs, true),
        retained_intron_run: compute_consensus_gap_run(seqs, false),
    }
}

// If `query_gap=true`, measure the longest run where query is gap and ≥70% of others are residues.
// If `query_gap=false`, measure the longest run where query is residue and ≥70% of others are gaps.
fn compute_consensus_gap_run(seqs: &[Vec<u8>], query_gap: bool) -> usize {
    if seqs.len() < 2 {
        return 0;
    }
    let n = seqs.len();
    let mut run = 0usize;
    let mut best = 0usize;
    for (c, &query_char) in seqs[0].iter().enumerate() {
        let q = query_char == b'-';
        let others_non_gap = (1..n).filter(|&r| seqs[r][c] != b'-').count();
        let others_gap = (1..n).filter(|&r| seqs[r][c] == b'-').count();
        let cond = if query_gap {
            q && (others_non_gap as f64) / (n as f64 - 1.0) >= 0.7
        } else {
            !q && (others_gap as f64) / (n as f64 - 1.0) >= 0.7
        };
        if cond {
            run += 1;
        } else {
            best = best.max(run);
            run = 0;
        }
    }
    best.max(run)
}
