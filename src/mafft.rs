use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
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

    let mut child = Command::new(mafft_bin)
        .arg("--auto")
        .arg("--thread")
        .arg("1")
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
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
    if n == 0 || cols == 0 {
        return AlignmentMetrics::default();
    }
    let mut conserved = 0usize;
    let mut pid_sum = 0usize;
    let mut valid_cols = 0usize;
    for c in 0..cols {
        let mut base = None;
        let mut all_same = true;
        let mut non_gap = 0usize;
        for r in 0..n {
            let ch = seqs[r][c];
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
        if non_gap >= 2 {
            valid_cols += 1;
        }
        if all_same && non_gap == n {
            conserved += 1;
        }
    }
    // pairwise identity: between first and others averaged
    let mut matches = 0usize;
    let mut compared = 0usize;
    for r in 1..n {
        for c in 0..cols {
            let a = seqs[0][c];
            let b = seqs[r][c];
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
    let query = &seqs[0];
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
    }
}
