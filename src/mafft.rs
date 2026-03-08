use std::collections::HashMap;
use std::io::Write;
use std::process::{Command, Stdio};

use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use spoa::{AlignmentEngine, AlignmentType, Graph};

use crate::taxonomy::canonical_accession;

/// Represents a single conserved block in the alignment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservedBlock {
    pub start: usize,
    pub end: usize,
    pub length: usize,
    pub conservation_fraction: f64,
    pub panel_count: usize,
    pub query_present: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignerBackend {
    Mafft,
    Spoa,
}

impl std::fmt::Display for AlignerBackend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignerBackend::Mafft => write!(f, "mafft"),
            AlignerBackend::Spoa => write!(f, "spoa"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct AlignerConfig {
    pub backend: AlignerBackend,
    /// MAFFT executable path or command name (resolved via PATH).
    pub mafft_bin: String,
    /// Use MAFFT fast mode (FFT-NS-1) instead of auto-tuned settings.
    pub mafft_fast: bool,
    /// MAFFT threads to use per alignment job (minimum 1).
    pub mafft_threads_per_job: usize,
    /// Maximum number of MAFFT jobs to run concurrently.
    pub mafft_max_jobs: usize,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AlignmentResult {
    pub gene_id: String,
    pub seq_ids: Vec<String>,
    pub aligned_seqs: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct AlignmentMetrics {
    pub mafft_enabled: bool,
    pub strategy_used: String,
    pub conserved_fraction: f64,
    pub pairwise_identity: f64,
    pub panel_pairwise_identity: f64,
    pub divergence_ratio: f64,
    pub sequences_aligned: usize,
    pub query_gap_fraction: f64,
    pub gap_run_count: usize,
    pub max_gap_run: usize,
    pub motif_mismatch_fraction: f64,
    pub start_concordance: f64,
    pub start_class: String,
    pub end_concordance: f64,
    pub end_class: String,
    pub missing_exon_run: usize,
    pub retained_intron_run: usize,
    // Block-level conservation metrics
    pub conserved_blocks: Vec<ConservedBlock>,
    pub missing_blocks: Vec<ConservedBlock>,
    pub extra_blocks: Vec<ConservedBlock>,
    pub total_conserved_query: usize,
    pub total_conserved_panel: usize,
    pub block_conservation_score: f64,
}

/// Result of block-level conservation analysis
#[derive(Debug, Clone, Default)]
struct BlockConservationResult {
    conserved_blocks: Vec<ConservedBlock>,
    missing_blocks: Vec<ConservedBlock>,
    extra_blocks: Vec<ConservedBlock>,
    total_conserved_query: usize,
    total_conserved_panel: usize,
    score: f64,
}

/// Analyze block-level conservation in alignment
/// Identifies conserved blocks and compares query vs panel
fn analyze_block_conservation(sequences: &[Vec<u8>], query_idx: usize) -> BlockConservationResult {
    if sequences.is_empty() {
        return BlockConservationResult::default();
    }

    let seq_len = sequences[0].len();
    let window_size = 5; // Minimum block size
    let min_conservation = 0.8; // 80% panel agreement
    let min_panel_support = 0.5; // 50% of sequences

    let mut conserved_blocks: Vec<ConservedBlock> = Vec::new();
    let mut in_block = false;
    let mut block_start = 0;
    let mut block_conservation_sum = 0.0;

    // Sliding window analysis
    for pos in 0..seq_len {
        // Calculate conservation at this position
        let (panel_cons, count) = calculate_position_conservation(sequences, query_idx, pos);

        if panel_cons >= min_conservation
            && count as f64 >= (sequences.len() as f64 * min_panel_support)
        {
            // Start or extend block
            if !in_block {
                block_start = pos;
                in_block = true;
            }
            block_conservation_sum += panel_cons;
        } else if in_block {
            // End block
            let block_len = pos - block_start;
            if block_len >= window_size {
                let avg_conservation = block_conservation_sum / block_len as f64;
                let query_has_block = sequences[query_idx][block_start..pos]
                    .iter()
                    .any(|&c| c != b'-');

                conserved_blocks.push(ConservedBlock {
                    start: block_start,
                    end: pos,
                    length: block_len,
                    conservation_fraction: avg_conservation,
                    panel_count: count,
                    query_present: query_has_block,
                });
            }
            in_block = false;
            block_conservation_sum = 0.0;
        }
    }

    // Handle block that extends to end
    if in_block {
        let block_len = seq_len - block_start;
        if block_len >= window_size {
            let avg_conservation = block_conservation_sum / block_len as f64;
            let query_has_block = sequences[query_idx][block_start..seq_len]
                .iter()
                .any(|&c| c != b'-');

            conserved_blocks.push(ConservedBlock {
                start: block_start,
                end: seq_len,
                length: block_len,
                conservation_fraction: avg_conservation,
                panel_count: sequences.len() - 1,
                query_present: query_has_block,
            });
        }
    }

    // Classify blocks
    let missing_blocks: Vec<ConservedBlock> = conserved_blocks
        .iter()
        .filter(|b| !b.query_present)
        .cloned()
        .collect();

    let extra_blocks: Vec<ConservedBlock> = conserved_blocks
        .iter()
        .filter(|b| b.query_present && b.panel_count < (sequences.len() / 3))
        .cloned()
        .collect();

    let total_conserved_query: usize = conserved_blocks
        .iter()
        .filter(|b| b.query_present)
        .map(|b| b.length)
        .sum();

    let total_conserved_panel: usize = conserved_blocks.iter().map(|b| b.length).sum();

    // Calculate score
    let score = if total_conserved_panel == 0 {
        1.0 // No blocks = no penalty
    } else {
        let query_coverage = total_conserved_query as f64 / total_conserved_panel as f64;
        let missing_penalty: f64 = missing_blocks
            .iter()
            .map(|b| {
                let weight = b.conservation_fraction; // High conservation = high penalty
                b.length as f64 * weight
            })
            .sum();
        let total_block_length = total_conserved_panel as f64;
        if total_block_length > 0.0 {
            let normalized_missing = missing_penalty / total_block_length;
            (query_coverage - normalized_missing).clamp(0.0, 1.0)
        } else {
            query_coverage.clamp(0.0, 1.0)
        }
    };

    BlockConservationResult {
        conserved_blocks,
        missing_blocks,
        extra_blocks,
        total_conserved_query,
        total_conserved_panel,
        score,
    }
}

/// Calculate conservation fraction at a specific alignment position
fn calculate_position_conservation(
    sequences: &[Vec<u8>],
    query_idx: usize,
    pos: usize,
) -> (f64, usize) {
    let panel_seqs: Vec<_> = sequences
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != query_idx)
        .map(|(_, seq)| seq.get(pos).copied().unwrap_or(b'-'))
        .collect();

    if panel_seqs.is_empty() {
        return (0.0, 0);
    }

    // Find consensus residue (most common non-gap character)
    let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
    for &c in &panel_seqs {
        if c != b'-' {
            *counts.entry(c).or_insert(0) += 1;
        }
    }

    let consensus = counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(c, _)| c)
        .unwrap_or(b'-');

    let count = panel_seqs
        .iter()
        .filter(|&&c| c == consensus || c == b'-')
        .count();

    let fraction = if panel_seqs.is_empty() {
        0.0
    } else {
        count as f64 / panel_seqs.len() as f64
    };

    (fraction, count)
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

/// Run MAFFT on the provided sequences (query expected first) and return the aligned rows.
pub fn run_mafft_alignment(
    cfg: &AlignerConfig,
    gene_id: &str,
    ids_and_seqs: &[(String, String)],
) -> Result<AlignmentResult, String> {
    if ids_and_seqs.is_empty() {
        return Err("no sequences provided for alignment".into());
    }

    #[cfg(test)]
    if cfg.mafft_bin.is_empty() {
        // In tests we allow bypassing the external MAFFT binary to keep
        // alignment unit tests self-contained.
        return Ok(AlignmentResult {
            gene_id: gene_id.to_string(),
            seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
            aligned_seqs: ids_and_seqs.iter().map(|(_, s)| s.clone()).collect(),
        });
    }

    if cfg.mafft_bin.trim().is_empty() {
        return Err("mafft_bin is empty; set a MAFFT executable path".into());
    }
    let bin_path = std::path::Path::new(&cfg.mafft_bin);
    if bin_path.parent().is_some() && !bin_path.exists() {
        return Err(format!("mafft_bin not found: {}", cfg.mafft_bin));
    }

    let threads = cfg.mafft_threads_per_job.max(1);
    let mut cmd = Command::new(&cfg.mafft_bin);
    if cfg.mafft_fast {
        // FFT-NS-1: fastest reasonable settings
        cmd.arg("--retree").arg("1").arg("--maxiterate").arg("0");
    } else {
        cmd.arg("--auto");
    }
    let mut child = cmd
        .arg("--quiet")
        .arg("--thread")
        .arg(threads.to_string())
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("failed to start mafft: {}", e))?;
    {
        let stdin = child.stdin.as_mut().ok_or("mafft stdin unavailable")?;
        for (id, seq) in ids_and_seqs.iter() {
            writeln!(stdin, ">{}", id).map_err(|e| e.to_string())?;
            writeln!(stdin, "{}", seq).map_err(|e| e.to_string())?;
        }
    }
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!("mafft exited with status {}", output.status));
    }
    let aln = String::from_utf8_lossy(&output.stdout);
    let seqs = parse_fasta_sequences(&aln)
        .map_err(|e| format!("mafft output parse failed: {}", e))?
        .into_iter()
        .map(|s| String::from_utf8_lossy(&s).into_owned())
        .collect::<Vec<_>>();

    if seqs.len() != ids_and_seqs.len() {
        return Err(format!(
            "mafft MSA row count mismatch: expected {} got {}",
            ids_and_seqs.len(),
            seqs.len()
        ));
    }

    Ok(AlignmentResult {
        gene_id: gene_id.to_string(),
        seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
        aligned_seqs: seqs,
    })
}

fn parse_fasta_sequences(s: &str) -> Result<Vec<Vec<u8>>, String> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    let mut saw_header = false;
    for line in s.lines() {
        if line.starts_with('>') {
            saw_header = true;
            if !cur.is_empty() {
                seqs.push(std::mem::take(&mut cur));
            }
        } else {
            if !saw_header {
                return Err("FASTA missing header before sequence data".into());
            }
            cur.extend_from_slice(line.as_bytes());
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    if seqs.is_empty() {
        return Err("no sequences parsed from FASTA output".into());
    }
    Ok(seqs)
}

pub const MIN_PANEL_FOR_TERMINI: usize = 5;
pub const MIN_PANEL_FOR_CONSERVED_REGIONS: usize = 10;

pub fn run_spoa_alignment(
    gene_id: &str,
    ids_and_seqs: &[(String, String)],
) -> Result<AlignmentResult, String> {
    if ids_and_seqs.is_empty() {
        return Err("no sequences provided for alignment".to_string());
    }
    if ids_and_seqs.len() < 2 {
        return Ok(AlignmentResult {
            gene_id: gene_id.to_string(),
            seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
            aligned_seqs: ids_and_seqs.iter().map(|(_, s)| s.clone()).collect(),
        });
    }

    // Simple protein scoring; can be tuned later.
    let mut engine = AlignmentEngine::new(AlignmentType::kNW, 2, -1, -2, -1, -2, -1);
    let mut graph = Graph::new();

    for (_, seq) in ids_and_seqs.iter() {
        let s = seq.as_bytes();
        let alignment = engine.align(s, &graph);
        graph.add_alignment(&alignment, s, 1);
    }

    let msa = graph.multiple_sequence_alignment(false);

    let aligned_seqs: Vec<String> = msa
        .iter()
        .map(|row| String::from_utf8_lossy(row).into_owned())
        .collect();

    if aligned_seqs.len() != ids_and_seqs.len() {
        return Err(format!(
            "spoa MSA row count mismatch: expected {} got {}",
            ids_and_seqs.len(),
            aligned_seqs.len()
        ));
    }

    validate_spoa_alignment(ids_and_seqs, &aligned_seqs)
        .map_err(|e| format!("parity check failed: {}", e))?;

    Ok(AlignmentResult {
        gene_id: gene_id.to_string(),
        seq_ids: ids_and_seqs.iter().map(|(id, _)| id.clone()).collect(),
        aligned_seqs,
    })
}

fn validate_spoa_alignment(
    ids_and_seqs: &[(String, String)],
    aligned_seqs: &[String],
) -> Result<(), String> {
    if aligned_seqs.is_empty() {
        return Err("alignment is empty".to_string());
    }
    let expected_len = aligned_seqs[0].len();
    if expected_len == 0 {
        return Err("alignment length is zero".to_string());
    }
    for (idx, aln) in aligned_seqs.iter().enumerate() {
        if aln.len() != expected_len {
            return Err(format!(
                "row {} length mismatch: expected {} got {}",
                idx,
                expected_len,
                aln.len()
            ));
        }
    }
    for ((id, seq), aln) in ids_and_seqs.iter().zip(aligned_seqs.iter()) {
        let non_gap = aln.as_bytes().iter().filter(|&&c| c != b'-').count();
        if non_gap != seq.len() {
            return Err(format!(
                "sequence {} residue count mismatch: expected {} got {}",
                id,
                seq.len(),
                non_gap
            ));
        }
    }
    Ok(())
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
    // pairwise identity: between first (query) and others averaged
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

    // panel pairwise identity: average identity among the reference sequences (1..n)
    // To keep this fast (O(N^2 * L)), we only compute if n > 2.
    // If n=2 (query + 1 ref), panel identity is undefined/1.0, effectively same as pairwise.
    let mut panel_matches = 0usize;
    let mut panel_compared = 0usize;
    if n > 2 {
        for i in 1..n {
            for j in (i + 1)..n {
                let s1 = &seqs[i];
                let s2 = &seqs[j];
                for (&a, &b) in s1.iter().zip(s2.iter()) {
                    if a == b && a != b'-' {
                        panel_matches += 1;
                    }
                    if a != b'-' && b != b'-' {
                        panel_compared += 1;
                    }
                }
            }
        }
    }

    let panel_pairwise_identity = if n <= 2 {
        pairwise_identity // Fallback for single reference
    } else if panel_compared > 0 {
        panel_matches as f64 / panel_compared as f64
    } else {
        0.0
    };

    // Divergence ratio: How similar is the query compared to how similar the panel is?
    // 1.0 = Query is as central as any panel member.
    // < 0.5 = Query is significantly divergent/outlier.
    let divergence_ratio = if panel_pairwise_identity <= 1.0e-6 {
        1.0
    } else if panel_pairwise_identity > 0.0 {
        pairwise_identity / panel_pairwise_identity
    } else {
        1.0 // If panel is garbage/diverse, we can't judge the query harshly
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
    let mut start_concordance = 0.0;
    let mut start_class = "InsufficientPanel".to_string();
    let mut end_concordance = 0.0;
    let mut end_class = "InsufficientPanel".to_string();
    if n >= MIN_PANEL_FOR_TERMINI {
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
            let mid = others.len() / 2;
            others[mid]
        };
        let diff = query_start as i64 - modal as i64;
        start_concordance = (1.0 - (diff.unsigned_abs() as f64 / 30.0)).clamp(0.0, 1.0);
        start_class = if diff.abs() <= 3 {
            "LikelyComplete"
        } else if diff > 3 {
            "LikelyNTruncated"
        } else {
            "LikelyNExtended"
        }
        .to_string();

        let ends: Vec<usize> = seqs
            .iter()
            .map(|s| s.iter().rposition(|&c| c != b'-').unwrap_or(0))
            .collect();
        let query_end = *ends.first().unwrap_or(&0);
        let mut others_end: Vec<usize> = ends.iter().cloned().skip(1).collect();
        others_end.sort_unstable();
        let modal_end = if others_end.is_empty() {
            query_end
        } else {
            let mid = others_end.len() / 2;
            others_end[mid]
        };
        let diff_end = query_end as i64 - modal_end as i64;
        end_concordance = (1.0 - (diff_end.unsigned_abs() as f64 / 30.0)).clamp(0.0, 1.0);
        end_class = if diff_end.abs() <= 3 {
            "LikelyComplete"
        } else if diff_end < -3 {
            "LikelyCTruncated"
        } else {
            "LikelyCExtended"
        }
        .to_string();
    }

    // Block-level conservation analysis
    let block_result = analyze_block_conservation(seqs, 0);

    AlignmentMetrics {
        mafft_enabled: true,
        strategy_used: "auto".to_string(),
        conserved_fraction: if cols > 0 {
            conserved as f64 / cols as f64
        } else {
            0.0
        },
        pairwise_identity,
        panel_pairwise_identity,
        divergence_ratio,
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
        end_concordance,
        end_class,
        missing_exon_run: compute_consensus_gap_run(seqs, true),
        retained_intron_run: compute_consensus_gap_run(seqs, false),
        conserved_blocks: block_result.conserved_blocks,
        missing_blocks: block_result.missing_blocks,
        extra_blocks: block_result.extra_blocks,
        total_conserved_query: block_result.total_conserved_query,
        total_conserved_panel: block_result.total_conserved_panel,
        block_conservation_score: block_result.score,
    }
}

/// Compute metrics for an alignment result and tag with the strategy used.
fn metrics_from_alignment(result: &AlignmentResult, strategy: &str) -> AlignmentMetrics {
    let seqs: Vec<Vec<u8>> = result
        .aligned_seqs
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let mut metrics = compute_alignment_metrics(&seqs);
    metrics.strategy_used = strategy.to_string();
    metrics
}

/// Run alignment according to the configured backend, falling back to MAFFT on SPOA failure.
pub fn run_alignment_for_panel(
    cfg: &AlignerConfig,
    gene_id: &str,
    ids_and_seqs: &[(String, String)],
) -> Result<AlignmentMetrics, String> {
    match cfg.backend {
        AlignerBackend::Mafft => {
            let aln = run_mafft_alignment(cfg, gene_id, ids_and_seqs)?;
            Ok(metrics_from_alignment(&aln, "mafft"))
        }
        AlignerBackend::Spoa => {
            #[cfg(test)]
            let force_error = std::env::var("FORCE_SPOA_ERROR")
                .map(|v| v == "1")
                .unwrap_or(false);

            #[cfg(not(test))]
            let force_error = false;

            let spoa_result = if force_error {
                Err("forced spoa error".to_string())
            } else {
                run_spoa_alignment(gene_id, ids_and_seqs)
            };
            match spoa_result {
                Ok(aln) => Ok(metrics_from_alignment(&aln, "spoa")),
                Err(e) => {
                    log::warn!("SPOA failed for {}: {}; falling back to MAFFT", gene_id, e);
                    let aln = run_mafft_alignment(cfg, gene_id, ids_and_seqs)?;
                    Ok(metrics_from_alignment(&aln, "mafft_fallback"))
                }
            }
        }
    }
}

/// If `query_gap=true`, measure the longest run where query is gap and >=70% of others are residues.
/// If `query_gap=false`, measure the longest run where query is residue and >=70% of others are gaps.
fn compute_consensus_gap_run(seqs: &[Vec<u8>], query_gap: bool) -> usize {
    if seqs.len() < 2 {
        return 0;
    }
    let n = seqs.len();
    let mut run = 0usize;
    let mut best = 0usize;
    for (c, &query_char) in seqs[0].iter().enumerate() {
        let q = query_char == b'-';
        let mut others_non_gap = 0usize;
        for other in seqs.iter().skip(1) {
            if other[c] != b'-' {
                others_non_gap += 1;
            }
        }
        let others_gap = (n - 1).saturating_sub(others_non_gap);
        let denom = (n - 1) as f64;
        let cond = if query_gap {
            q && (others_non_gap as f64) / denom >= 0.7
        } else {
            !q && (others_gap as f64) / denom >= 0.7
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

#[cfg(test)]
mod tests {
    use super::{
        compute_alignment_metrics, run_alignment_for_panel, run_spoa_alignment,
        validate_spoa_alignment, AlignerBackend, AlignerConfig,
    };
    use proptest::prelude::*;

    #[test]
    fn start_concordance_flags_truncation() {
        let query: [u8; 12] = [
            b'-', b'-', b'-', b'-', b'M', b'K', b'T', b'A', b'A', b'-', b'-', b'-',
        ];
        let ref_row: [u8; 12] = [
            b'M', b'K', b'T', b'A', b'A', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
        ];
        let mut seqs = vec![query.to_vec()];
        for _ in 0..4 {
            seqs.push(ref_row.to_vec());
        }
        let metrics = compute_alignment_metrics(&seqs);
        assert_eq!(metrics.start_class, "LikelyNTruncated");
        assert!(metrics.start_concordance < 1.0);
    }

    #[test]
    fn end_concordance_flags_extension() {
        let query: [u8; 12] = [
            b'M', b'K', b'T', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A',
        ];
        let ref_row: [u8; 12] = [
            b'M', b'K', b'T', b'A', b'A', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
        ];
        let mut seqs = vec![query.to_vec()];
        for _ in 0..4 {
            seqs.push(ref_row.to_vec());
        }
        let metrics = compute_alignment_metrics(&seqs);
        assert_eq!(metrics.end_class, "LikelyCExtended");
        assert!(metrics.end_concordance < 1.0);
    }

    #[test]
    fn spoa_round_trip_rows_and_lengths() {
        let ids_and_seqs = vec![
            ("query".to_string(), "MKTAA".to_string()),
            ("hit1".to_string(), "MKAAA".to_string()),
            ("hit2".to_string(), "MKTAG".to_string()),
        ];
        let result = run_spoa_alignment("gene1", &ids_and_seqs).expect("spoa succeeds");
        assert_eq!(result.seq_ids.len(), ids_and_seqs.len());
        assert_eq!(result.aligned_seqs.len(), ids_and_seqs.len());
        let lens: Vec<usize> = result.aligned_seqs.iter().map(|s| s.len()).collect();
        assert!(lens.windows(2).all(|w| w[0] == w[1]));
    }

    #[test]
    fn spoa_falls_back_to_mafft_with_empty_bin_in_tests() {
        std::env::set_var("FORCE_SPOA_ERROR", "1");
        let cfg = AlignerConfig {
            backend: AlignerBackend::Spoa,
            mafft_bin: String::new(),
            mafft_fast: false,
            mafft_threads_per_job: 1,
            mafft_max_jobs: 1,
        };
        let ids_and_seqs = vec![
            ("query".to_string(), "ACDE".to_string()),
            ("hit1".to_string(), "ACDE".to_string()),
        ];
        let metrics = run_alignment_for_panel(&cfg, "gene1", &ids_and_seqs)
            .expect("fallback alignment succeeds");
        assert!(metrics.mafft_enabled);
        assert_eq!(metrics.strategy_used, "mafft_fallback");
        std::env::remove_var("FORCE_SPOA_ERROR");
    }

    #[test]
    fn spoa_parity_check_flags_residue_mismatch() {
        let ids_and_seqs = vec![
            ("query".to_string(), "ACDE".to_string()),
            ("hit1".to_string(), "ACDE".to_string()),
        ];
        let aligned = vec!["ACDE".to_string(), "ACD-".to_string()];
        let err = validate_spoa_alignment(&ids_and_seqs, &aligned)
            .expect_err("should detect residue mismatch");
        assert!(err.contains("residue count mismatch"));
    }

    proptest! {
        #[test]
        fn alignment_metrics_invariants(
            len in 1usize..60,
            nrefs in 1usize..4,
            seqs in prop::collection::vec(
                prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T', b'-']), 1..60),
                2..6
            ),
        ) {
            let total = nrefs + 1;
            let mut aligned: Vec<Vec<u8>> = Vec::with_capacity(total);
            for i in 0..total {
                let mut s = seqs[i % seqs.len()].clone();
                s.truncate(len);
                if s.len() < len {
                    s.resize(len, b'A');
                }
                aligned.push(s);
            }
            let metrics = compute_alignment_metrics(&aligned);
            assert!(metrics.start_concordance >= 0.0 && metrics.start_concordance <= 1.0);
            assert!(metrics.end_concordance >= 0.0 && metrics.end_concordance <= 1.0);
            assert!(metrics.divergence_ratio.is_finite());
            let query = &aligned[0];
            let gap_count = query.iter().filter(|&&c| c == b'-').count();
            if gap_count == 0 {
                assert_eq!(metrics.gap_run_count, 0);
                assert_eq!(metrics.max_gap_run, 0);
            } else {
                assert!(metrics.gap_run_count <= gap_count);
                assert!(metrics.max_gap_run <= gap_count);
            }
        }
    }
}
