use std::collections::HashMap;

use crate::ecs::{AlignmentPipelineConfig, GeneMetrics};
use crate::mafft::{load_sequences_by_ids, AlignerBackend, AlignerConfig};
use crate::metrics;

#[derive(Clone, Debug, Default)]
pub struct AlignmentOverrideInputs {
    pub threads_per_job: Option<usize>,
    pub max_jobs: Option<usize>,
    pub fast: Option<bool>,
}

pub struct AlignmentSetup {
    pub backend: AlignerBackend,
    pub job_count: usize,
    pub pipeline: Option<AlignmentPipelineConfig>,
}

pub fn build_alignment_setup(
    metrics: &[GeneMetrics],
    panel_map: &HashMap<String, Vec<String>>,
    intrinsic_map: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    ref_fasta: &str,
    mafft_bin: &str,
    backend: AlignerBackend,
    total_threads: usize,
    min_panel_hits: usize,
    overrides: &AlignmentOverrideInputs,
) -> AlignmentSetup {
    let total_threads = total_threads.max(1);
    let default_per_job = if total_threads >= 16 {
        8
    } else if total_threads >= 8 {
        4
    } else if total_threads >= 4 {
        2
    } else {
        1
    };

    let mafft_threads_per_job = overrides
        .threads_per_job
        .unwrap_or(default_per_job)
        .clamp(1, total_threads);
    let default_workers = (total_threads / mafft_threads_per_job).max(1);
    let requested_workers = overrides.max_jobs.unwrap_or(default_workers).max(1);
    let mafft_workers = requested_workers.min(default_workers).max(1);
    let mafft_fast = overrides.fast.unwrap_or(false);

    std::env::set_var("MAFFT_THREADS", mafft_threads_per_job.to_string());
    if mafft_fast {
        std::env::set_var("MAFFT_FAST", "1");
    } else {
        std::env::remove_var("MAFFT_FAST");
    }

    let all_ids: Vec<String> = panel_map.values().flat_map(|v| v.clone()).collect();
    let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
    let query_seq_map: HashMap<String, Vec<u8>> = intrinsic_map
        .iter()
        .map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone()))
        .collect();
    let jobs: Vec<(String, Vec<String>)> = metrics
        .iter()
        .filter_map(|g| {
            let ids = panel_map.get(&g.gene_id)?.clone();
            if ids.len() < min_panel_hits {
                return None;
            }
            if !query_seq_map.contains_key(&g.gene_id) {
                return None;
            }
            Some((g.gene_id.clone(), ids))
        })
        .collect();

    let job_count = jobs.len();
    let pipeline = if jobs.is_empty() {
        None
    } else {
        Some(AlignmentPipelineConfig {
            aligner: AlignerConfig {
                backend,
                mafft_bin: mafft_bin.to_string(),
                mafft_fast,
                mafft_threads_per_job,
                mafft_max_jobs: mafft_workers,
            },
            jobs,
            query_map: query_seq_map,
            ref_map: ref_seqs,
        })
    };

    AlignmentSetup {
        backend,
        job_count,
        pipeline,
    }
}
