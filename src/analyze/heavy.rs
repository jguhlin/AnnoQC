use crate::config_types::{FileConfig, ScoringConfigOverride};
use crate::*;
use std::collections::HashMap;

pub(in crate::analyze) struct HeavyStageResult {
    pub(in crate::analyze) mafft_requested: bool,
    pub(in crate::analyze) hmmer_requested: bool,
    pub(in crate::analyze) orphan_analysis_enabled: bool,
    pub(in crate::analyze) alignment_map: HashMap<String, AlignmentMetrics>,
    pub(in crate::analyze) alignment_secs: HashMap<String, f64>,
    pub(in crate::analyze) hmmsum_map: HashMap<String, HmmscanSummary>,
    pub(in crate::analyze) hmmer_secs: HashMap<String, f64>,
    pub(in crate::analyze) domains_arch_map: HashMap<String, f64>,
    pub(in crate::analyze) domains_arch_dbg: Vec<DomainsArchDebugRow>,
    pub(in crate::analyze) orphan_map: HashMap<String, hmmer::OrphanAnalysis>,
}

pub(in crate::analyze) fn run_heavy_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    metrics: &[GeneMetrics],
    panel_map: &HashMap<String, Vec<String>>,
    intrinsic_map: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    min_hits: usize,
    log_json: bool,
) -> HeavyStageResult {
    let mut alignment_map = HashMap::new();
    let mut alignment_secs = HashMap::new();
    let mut align_cfg_opt = None;
    let mut t_align_opt = None;
    if let (Some(ref_fasta), Some(mafft_bin)) =
        (cfg.reference_fasta.as_ref(), args.mafft_bin.as_ref())
    {
        let backend = if let Some(b) = file_cfg.mafft_backend.as_deref() {
            match b.to_ascii_lowercase().as_str() {
                "spoa" => AlignerBackend::Spoa,
                _ => AlignerBackend::Mafft,
            }
        } else {
            args.aligner.into()
        };
        t_align_opt = Some(step_start("alignment", log_json));
        let alignment_setup = build_alignment_setup(
            metrics,
            panel_map,
            intrinsic_map,
            ref_fasta,
            mafft_bin,
            backend,
            cfg.threads,
            min_hits,
            &AlignmentOverrideInputs {
                threads_per_job: args
                    .mafft_threads_per_job
                    .or(file_cfg.mafft_threads_per_job),
                max_jobs: args.mafft_max_jobs.or(file_cfg.mafft_max_jobs),
                fast: Some(if args.mafft_fast {
                    true
                } else {
                    file_cfg.mafft_fast.unwrap_or(false)
                }),
            },
        );
        log::info!("alignment backend: {}", alignment_setup.backend);
        log::info!(
            "alignment jobs queued (backend: {}, jobs: {})",
            alignment_setup.backend,
            alignment_setup.job_count
        );
        align_cfg_opt = alignment_setup.pipeline;
    }

    let mut hmmsum_map = HashMap::new();
    let mut hmmer_secs = HashMap::new();
    let mut hmmer_cfg_opt;
    let mut hmmer_timer = None;
    let hmmer_items: Vec<(String, Vec<u8>)> = intrinsic_map
        .iter()
        .map(|(gid, (_im, qseq))| (gid.clone(), qseq.clone()))
        .collect();
    let hmmer_setup = build_hmmer_setup(
        args.hmmscan_bin.as_deref(),
        args.pfam_db
            .as_deref()
            .or(file_cfg.pfam_db.as_deref())
            .or(file_cfg.pfam_db.as_deref()),
        hmmer_items,
        cfg.threads,
        &HmmerOverrideInputs {
            top_n: args.hmmer_top_n,
            threads: args.hmmer_threads,
            query_ievalue: args.hmmer_ievalue,
            ref_ievalue: args.hmmer_ref_ievalue,
            orphan_analysis: None,
        },
        &HmmerOverrideInputs {
            top_n: file_cfg.hmmer.as_ref().and_then(|h| h.top_n),
            threads: file_cfg.hmmer.as_ref().and_then(|h| h.threads),
            query_ievalue: file_cfg.hmmer.as_ref().and_then(|h| h.ievalue),
            ref_ievalue: file_cfg.hmmer.as_ref().and_then(|h| h.ref_ievalue),
            orphan_analysis: file_cfg.hmmer.as_ref().and_then(|h| h.orphan_analysis),
        },
        args.disable_orphan_analysis,
    );
    let hmmer_top_n = hmmer_setup.top_n;
    let hmmer_thread_cap = hmmer_setup.thread_cap;
    let hmmer_ref_ievalue = hmmer_setup.ref_ievalue;
    let orphan_analysis_enabled = hmmer_setup.orphan_analysis_enabled;
    let hmmer_bin_owned = hmmer_setup.bin_owned;
    let hmmer_db_owned = hmmer_setup.db_owned;
    hmmer_cfg_opt = hmmer_setup.config;
    log::info!("hmmscan jobs queued: {}", hmmer_setup.job_count);
    if hmmer_cfg_opt.is_some() {
        hmmer_timer = Some(step_start("hmmer", log_json));
    }

    let mafft_requested = align_cfg_opt.is_some();
    let hmmer_requested = hmmer_cfg_opt.is_some();
    if mafft_requested || hmmer_requested {
        let reserve_threads = cfg.threads.min(2);
        let heavy_results = run_heavy_pipelines(HeavyPipelineConfig {
            cpu_threads: cfg.threads,
            reserve_threads,
            log_json,
            alignment: align_cfg_opt.take(),
            hmmer: hmmer_cfg_opt.take(),
        });
        alignment_map = heavy_results.alignment_map;
        hmmsum_map = heavy_results.hmmer_map;
        alignment_secs = heavy_results.alignment_secs;
        hmmer_secs = heavy_results.hmmer_secs;
    }
    if let Some(t) = t_align_opt {
        let _ = step_finish("alignment", t, log_json);
    }
    if let Some(t) = hmmer_timer {
        let _ = step_finish("hmmer", t, log_json);
    }

    let domain_order_weight: f64 = file_cfg
        .scoring
        .as_ref()
        .and_then(|s: &ScoringConfigOverride| s.weights.get("domains_order"))
        .copied()
        .unwrap_or(0.0);
    let hmmer_post = postprocess_hmmer(HmmerPostProcessInputs {
        hmm_bin: hmmer_bin_owned.as_deref(),
        db_path: hmmer_db_owned.as_deref(),
        reference_fasta: cfg.reference_fasta.as_deref(),
        panel_map,
        hmmsum_map: &hmmsum_map,
        metrics,
        cfg_threads: cfg.threads,
        log_json,
        hmmer_thread_cap,
        hmmer_top_n,
        hmmer_ref_ievalue,
        orphan_analysis_enabled,
        pfam_clans_path: args
            .pfam_clans
            .as_deref()
            .or(file_cfg.pfam_clans.as_deref()),
        domain_order_weight,
    });

    HeavyStageResult {
        mafft_requested,
        hmmer_requested,
        orphan_analysis_enabled,
        alignment_map,
        alignment_secs,
        hmmsum_map,
        hmmer_secs,
        domains_arch_map: hmmer_post.domains_arch_map,
        domains_arch_dbg: hmmer_post.domains_arch_dbg,
        orphan_map: hmmer_post.orphan_map,
    }
}
