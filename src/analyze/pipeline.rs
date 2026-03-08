use crate::scoring_support::{
    build_component_scores, build_scores_map, export_high_sequences, scoring_thresholds,
};
use crate::*;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::manifest::{
    build_config_snapshot, collect_plugin_manifest, collect_rule_manifest, write_run_manifest,
};

pub(super) struct DiamondStageResult {
    pub(super) diamond_tsv: PathBuf,
    pub(super) diamond_secs: f64,
    pub(super) refprot_grouped: HashMap<String, Vec<diamond::DiamondHitRow>>,
}

pub(super) fn run_diamond_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    log_json: bool,
) -> Result<DiamondStageResult, Box<dyn std::error::Error>> {
    let dia_cfg = DiamondConfig {
        bin: cfg.diamond_bin.clone(),
        db: cfg.db.clone(),
        query_fasta: cfg.fasta.clone(),
        threads: cfg.threads,
        out_dir: cfg.out.clone(),
        out_name: "diamond.blastp.tsv".to_string(),
        retries: 2,
        max_hsps: args
            .diamond_max_hsps
            .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.max_hsps))
            .unwrap_or(5),
    };
    log::info!(
        "diamond: started mode={:?} out_dir={} out_name={}",
        args.diamond_mode,
        cfg.out,
        dia_cfg.out_name
    );
    let t_diamond = step_start("diamond", log_json);
    let diamond_tsv = match args.diamond_mode {
        DiamondMode::Single => blastp_once(&dia_cfg)?,
        DiamondMode::Batch => diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?,
        DiamondMode::Auto => {
            let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
            let threshold = args
                .diamond_auto_threshold
                .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold))
                .unwrap_or(200_000usize);
            if n > threshold {
                diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?
            } else {
                blastp_once(&dia_cfg)?
            }
        }
    };
    let diamond_secs = step_finish("diamond", t_diamond, log_json);
    log::info!(
        "diamond: finished out={} seconds={:.3}",
        diamond_tsv.display(),
        diamond_secs
    );

    let refprot_map_path = Path::new("share/refprot/aves/refprot_proteome_map.tsv");
    let refprot_proteome_map = load_refprot_proteome_map(refprot_map_path);
    let mut refprot_grouped: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
    if let Some(ref_db) = cfg.refprot_db.as_ref() {
        if Path::new(ref_db).exists() {
            let refprot_ok = if let Err(e) = preflight::check_diamond_db(&cfg.diamond_bin, ref_db) {
                log::warn!("refprot db preflight failed; skipping refprot: {}", e);
                false
            } else {
                true
            };
            if refprot_ok {
                let ref_db_size = fs::metadata(ref_db).map(|m| m.len()).unwrap_or(0);
                log::info!(
                    "refprot: db={} size_bytes={} fasta={} mode={:?} out_dir={}",
                    ref_db,
                    ref_db_size,
                    cfg.fasta,
                    args.diamond_mode,
                    cfg.out
                );
                let t_ref = step_start("diamond_refprot", log_json);
                let mut ref_cfg = dia_cfg.clone();
                ref_cfg.out_name = "diamond.refprot.tsv".into();
                ref_cfg.db = ref_db.clone();
                let mut used_chunked = false;
                let ref_tsv_result = match args.diamond_mode {
                    DiamondMode::Batch => {
                        used_chunked = true;
                        diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json)
                    }
                    DiamondMode::Auto => {
                        let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
                        let threshold = args
                            .diamond_auto_threshold
                            .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold))
                            .unwrap_or(200_000usize);
                        if n > threshold {
                            used_chunked = true;
                            diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json)
                        } else {
                            blastp_once(&ref_cfg)
                        }
                    }
                    DiamondMode::Single => blastp_once(&ref_cfg),
                };
                if let Ok(ref_tsv) = ref_tsv_result {
                    let mut ref_bytes = fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                    if ref_bytes == 0 && !used_chunked {
                        log::warn!(
                            "refprot blastp produced empty output; retrying in chunked mode"
                        );
                        let _ = diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json);
                        ref_bytes = fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                    }
                    if ref_bytes == 0 {
                        log::warn!("refprot blastp output is empty; continuing without refprot");
                    } else {
                        refprot_grouped = diamond::parse_tsv_grouped(&ref_tsv, None, Some(100))
                            .unwrap_or_default();
                        annotate_refprot_hits(&mut refprot_grouped, &refprot_proteome_map);
                    }
                } else if let Err(e) = ref_tsv_result {
                    log::warn!("refprot blastp failed; continuing without refprot: {}", e);
                }
                let _ = step_finish("diamond_refprot", t_ref, log_json);
            }
        } else {
            log::warn!("refprot db '{}' not found; skipping", ref_db);
        }
    }

    Ok(DiamondStageResult {
        diamond_tsv,
        diamond_secs,
        refprot_grouped,
    })
}

pub(super) struct ConsensusStageResult {
    pub(super) grouped: HashMap<String, Vec<diamond::DiamondHitRow>>,
    pub(super) cons_cfg: consensus::ConsensusConfig,
    pub(super) panel_map: HashMap<String, Vec<String>>,
    pub(super) len_map: HashMap<String, LengthSummary>,
    pub(super) panel_prov_map: HashMap<String, PanelProvenanceCounts>,
    pub(super) panel_prov_rows: Vec<(String, PanelProvenanceCounts)>,
    pub(super) refprot_used_panels: u64,
    pub(super) taxonomy_hits_map: HashMap<String, Vec<diamond::DiamondHitRow>>,
    pub(super) consensus_secs: f64,
    pub(super) backfill_used: u64,
    pub(super) backfill_added_total: u64,
}

pub(super) fn build_consensus_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    metrics: &[GeneMetrics],
    diamond: &DiamondStageResult,
    log_json: bool,
) -> Result<ConsensusStageResult, Box<dyn std::error::Error>> {
    let t_consensus = step_start("consensus", log_json);
    let qlen_map: HashMap<String, usize> = metrics
        .iter()
        .map(|m| (m.gene_id.clone(), m.length))
        .collect();
    let grouped = diamond::parse_tsv_grouped(&diamond.diamond_tsv, Some(&qlen_map), Some(1000))
        .unwrap_or_default();
    let (cons_cfg, refprot_fallback_cfg) = resolve_consensus_and_refprot_config(
        &ConsensusOverrideInputs {
            min_hits: file_cfg.consensus.as_ref().and_then(|c| c.min_hits),
            max_panel: file_cfg.consensus.as_ref().and_then(|c| c.max_panel),
            filt_qcov: file_cfg.consensus.as_ref().and_then(|c| c.filt_qcov),
            filt_scov: file_cfg.consensus.as_ref().and_then(|c| c.filt_scov),
            filt_evalue: file_cfg.consensus.as_ref().and_then(|c| c.filt_evalue),
            filt_pident: file_cfg.consensus.as_ref().and_then(|c| c.filt_pident),
            redundancy_pident: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.redundancy_pident),
            max_high_identity: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.max_high_identity),
            len_ratio_tolerance: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.len_ratio_tolerance),
            backfill_enabled: file_cfg.consensus.as_ref().and_then(|c| c.backfill_enabled),
            backfill_min_primary_hits: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.backfill_min_primary_hits),
            backfill_max_added: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.backfill_max_added),
            refprot_proteome_cap: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.refprot_proteome_cap),
            diversity_rank_index: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.diversity_rank_index),
            diversity_rank_cap: file_cfg
                .consensus
                .as_ref()
                .and_then(|c| c.diversity_rank_cap),
        },
        &RefProtOverrideInputs {
            trigger_k: file_cfg.refprot.as_ref().and_then(|r| r.trigger_k),
            min_qcov: file_cfg.refprot.as_ref().and_then(|r| r.min_qcov),
            min_scov: file_cfg.refprot.as_ref().and_then(|r| r.min_scov),
            max_evalue: file_cfg.refprot.as_ref().and_then(|r| r.max_evalue),
            min_pident: file_cfg.refprot.as_ref().and_then(|r| r.min_pident),
            max_hits: file_cfg.refprot.as_ref().and_then(|r| r.max_hits),
            proteome_cap: file_cfg.refprot.as_ref().and_then(|r| r.proteome_cap),
        },
        &RefProtOverrideInputs {
            trigger_k: args.refprot_trigger_k,
            min_qcov: args.refprot_min_qcov,
            min_scov: args.refprot_min_scov,
            max_evalue: args.refprot_max_evalue,
            min_pident: args.refprot_min_pident,
            max_hits: args.refprot_max_hits,
            proteome_cap: args.refprot_proteome_cap,
        },
        Path::new("clusters.recluster"),
    );
    let mut refprot_grouped = diamond.refprot_grouped.clone();
    if !refprot_grouped.is_empty() {
        filter_refprot_hits(&mut refprot_grouped, &refprot_fallback_cfg);
    }
    let consensus_output = build_consensus_outputs(
        metrics,
        &grouped,
        &refprot_grouped,
        &cons_cfg,
        &refprot_fallback_cfg,
    );
    let consensus_secs = step_finish("consensus", t_consensus, log_json);

    let panel_stats_rows = consensus_output.panel_stats_rows;
    let panel_len_stats = consensus_output.panel_len_stats;
    let panel_agg_rows = consensus_output.panel_agg_rows;
    let mut backfill_used = 0u64;
    let mut backfill_added_total = 0u64;
    if !panel_stats_rows.is_empty() {
        let panel_dbg = Path::new(&cfg.out).join("panel_debug.csv");
        let mut f = File::create(panel_dbg)?;
        writeln!(f, "gene_id,total_hits,phase,selected,filtered_hits,len_ratio_window_min,len_ratio_window_max,median_len_ratio,len_ratio_min,len_ratio_max,median_pident,high_identity_dropped,backfill_from_clusters,diversity_rank_index,diversity_cap,diversity_keys_used,diversity_skipped")?;
        for (gid, st) in &panel_stats_rows {
            writeln!(
                f,
                "{},{},{},{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{},{},{},{}",
                gid,
                st.total_hits,
                st.phase,
                st.selected,
                st.filtered_hits,
                st.len_ratio_window_min,
                st.len_ratio_window_max,
                st.median_len_ratio,
                st.len_ratio_min,
                st.len_ratio_max,
                st.median_pident,
                st.high_identity_dropped,
                st.backfill_from_clusters,
                st.diversity_rank_index,
                st.diversity_cap,
                st.diversity_keys_used,
                st.diversity_skipped,
            )?;
        }
        if !panel_agg_rows.is_empty() {
            let agg_path = Path::new(&cfg.out).join("panel_agg_debug.csv");
            let mut af = File::create(agg_path)?;
            writeln!(af, "gene_id,subject_id,qcov,scov,len_ratio,bitscore,hit_count,selected,source,quality,diversity_key,len_median,len_mad,expected_len,qlen")?;
            for row in &panel_agg_rows {
                let stats = panel_len_stats
                    .get(&row.gene_id)
                    .cloned()
                    .unwrap_or_default();
                writeln!(
                    af,
                    "{},{},{:.3},{:.3},{:.3},{:.1},{},{},{:?},{:.3},{},{:.3},{},{},{}",
                    row.gene_id,
                    row.subject_id,
                    row.qcov,
                    row.scov,
                    row.len_ratio,
                    row.bitscore,
                    row.hit_count,
                    row.selected,
                    row.source,
                    row.quality,
                    row.diversity_key.clone().unwrap_or_default(),
                    stats.median,
                    stats.mad,
                    stats.expected_len,
                    stats.qlen,
                )?;
            }
        }
        for (_gid, st) in &panel_stats_rows {
            if st.backfill_from_clusters > 0 {
                backfill_used += 1;
                backfill_added_total += st.backfill_from_clusters as u64;
            }
        }
    }

    Ok(ConsensusStageResult {
        grouped,
        cons_cfg,
        panel_map: consensus_output.panel_map,
        len_map: consensus_output.len_map,
        panel_prov_map: consensus_output.panel_prov_map,
        panel_prov_rows: consensus_output.panel_prov_rows,
        refprot_used_panels: consensus_output.refprot_used_panels,
        taxonomy_hits_map: consensus_output.taxonomy_hits_map,
        consensus_secs,
        backfill_used,
        backfill_added_total,
    })
}

pub(super) struct HeavyStageResult {
    pub(super) mafft_requested: bool,
    pub(super) hmmer_requested: bool,
    pub(super) orphan_analysis_enabled: bool,
    pub(super) alignment_map: HashMap<String, AlignmentMetrics>,
    pub(super) alignment_secs: HashMap<String, f64>,
    pub(super) hmmsum_map: HashMap<String, HmmscanSummary>,
    pub(super) hmmer_secs: HashMap<String, f64>,
    pub(super) domains_arch_map: HashMap<String, f64>,
    pub(super) domains_arch_dbg: Vec<DomainsArchDebugRow>,
    pub(super) orphan_map: HashMap<String, hmmer::OrphanAnalysis>,
}

pub(super) fn run_heavy_stage(
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

pub(super) struct EmitStageResult {
    pub(super) checksums: Checksums,
    pub(super) stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    pub(super) taxonomy_setup: taxonomy_support::TaxonomySetup,
    pub(super) taxonomy_enabled_effective: bool,
}

pub(super) fn prepare_emit_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    metrics: &[GeneMetrics],
    diamond: &DiamondStageResult,
    tools: &preflight::ToolVersions,
    calibration: CalibrationSettings,
    report_format: ReportFormat,
    rhai_paths: &[String],
    taxonomy_hits_map: Arc<HashMap<String, Vec<diamond::DiamondHitRow>>>,
) -> Result<EmitStageResult, Box<dyn std::error::Error>> {
    let qlen_map: HashMap<String, usize> = metrics
        .iter()
        .map(|m| (m.gene_id.clone(), m.length))
        .collect();
    let stats =
        Arc::new(parse_tsv_stats(&diamond.diamond_tsv, Some(&qlen_map)).unwrap_or_default());
    let checksums = collect_checksums(cfg)?;
    let snapshot = build_config_snapshot(cfg, args, file_cfg, &calibration, report_format);
    let plugin_manifest = collect_plugin_manifest(&args.plugin);
    let rule_manifest = collect_rule_manifest(rhai_paths);
    write_run_manifest(
        cfg,
        tools,
        &checksums,
        &snapshot,
        &plugin_manifest,
        &rule_manifest,
    )?;

    let mut taxonomy_enabled_effective = args.enable_taxonomy
        || file_cfg
            .taxonomy
            .as_ref()
            .and_then(|t| t.enabled)
            .unwrap_or(false);
    let cache_path = args
        .taxonomy_cache
        .as_deref()
        .or_else(|| {
            file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.cache_path.as_deref())
        })
        .or(file_cfg.taxonomy_cache.as_deref());
    let taxdump_dir = args
        .taxonomy_taxdump_dir
        .as_deref()
        .or_else(|| {
            file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.taxdump_dir.as_deref())
        })
        .or(file_cfg.taxonomy_taxdump_dir.as_deref());
    let tax_cfg = TaxonomyConsensusConfig {
        min_hits: args
            .taxonomy_min_consensus
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_consensus))
            .unwrap_or(5)
            .max(1),
        top_hits: args
            .taxonomy_top_hits
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.top_hits))
            .unwrap_or(20)
            .max(1),
        min_support: args
            .taxonomy_min_support
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.min_support))
            .unwrap_or(0.75),
        coarse_rank_index: args
            .taxonomy_coarse_rank_index
            .or_else(|| file_cfg.taxonomy.as_ref().and_then(|t| t.coarse_rank_index))
            .unwrap_or(1),
        coarse_min_support: args
            .taxonomy_coarse_min_support
            .or_else(|| {
                file_cfg
                    .taxonomy
                    .as_ref()
                    .and_then(|t| t.coarse_min_support)
            })
            .unwrap_or(0.6),
    };
    let taxonomy_setup = build_taxonomy_setup(
        metrics,
        &taxonomy_hits_map,
        stats.as_ref(),
        taxonomy_enabled_effective,
        cache_path,
        cfg.reference_fasta.as_deref(),
        taxdump_dir,
        tax_cfg,
        TaxonomyWarnings {
            expected_domain: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.expected_domain.clone()),
            warn_non_target_min_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_frac)
                .unwrap_or(0.05),
            warn_non_target_min_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_min_hits)
                .unwrap_or(200),
            warn_non_target_strong_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_frac)
                .unwrap_or(0.10),
            warn_non_target_strong_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_non_target_strong_hits)
                .unwrap_or(500),
            warn_genus_min_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_frac)
                .unwrap_or(0.15),
            warn_genus_min_hits: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.warn_genus_min_hits)
                .unwrap_or(300),
            low_coverage_frac: file_cfg
                .taxonomy
                .as_ref()
                .and_then(|t| t.low_coverage_frac)
                .unwrap_or(0.05),
        },
    )?;
    taxonomy_enabled_effective = taxonomy_setup.enabled;

    Ok(EmitStageResult {
        checksums,
        stats,
        taxonomy_setup,
        taxonomy_enabled_effective,
    })
}

pub(super) struct ScoringRenderInputs<'a> {
    pub(super) cfg: &'a EffectiveConfig,
    pub(super) file_cfg: &'a FileConfig,
    pub(super) args: &'a AnalyzeArgs,
    pub(super) metrics: &'a [GeneMetrics],
    pub(super) report_format: ReportFormat,
    pub(super) tools: &'a preflight::ToolVersions,
    pub(super) checksums: Checksums,
    pub(super) stats: Arc<HashMap<String, diamond::DiamondHitStats>>,
    pub(super) taxonomy_setup: taxonomy_support::TaxonomySetup,
    pub(super) intrinsic_map: Arc<HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>>,
    pub(super) alignment_map: Arc<HashMap<String, AlignmentMetrics>>,
    pub(super) hmmsum_map: Arc<HashMap<String, HmmscanSummary>>,
    pub(super) domains_arch_map: Arc<HashMap<String, f64>>,
    pub(super) len_map: Arc<HashMap<String, LengthSummary>>,
    pub(super) orphan_map: Arc<HashMap<String, hmmer::OrphanAnalysis>>,
    pub(super) structvar_map: Arc<HashMap<String, structvar::StructVar>>,
    pub(super) panel_prov_map: HashMap<String, PanelProvenanceCounts>,
    pub(super) domains_arch_dbg: &'a [DomainsArchDebugRow],
    pub(super) rhai_paths: &'a [String],
    pub(super) taxonomy_enabled_effective: bool,
    pub(super) hmmer_requested: bool,
    pub(super) mafft_requested: bool,
    pub(super) orphan_analysis_enabled: bool,
    pub(super) calibration: CalibrationSettings,
    pub(super) log_json: bool,
}

pub(super) struct ScoringRenderResult {
    pub(super) render_summary: ecs::RenderSummary,
    pub(super) emit_secs: f64,
    pub(super) taxonomy_enabled_effective: bool,
    pub(super) plugin_count: usize,
    pub(super) rule_count: usize,
    pub(super) taxsum_map: Arc<HashMap<String, Option<TaxonomyEvidence>>>,
    pub(super) resolver: Option<Arc<taxonomy::TaxonomyResolver>>,
}

pub(super) fn run_scoring_and_render_stage(
    input: ScoringRenderInputs<'_>,
) -> Result<ScoringRenderResult, Box<dyn std::error::Error>> {
    let t_emit = step_start("emit_outputs", input.log_json);
    let resolver = input.taxonomy_setup.resolver.clone();
    let taxsum_map = Arc::clone(&input.taxonomy_setup.taxsum_map);
    let genomic_map = load_genomic_context(
        input.args.gff.as_deref(),
        input.args.genome.as_deref(),
        input.log_json,
    );
    let (rnaseq_enabled, rnaseq_map) = load_rnaseq_data(
        input
            .args
            .rnaseq_file
            .clone()
            .or_else(|| input.file_cfg.rnaseq.as_ref().and_then(|r| r.file.clone())),
    );
    let t_scoring = step_start("scoring", input.log_json);
    let (scores_map_inner, raw_scores_inner) = build_scores_map(
        input.metrics,
        input.stats.as_ref(),
        input.intrinsic_map.as_ref(),
        &input.file_cfg.scoring,
        input.taxonomy_enabled_effective,
        input.hmmer_requested,
        if input.hmmer_requested {
            Some(input.hmmsum_map.as_ref())
        } else {
            None
        },
        Some(input.domains_arch_map.as_ref()),
        Some(input.len_map.as_ref()),
        if input.orphan_analysis_enabled {
            Some(input.orphan_map.as_ref())
        } else {
            None
        },
        if input.taxonomy_enabled_effective {
            Some(taxsum_map.as_ref())
        } else {
            None
        },
        Some(input.alignment_map.as_ref()),
        Some(input.structvar_map.as_ref()),
        genomic_map.as_ref().map(|gm| gm.as_ref()),
        Some(rnaseq_map.as_ref()),
        rnaseq_enabled,
        input.calibration,
        input.args.classify_no_data,
    );
    let comp_map = build_component_scores(
        input.metrics,
        input.stats.as_ref(),
        input.intrinsic_map.as_ref(),
        input.taxonomy_enabled_effective,
        if input.orphan_analysis_enabled {
            Some(input.orphan_map.as_ref())
        } else {
            None
        },
        if input.args.enable_taxonomy {
            Some(taxsum_map.as_ref())
        } else {
            None
        },
        Some(input.alignment_map.as_ref()),
        Some(input.len_map.as_ref()),
        genomic_map.as_ref().map(|gm| gm.as_ref()),
        Some(rnaseq_map.as_ref()),
        rnaseq_enabled,
        &input.file_cfg.scoring,
    );
    let _ = step_finish("scoring", t_scoring, input.log_json);
    if input.args.export_high || input.file_cfg.export_high.unwrap_or(false) {
        let export_high_path = input
            .args
            .export_high_path
            .clone()
            .or_else(|| input.file_cfg.export_high_path.clone());
        if let Some((path, count)) = export_high_sequences(
            &input.cfg.out,
            export_high_path.as_deref(),
            input.metrics,
            input.intrinsic_map.as_ref(),
            &scores_map_inner,
        )? {
            log::info!(
                "exported {} high-scoring genes to {}",
                count,
                path.display()
            );
        } else {
            log::info!("no high-scoring genes to export");
        }
    }
    let scores_map = Arc::new(scores_map_inner);
    let raw_scores_map = Arc::new(raw_scores_inner);
    let comp_map = Arc::new(comp_map);
    let panel_prov_map = Arc::new(input.panel_prov_map);
    let extensions = load_extensions(&input.args.plugin, input.rhai_paths);
    let plugin_count = extensions.plugin_count;
    let rule_count = extensions.rule_count;
    let (th_high, th_med) = scoring_thresholds(&input.file_cfg.scoring);
    let render_ctx = build_render_context(RenderContextInputs {
        stats: Arc::clone(&input.stats),
        intrinsic_map: Arc::clone(&input.intrinsic_map),
        alignment_map: Arc::clone(&input.alignment_map),
        hmmsum_map: Arc::clone(&input.hmmsum_map),
        taxsum_map: Arc::clone(&taxsum_map),
        taxonomy_resolver: resolver.clone(),
        genomic_map: genomic_map.clone(),
        scores_map: Arc::clone(&scores_map),
        raw_scores_map: Arc::clone(&raw_scores_map),
        comp_map: Arc::clone(&comp_map),
        arch_map: Arc::clone(&input.domains_arch_map),
        len_map: Arc::clone(&input.len_map),
        orphan_map: Arc::clone(&input.orphan_map),
        structvar_map: Arc::clone(&input.structvar_map),
        panel_prov_map,
        cov_delta_thresh: input.args.coverage_delta_threshold,
        taxonomy_enabled: input.taxonomy_enabled_effective,
        taxonomy_warnings: input.taxonomy_setup.warnings.clone(),
        orphan_analysis_enabled: input.orphan_analysis_enabled,
        csv_verbose: input.args.csv_verbose,
        mafft_missing_exon_thresh: input
            .args
            .alignment_missing_exon
            .or(input.file_cfg.alignment_missing_exon)
            .unwrap_or(30),
        mafft_retained_intron_thresh: input
            .args
            .alignment_retained_intron
            .or(input.file_cfg.alignment_retained_intron)
            .unwrap_or(30),
        features_string: build_features_string(&FeatureSummary {
            profile_name: &format!("{:?}", input.args.profile),
            has_genomic_input: input.args.gff.is_some() || input.args.genome.is_some(),
            nucleotide_input: input.args.nucleotide,
            taxonomy_enabled: input.taxonomy_enabled_effective,
            genomic_enabled: genomic_map.is_some(),
            plugins_enabled: plugin_count > 0,
            rules_enabled: rule_count > 0,
            hmmer_enabled: input.hmmer_requested,
            orphan_enabled: input.orphan_analysis_enabled,
            alignment_enabled: input.mafft_requested,
        }),
        extensions,
        rnaseq_map: Arc::clone(&rnaseq_map),
        rnaseq_enabled,
        th_high,
        th_med,
    });
    let render_summary = run_render_pipeline(
        &input.cfg.out,
        input.metrics,
        render_ctx,
        resolve_render_max_jobs(
            input.args.render_max_jobs,
            input.file_cfg.render_max_jobs,
            input.cfg.threads,
        ),
        input.tools,
        &input.checksums,
        input.report_format.output_config(input.args.resume),
    )?;
    write_domains_arch_debug(&input.cfg.out, input.domains_arch_dbg)?;
    write_structvar_summary(&input.cfg.out, input.structvar_map.as_ref())?;
    let emit_secs = step_finish("emit_outputs", t_emit, input.log_json);

    Ok(ScoringRenderResult {
        render_summary,
        emit_secs,
        taxonomy_enabled_effective: input.taxonomy_enabled_effective,
        plugin_count,
        rule_count,
        taxsum_map,
        resolver,
    })
}

pub(super) struct TimingInputs {
    pub(super) diamond_secs: f64,
    pub(super) ecs_secs: f64,
    pub(super) intrinsic_secs: f64,
    pub(super) consensus_secs: f64,
    pub(super) emit_secs: f64,
}

pub(super) struct FinalizeInputs<'a> {
    pub(super) cfg: &'a EffectiveConfig,
    pub(super) file_cfg: &'a FileConfig,
    pub(super) args: &'a AnalyzeArgs,
    pub(super) report_format: ReportFormat,
    pub(super) calibration: CalibrationSettings,
    pub(super) metrics: &'a [GeneMetrics],
    pub(super) render_summary: &'a ecs::RenderSummary,
    pub(super) taxonomy_enabled_effective: bool,
    pub(super) hmmer_requested: bool,
    pub(super) mafft_requested: bool,
    pub(super) plugin_count: usize,
    pub(super) rule_count: usize,
    pub(super) taxsum_map: &'a HashMap<String, Option<TaxonomyEvidence>>,
    pub(super) resolver: Option<&'a taxonomy::TaxonomyResolver>,
    pub(super) hmmsum_map: &'a HashMap<String, HmmscanSummary>,
    pub(super) alignment_map: &'a HashMap<String, AlignmentMetrics>,
    pub(super) alignment_secs: &'a HashMap<String, f64>,
    pub(super) hmmer_secs: &'a HashMap<String, f64>,
    pub(super) backfill_used: u64,
    pub(super) backfill_added_total: u64,
    pub(super) refprot_used_panels: u64,
    pub(super) timings: TimingInputs,
}

pub(super) fn finalize_outputs(
    input: FinalizeInputs<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let hmmer_genes = input.hmmsum_map.len() as u64;
    let hmmer_hits_total: u64 = input.hmmsum_map.values().map(|s| s.hits_count as u64).sum();
    let mafft_alignments = input.alignment_map.len() as u64;
    let mut totals = build_base_totals(
        input.metrics.len(),
        hmmer_genes,
        hmmer_hits_total,
        mafft_alignments,
        &format!("{:?}", input.args.diamond_mode),
        input.backfill_used,
        input.backfill_added_total,
        input.refprot_used_panels,
    );
    let steps = build_steps(
        input.timings.diamond_secs,
        input.timings.ecs_secs,
        input.timings.intrinsic_secs,
        input.timings.consensus_secs,
        input.timings.emit_secs,
    );
    if input.taxonomy_enabled_effective {
        augment_totals_from_taxonomy(
            &input.cfg.out,
            &mut totals,
            input.metrics,
            input.taxsum_map,
            input.resolver,
            input
                .file_cfg
                .refprot
                .as_ref()
                .map(|cfg_rp| RefprotSelectionConfig {
                    enabled: cfg_rp.enabled.unwrap_or(false),
                    readme_path: cfg_rp.readme_path.as_deref(),
                    max_scopes: cfg_rp.max_scopes.unwrap_or(2),
                    max_proteomes: cfg_rp.max_proteomes.unwrap_or(10),
                }),
        );
    }
    write_slowest_genes(&input.cfg.out, input.alignment_secs, input.hmmer_secs)?;
    let run_summary = build_run_summary_obj(RunSummaryInputs {
        schema_version: OUTPUT_SCHEMA_VERSION,
        total_genes: input.metrics.len(),
        render_summary: input.render_summary,
        aligner: format!("{:?}", input.args.aligner).to_lowercase(),
        report_format: format!("{:?}", input.report_format).to_lowercase(),
        resume: input.args.resume,
        taxonomy_enabled: input.taxonomy_enabled_effective,
        hmmer_enabled: input.hmmer_requested,
        alignment_enabled: input.mafft_requested,
        genomic_enabled: input.args.gff.is_some(),
        plugins: input.plugin_count,
        rules: input.rule_count,
        nucleotide: input.args.nucleotide,
        calibration: format!("{:?}", input.calibration.mode).to_lowercase(),
        timings: steps.clone(),
    });
    write_run_summary(&input.cfg.out, &run_summary)?;
    let runm = build_run_metrics_obj(OUTPUT_SCHEMA_VERSION, steps, totals);
    write_run_metrics(&input.cfg.out, &runm)?;
    Ok(())
}
