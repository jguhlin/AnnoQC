use crate::config_types::FileConfig;
use crate::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use super::diamond::DiamondStageResult;

pub(in crate::analyze) struct ConsensusStageResult {
    pub(in crate::analyze) grouped: HashMap<String, Vec<diamond::DiamondHitRow>>,
    pub(in crate::analyze) cons_cfg: consensus::ConsensusConfig,
    pub(in crate::analyze) panel_map: HashMap<String, Vec<String>>,
    pub(in crate::analyze) len_map: HashMap<String, LengthSummary>,
    pub(in crate::analyze) panel_prov_map: HashMap<String, PanelProvenanceCounts>,
    pub(in crate::analyze) panel_prov_rows: Vec<(String, PanelProvenanceCounts)>,
    pub(in crate::analyze) refprot_used_panels: u64,
    pub(in crate::analyze) taxonomy_hits_map: HashMap<String, Vec<diamond::DiamondHitRow>>,
    pub(in crate::analyze) consensus_secs: f64,
    pub(in crate::analyze) backfill_used: u64,
    pub(in crate::analyze) backfill_added_total: u64,
}

pub(in crate::analyze) fn build_consensus_stage(
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
