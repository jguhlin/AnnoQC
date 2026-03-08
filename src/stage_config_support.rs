use std::path::Path;
use std::sync::Arc;

use crate::clusters;
use crate::consensus;
use crate::consensus_support::RefProtFallbackConfig;

#[derive(Clone, Debug, Default)]
pub struct ConsensusOverrideInputs {
    pub min_hits: Option<usize>,
    pub max_panel: Option<usize>,
    pub filt_qcov: Option<f64>,
    pub filt_scov: Option<f64>,
    pub filt_evalue: Option<f64>,
    pub filt_pident: Option<f64>,
    pub redundancy_pident: Option<f64>,
    pub max_high_identity: Option<usize>,
    pub len_ratio_tolerance: Option<f64>,
    pub backfill_enabled: Option<bool>,
    pub backfill_min_primary_hits: Option<usize>,
    pub backfill_max_added: Option<usize>,
    pub refprot_proteome_cap: Option<usize>,
    pub diversity_rank_index: Option<usize>,
    pub diversity_rank_cap: Option<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct RefProtOverrideInputs {
    pub trigger_k: Option<usize>,
    pub min_qcov: Option<f64>,
    pub min_scov: Option<f64>,
    pub max_evalue: Option<f64>,
    pub min_pident: Option<f64>,
    pub max_hits: Option<usize>,
    pub proteome_cap: Option<usize>,
}

pub fn resolve_consensus_and_refprot_config(
    consensus_inputs: &ConsensusOverrideInputs,
    refprot_config_inputs: &RefProtOverrideInputs,
    refprot_cli_inputs: &RefProtOverrideInputs,
    cluster_map_path: &Path,
) -> (consensus::ConsensusConfig, RefProtFallbackConfig) {
    let mut cons_cfg = consensus::ConsensusConfig::default();
    if let Some(v) = consensus_inputs.min_hits {
        cons_cfg.min_hits = v;
    }
    if let Some(v) = consensus_inputs.max_panel {
        cons_cfg.max_panel = v;
    }
    if let Some(v) = consensus_inputs.filt_qcov {
        cons_cfg.filt_qcov = v;
    }
    if let Some(v) = consensus_inputs.filt_scov {
        cons_cfg.filt_scov = v;
    }
    if let Some(v) = consensus_inputs.filt_evalue {
        cons_cfg.filt_evalue = v;
    }
    if let Some(v) = consensus_inputs.filt_pident {
        cons_cfg.filt_pident = v;
    }
    if let Some(v) = consensus_inputs.redundancy_pident {
        cons_cfg.redundancy_pident = v;
    }
    if let Some(v) = consensus_inputs.max_high_identity {
        cons_cfg.max_high_identity = v;
    }
    if let Some(v) = consensus_inputs.len_ratio_tolerance {
        cons_cfg.len_ratio_tolerance = v;
    }
    if let Some(v) = consensus_inputs.backfill_enabled {
        cons_cfg.backfill_enabled = v;
    }
    if let Some(v) = consensus_inputs.backfill_min_primary_hits {
        cons_cfg.backfill_min_primary_hits = v;
    }
    if let Some(v) = consensus_inputs.backfill_max_added {
        cons_cfg.backfill_max_added = v;
    }
    if let Some(v) = consensus_inputs.refprot_proteome_cap {
        cons_cfg.refprot_proteome_cap = v;
    }
    if let Some(v) = consensus_inputs.diversity_rank_index {
        cons_cfg.diversity_rank_index = v;
    }
    if let Some(v) = consensus_inputs.diversity_rank_cap {
        cons_cfg.diversity_rank_cap = v;
    }

    let mut refprot_fallback_cfg = RefProtFallbackConfig {
        trigger_k: cons_cfg.min_hits,
        ..Default::default()
    };
    if let Some(v) = refprot_config_inputs.trigger_k {
        refprot_fallback_cfg.trigger_k = v;
    }
    if let Some(v) = refprot_config_inputs.min_qcov {
        refprot_fallback_cfg.min_qcov = v;
    }
    if let Some(v) = refprot_config_inputs.min_scov {
        refprot_fallback_cfg.min_scov = v;
    }
    if let Some(v) = refprot_config_inputs.max_evalue {
        refprot_fallback_cfg.max_evalue = v;
    }
    if let Some(v) = refprot_config_inputs.min_pident {
        refprot_fallback_cfg.min_pident = v;
    }
    if let Some(v) = refprot_config_inputs.max_hits {
        refprot_fallback_cfg.max_hits = v.max(1);
    }
    if let Some(v) = refprot_config_inputs.proteome_cap {
        cons_cfg.refprot_proteome_cap = v;
    }

    if let Some(v) = refprot_cli_inputs.trigger_k {
        refprot_fallback_cfg.trigger_k = v;
    }
    if let Some(v) = refprot_cli_inputs.min_qcov {
        refprot_fallback_cfg.min_qcov = v;
    }
    if let Some(v) = refprot_cli_inputs.min_scov {
        refprot_fallback_cfg.min_scov = v;
    }
    if let Some(v) = refprot_cli_inputs.max_evalue {
        refprot_fallback_cfg.max_evalue = v;
    }
    if let Some(v) = refprot_cli_inputs.min_pident {
        refprot_fallback_cfg.min_pident = v;
    }
    if let Some(v) = refprot_cli_inputs.max_hits {
        refprot_fallback_cfg.max_hits = v.max(1);
    }
    if let Some(v) = refprot_cli_inputs.proteome_cap {
        cons_cfg.refprot_proteome_cap = v;
    }

    if cluster_map_path.exists() {
        match clusters::load_cluster_map(cluster_map_path.to_str().unwrap()) {
            Ok(map) => {
                cons_cfg.clusters = Some(Arc::new(map));
                log::info!(
                    "consensus: loaded cluster map from {}",
                    cluster_map_path.display()
                );
            }
            Err(e) => log::warn!("consensus: failed to load cluster map: {}", e),
        }
    }

    (cons_cfg, refprot_fallback_cfg)
}
