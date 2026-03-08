use std::collections::HashMap;

use crate::diamond;
use crate::ecs::GeneMetrics;
use crate::structvar;

#[derive(Clone, Debug, Default)]
pub struct StructVarOverrideInputs {
    pub min_hsp_len: Option<usize>,
    pub min_hsp_frac: Option<f64>,
    pub fusion_min_gap: Option<usize>,
    pub dup_max_gap: Option<usize>,
    pub split_delta: Option<f64>,
    pub min_subject_cov: Option<f64>,
    pub orient_majority: Option<f64>,
}

pub fn resolve_thresholds(
    config: &StructVarOverrideInputs,
    cli: &StructVarOverrideInputs,
) -> structvar::StructVarThresholds {
    let mut thresholds = structvar::StructVarThresholds::default();

    if let Some(v) = config.min_hsp_len {
        thresholds.min_strong_len = v;
    }
    if let Some(v) = config.min_hsp_frac {
        thresholds.min_strong_frac = v;
    }
    if let Some(v) = config.fusion_min_gap {
        thresholds.fusion_min_gap = v;
    }
    if let Some(v) = config.dup_max_gap {
        thresholds.dup_max_gap = v;
    }
    if let Some(v) = config.split_delta {
        thresholds.split_delta = v;
    }
    if let Some(v) = config.min_subject_cov {
        thresholds.min_subject_cov = v;
    }
    if let Some(v) = config.orient_majority {
        thresholds.orient_majority = v;
    }

    if let Some(v) = cli.min_hsp_len {
        thresholds.min_strong_len = v;
    }
    if let Some(v) = cli.min_hsp_frac {
        thresholds.min_strong_frac = v;
    }
    if let Some(v) = cli.fusion_min_gap {
        thresholds.fusion_min_gap = v;
    }
    if let Some(v) = cli.dup_max_gap {
        thresholds.dup_max_gap = v;
    }
    if let Some(v) = cli.split_delta {
        thresholds.split_delta = v;
    }
    if let Some(v) = cli.min_subject_cov {
        thresholds.min_subject_cov = v;
    }
    if let Some(v) = cli.orient_majority {
        thresholds.orient_majority = v;
    }

    thresholds
}

pub fn analyze_structvar(
    metrics: &[GeneMetrics],
    grouped: &HashMap<String, Vec<diamond::DiamondHitRow>>,
    thresholds: &structvar::StructVarThresholds,
) -> HashMap<String, structvar::StructVar> {
    let mut out = HashMap::new();
    for m in metrics {
        if let Some(hits) = grouped.get(&m.gene_id) {
            out.insert(m.gene_id.clone(), structvar::analyze(hits, thresholds));
        }
    }
    out
}
