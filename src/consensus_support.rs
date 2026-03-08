use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use crate::consensus;
use crate::diamond;
use crate::ecs::GeneMetrics;
use crate::length;
use crate::HitSource;

pub type LengthSummary = (f64, f64, f64, String, f64, f64, bool, usize);

#[derive(Clone)]
pub struct PanelAggDebugRow {
    pub gene_id: String,
    pub subject_id: String,
    pub qcov: f64,
    pub scov: f64,
    pub len_ratio: f64,
    pub bitscore: f64,
    pub hit_count: usize,
    pub selected: bool,
    pub source: HitSource,
    pub quality: f64,
    pub diversity_key: Option<String>,
}

#[derive(Clone, Default)]
pub struct PanelProvenanceCounts {
    pub swissprot: usize,
    pub refprot: usize,
    pub cluster: usize,
}

#[derive(Clone)]
pub struct RefProtFallbackConfig {
    pub trigger_k: usize,
    pub min_qcov: f64,
    pub min_scov: f64,
    pub max_evalue: f64,
    pub min_pident: f64,
    pub max_hits: usize,
}

impl Default for RefProtFallbackConfig {
    fn default() -> Self {
        Self {
            trigger_k: 5,
            min_qcov: 0.50,
            min_scov: 0.25,
            max_evalue: 1e-10,
            min_pident: 30.0,
            max_hits: 100,
        }
    }
}

pub struct ConsensusBuildOutput {
    pub panel_map: HashMap<String, Vec<String>>,
    pub panel_stats_rows: Vec<(String, consensus::PanelStats)>,
    pub panel_len_stats: HashMap<String, consensus::LenStats>,
    pub panel_agg_rows: Vec<PanelAggDebugRow>,
    pub len_map: HashMap<String, LengthSummary>,
    pub panel_prov_map: HashMap<String, PanelProvenanceCounts>,
    pub panel_prov_rows: Vec<(String, PanelProvenanceCounts)>,
    pub taxonomy_hits_map: HashMap<String, Vec<diamond::DiamondHitRow>>,
    pub refprot_used_panels: u64,
}

pub fn load_refprot_proteome_map(path: &Path) -> HashMap<String, String> {
    let mut map = HashMap::new();
    if let Ok(text) = fs::read_to_string(path) {
        for line in text.lines() {
            let mut parts = line.split('\t');
            if let (Some(acc), Some(pid)) = (parts.next(), parts.next()) {
                map.insert(acc.to_string(), pid.to_string());
            }
        }
    }
    map
}

pub fn annotate_refprot_hits(
    grouped: &mut HashMap<String, Vec<diamond::DiamondHitRow>>,
    map: &HashMap<String, String>,
) {
    for hits in grouped.values_mut() {
        for hit in hits.iter_mut() {
            let proteome_id = map
                .get(&hit.sseqid)
                .cloned()
                .unwrap_or_else(|| "refprot".to_string());
            hit.source = HitSource::RefProt(proteome_id);
        }
    }
}

fn parse_evalue_to_f64(value: &str) -> f64 {
    if value.trim().is_empty() {
        return 1.0;
    }
    value.trim().parse::<f64>().unwrap_or_else(|_| {
        match value.trim().to_ascii_lowercase().as_str() {
            "inf" => f64::INFINITY,
            _ => 1.0,
        }
    })
}

pub fn filter_refprot_hits(
    grouped: &mut HashMap<String, Vec<diamond::DiamondHitRow>>,
    cfg: &RefProtFallbackConfig,
) {
    let mut empty_keys = Vec::new();
    for (gene, hits) in grouped.iter_mut() {
        hits.retain(|row| {
            if row.qcov < cfg.min_qcov {
                return false;
            }
            if row.scov < cfg.min_scov {
                return false;
            }
            if row.pident < cfg.min_pident {
                return false;
            }
            let eval = parse_evalue_to_f64(&row.evalue);
            if eval.is_nan() || eval > cfg.max_evalue {
                return false;
            }
            true
        });
        hits.sort_by(|a, b| b.bitscore.total_cmp(&a.bitscore));
        if hits.len() > cfg.max_hits {
            hits.truncate(cfg.max_hits);
        }
        if hits.is_empty() {
            empty_keys.push(gene.clone());
        }
    }
    for key in empty_keys {
        grouped.remove(&key);
    }
}

pub fn build_consensus_outputs(
    metrics: &[GeneMetrics],
    grouped: &HashMap<String, Vec<diamond::DiamondHitRow>>,
    refprot_grouped: &HashMap<String, Vec<diamond::DiamondHitRow>>,
    cons_cfg: &consensus::ConsensusConfig,
    refprot_fallback_cfg: &RefProtFallbackConfig,
) -> ConsensusBuildOutput {
    let mut panel_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut panel_stats_rows: Vec<(String, consensus::PanelStats)> = Vec::new();
    let mut panel_len_stats: HashMap<String, consensus::LenStats> = HashMap::new();
    let mut panel_agg_rows: Vec<PanelAggDebugRow> = Vec::new();
    let mut len_map: HashMap<String, LengthSummary> = HashMap::new();
    let mut panel_prov_map: HashMap<String, PanelProvenanceCounts> = HashMap::new();
    let mut panel_prov_rows: Vec<(String, PanelProvenanceCounts)> = Vec::new();
    let mut refprot_used_panels: u64 = 0;
    let mut taxonomy_hits_map: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
    let trigger_k = refprot_fallback_cfg.trigger_k.max(0);

    for m in metrics {
        let primary_hits = grouped.get(&m.gene_id);
        let fallback_hits = refprot_grouped.get(&m.gene_id);
        if primary_hits.is_some() || fallback_hits.is_some() {
            let mut panel_input: Vec<diamond::DiamondHitRow> =
                primary_hits.cloned().unwrap_or_default();
            if panel_input.len() < trigger_k {
                if let Some(extra) = fallback_hits {
                    panel_input.extend_from_slice(extra);
                }
            }
            let mut taxonomy_hits: Vec<diamond::DiamondHitRow> =
                primary_hits.cloned().unwrap_or_default();
            if let Some(extra) = fallback_hits {
                taxonomy_hits.extend_from_slice(extra);
            }
            taxonomy_hits.sort_by(|a, b| {
                b.bitscore
                    .partial_cmp(&a.bitscore)
                    .unwrap_or(Ordering::Equal)
            });
            taxonomy_hits_map.insert(m.gene_id.clone(), taxonomy_hits);
            if panel_input.is_empty() {
                panel_stats_rows.push((
                    m.gene_id.clone(),
                    consensus::PanelStats {
                        total_hits: 0,
                        ..Default::default()
                    },
                ));
                panel_prov_map.insert(m.gene_id.clone(), PanelProvenanceCounts::default());
                panel_prov_rows.push((m.gene_id.clone(), PanelProvenanceCounts::default()));
                continue;
            }
            let panel_res = consensus::select_panel_with_result(&panel_input, cons_cfg);
            let selection = panel_res.selection.clone();
            panel_len_stats.insert(m.gene_id.clone(), panel_res.len_stats.clone());
            let selected_set: HashSet<_> = selection.ids.iter().cloned().collect();
            let mut source_by_id: HashMap<String, HitSource> = HashMap::new();
            for agg in &panel_res.aggregated_hits {
                source_by_id
                    .entry(agg.sseqid.clone())
                    .or_insert_with(|| agg.source.clone());
                panel_agg_rows.push(PanelAggDebugRow {
                    gene_id: m.gene_id.clone(),
                    subject_id: agg.sseqid.clone(),
                    qcov: agg.qcov,
                    scov: agg.scov,
                    len_ratio: agg.len_ratio,
                    bitscore: agg.bitscore,
                    hit_count: agg.hit_count,
                    selected: selected_set.contains(&agg.sseqid),
                    source: agg.source.clone(),
                    quality: agg.quality,
                    diversity_key: agg.diversity_key(cons_cfg.diversity_rank_index),
                });
            }
            let mut prov_counts = PanelProvenanceCounts::default();
            for sid in &selection.ids {
                match source_by_id.get(sid) {
                    Some(HitSource::SwissProt) => prov_counts.swissprot += 1,
                    Some(HitSource::RefProt(_)) => prov_counts.refprot += 1,
                    Some(HitSource::Cluster) => prov_counts.cluster += 1,
                    None => prov_counts.cluster += 1,
                }
            }
            let panel_ids = selection.ids.clone();
            panel_stats_rows.push((m.gene_id.clone(), selection.stats.clone()));
            if !panel_ids.is_empty() {
                panel_map.insert(m.gene_id.clone(), panel_ids.clone());
                let id_set: HashSet<String> = panel_ids.iter().cloned().collect();
                let slens: Vec<usize> = panel_input
                    .iter()
                    .filter(|h| id_set.contains(&h.sseqid))
                    .map(|h| h.slen)
                    .filter(|&x| x > 0)
                    .collect();
                if slens.len() >= cons_cfg.min_hits {
                    if let Some(lc) = length::compute_length_consistency(m.length, &slens) {
                        len_map.insert(
                            m.gene_id.clone(),
                            (
                                lc.score,
                                lc.z,
                                lc.ratio,
                                lc.class_,
                                lc.expected_min,
                                lc.expected_max,
                                lc.in_expected_range,
                                lc.panel_n,
                            ),
                        );
                    }
                }
            }
            if prov_counts.refprot > 0 {
                refprot_used_panels += 1;
            }
            panel_prov_map.insert(m.gene_id.clone(), prov_counts.clone());
            panel_prov_rows.push((m.gene_id.clone(), prov_counts));
        } else {
            taxonomy_hits_map.insert(m.gene_id.clone(), Vec::new());
            panel_stats_rows.push((
                m.gene_id.clone(),
                consensus::PanelStats {
                    total_hits: 0,
                    ..Default::default()
                },
            ));
            panel_prov_map.insert(m.gene_id.clone(), PanelProvenanceCounts::default());
            panel_prov_rows.push((m.gene_id.clone(), PanelProvenanceCounts::default()));
        }
    }

    ConsensusBuildOutput {
        panel_map,
        panel_stats_rows,
        panel_len_stats,
        panel_agg_rows,
        len_map,
        panel_prov_map,
        panel_prov_rows,
        taxonomy_hits_map,
        refprot_used_panels,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_ref_hit(
        gene: &str,
        subj: &str,
        qcov: f64,
        scov: f64,
        pident: f64,
        evalue: &str,
        bitscore: f64,
    ) -> diamond::DiamondHitRow {
        diamond::DiamondHitRow {
            qseqid: gene.into(),
            sseqid: subj.into(),
            bitscore,
            evalue: evalue.into(),
            length: 120,
            qcov,
            scov,
            pident,
            qstart: 1,
            qend: 100,
            sstart: 5,
            send: 105,
            qlen: 100,
            slen: 110,
            source: HitSource::RefProt("P0001".into()),
            staxid: Some(1),
            lineage: vec![],
        }
    }

    #[test]
    fn filter_refprot_hits_respects_thresholds_and_cap() {
        let mut grouped: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
        grouped.insert(
            "gene1".into(),
            vec![
                build_ref_hit("gene1", "bad_qcov", 0.2, 0.4, 40.0, "1e-20", 210.0),
                build_ref_hit("gene1", "good1", 0.8, 0.6, 45.0, "1e-20", 260.0),
                build_ref_hit("gene1", "good2", 0.9, 0.7, 50.0, "1e-15", 240.0),
                build_ref_hit("gene1", "high_eval", 0.9, 0.7, 50.0, "1e-2", 230.0),
            ],
        );
        grouped.insert(
            "gene2".into(),
            vec![build_ref_hit(
                "gene2", "fail_all", 0.1, 0.1, 10.0, "1", 200.0,
            )],
        );
        let mut cfg = RefProtFallbackConfig {
            min_qcov: 0.5,
            min_scov: 0.4,
            min_pident: 30.0,
            max_evalue: 1e-10,
            max_hits: 1,
            ..Default::default()
        };
        cfg.trigger_k = 5;
        filter_refprot_hits(&mut grouped, &cfg);
        let kept = grouped.get("gene1").unwrap();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].sseqid, "good1");
        assert!(!grouped.contains_key("gene2"));
    }
}
