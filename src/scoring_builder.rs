use std::collections::HashMap;

use crate::config_types::ScoringConfigOverride;
use crate::consensus_support::LengthSummary;
use crate::rnaseq::RnaseqMetrics;
use crate::scoring_thresholds::{
    format_classification, scoring_thresholds, scoring_weights, PillarPresence, ScoreTemp,
};
use crate::taxonomy::TaxonomyEvidence;
use crate::*;

use crate::scoring_calibration::{
    apply_isotonic_calibration, apply_percentile_calibration, calibration_has_min_samples,
};

#[derive(Clone, Debug, Default)]
pub(crate) struct ComponentScores {
    pub(crate) homology: f64,
    pub(crate) intrinsic: f64,
    pub(crate) taxonomy: Option<f64>,
    pub(crate) orphan: f64,
    pub(crate) subject_cov: f64,
    pub(crate) termini: Option<f64>,
    pub(crate) divergence: Option<f64>,
    pub(crate) conserved_regions: Option<f64>,
    pub(crate) genomic: f64,
    pub(crate) rnaseq: f64,
}

fn adjust_homology_score(
    raw: f64,
    divergence_score: f64,
    divergence_present: bool,
    length_score: f64,
    length_present: bool,
) -> f64 {
    let mut score = raw;
    if divergence_present {
        let adj = (0.5 + 0.5 * divergence_score.clamp(0.0, 1.0)).clamp(0.0, 1.0);
        score *= adj;
    }
    if length_present {
        let adj = (0.5 + 0.5 * length_score.clamp(0.0, 1.0)).clamp(0.0, 1.0);
        score *= adj;
    }
    score.clamp(0.0, 1.0)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_scores_map(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scoring: &Option<ScoringConfigOverride>,
    taxonomy_enabled: bool,
    hmmer_enabled: bool,
    hmmsum_map: Option<&HashMap<String, hmmer::HmmscanSummary>>,
    arch_map: Option<&HashMap<String, f64>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
    alignment_map: Option<&HashMap<String, mafft::AlignmentMetrics>>,
    structvar_map: Option<&HashMap<String, structvar::StructVar>>,
    genomic_map: Option<&HashMap<String, genomic::GenomicMetrics>>,
    rnaseq_map: Option<&HashMap<String, RnaseqMetrics>>,
    rnaseq_enabled: bool,
    calibration: CalibrationSettings,
    classify_no_data: bool,
) -> (HashMap<String, (f64, String)>, HashMap<String, f64>) {
    let mut weights = scoring_weights(scoring);
    if !hmmer_enabled {
        if weights.domains > 0.0 || weights.domains_strength > 0.0 || weights.orphan > 0.0 {
            log::warn!("hmmer disabled; ignoring domains/domains_strength/orphan weights");
        }
        weights.domains = 0.0;
        weights.domains_strength = 0.0;
        weights.orphan = 0.0;
    }
    let (th_high, th_med) = scoring_thresholds(scoring);
    let genomic_cfg = scoring.as_ref().and_then(|s| s.genomic.as_ref());
    let caps_cfg = scoring.as_ref().and_then(|s| s.caps.as_ref());
    let mut staging = Vec::with_capacity(metrics.len());
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let homology_present = s.map(|s| s.count > 0).unwrap_or(false);
        let intrinsic_present = intrinsic.contains_key(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h_raw = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let (t, taxonomy_present) = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            let present = evidence
                .map(|ev| ev.considered > 0 || ev.top_hit.is_some())
                .unwrap_or(false);
            (compute_taxonomy_score(evidence), present)
        } else {
            (0.0, false)
        };
        let d = arch_map
            .and_then(|am| am.get(&m.gene_id))
            .copied()
            .unwrap_or(0.0);
        let domains_present = arch_map
            .map(|am| am.contains_key(&m.gene_id))
            .unwrap_or(false);
        let domains_strength_score =
            compute_domains_strength_score(hmmsum_map.and_then(|hm| hm.get(&m.gene_id)));
        let domains_strength_present = hmmsum_map
            .map(|hm| hm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let l = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let length_present = len_map
            .map(|lm| lm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let subject_cov_score = compute_subject_cov_score(s);
        let subject_cov_present = homology_present;
        let o = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        let orphan_present = orphan_map
            .map(|om| om.contains_key(&m.gene_id))
            .unwrap_or(false);
        let align_entry = alignment_map.and_then(|am| am.get(&m.gene_id));
        let termini_present = align_entry
            .map(|a| a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI)
            .unwrap_or(false);
        let termini_score = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            })
            .unwrap_or(0.0);
        let divergence_present = align_entry
            .map(|a| a.sequences_aligned > 0)
            .unwrap_or(false);
        let divergence_score =
            scoring::compute_divergence_score(alignment_map.and_then(|am| am.get(&m.gene_id)));
        let conserved_regions_present = align_entry
            .map(|a| {
                a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
            })
            .unwrap_or(false);
        let conserved_regions_score =
            compute_conserved_regions_score(alignment_map.and_then(|am| am.get(&m.gene_id)));
        let block_conservation_present = align_entry
            .map(|a| a.mafft_enabled && !a.conserved_blocks.is_empty())
            .unwrap_or(false);
        let block_conservation_score = scoring::compute_block_conservation_score(
            alignment_map.and_then(|am| am.get(&m.gene_id)),
        );
        let genomic_score = compute_genomic_score_with_cfg(
            genomic_map.and_then(|gm| gm.get(&m.gene_id)),
            genomic_cfg.and_then(|g| g.min_canonical),
            genomic_cfg.and_then(|g| g.max_noncanonical),
            genomic_cfg.and_then(|g| g.max_weird),
        );
        let genomic_present = genomic_map
            .map(|gm| gm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let (rnaseq_score, rnaseq_present) = if rnaseq_enabled {
            let evidence = rnaseq_map.and_then(|rm| rm.get(&m.gene_id));
            let score = compute_rnaseq_score(evidence);
            let present = rnaseq_map.map_or(false, |rm| rm.contains_key(&m.gene_id));
            (score, present)
        } else {
            (0.0, false)
        };
        let h = adjust_homology_score(
            h_raw,
            divergence_score,
            divergence_present,
            l,
            length_present,
        );
        let presence = PillarPresence {
            homology: homology_present,
            intrinsic: intrinsic_present,
            taxonomy: taxonomy_present,
            domains: domains_present,
            domains_strength: domains_strength_present,
            length: length_present,
            orphan: orphan_present,
            subject_cov: subject_cov_present,
            termini: termini_present,
            divergence: divergence_present,
            conserved_regions: conserved_regions_present,
            block_conservation: block_conservation_present,
            genomic: genomic_present,
            rnaseq: rnaseq_present,
        };
        let missing = presence.missing_with_weights(&weights);
        let available_weight = weights.sum_available(&presence);
        let mut total_weight = available_weight;
        if !presence.homology {
            total_weight += weights.homology;
        }
        if !presence.domains {
            total_weight += weights.domains;
        }
        if hmmer_enabled && !presence.domains_strength {
            total_weight += weights.domains_strength;
        }
        let mut numerator = 0.0;
        if presence.homology {
            numerator += weights.homology * h;
        }
        if presence.intrinsic {
            numerator += weights.intrinsic * i;
        }
        if presence.taxonomy {
            numerator += weights.taxonomy * t;
        }
        if presence.domains {
            numerator += weights.domains * d;
        }
        if presence.domains_strength {
            numerator += weights.domains_strength * domains_strength_score;
        }
        if presence.length {
            numerator += weights.length * l;
        }
        if presence.orphan {
            numerator += weights.orphan * o;
        }
        if presence.subject_cov {
            numerator += weights.subject_cov * subject_cov_score;
        }
        if presence.termini {
            numerator += weights.termini * termini_score;
        }
        if presence.divergence {
            numerator += weights.divergence * divergence_score;
        }
        if presence.conserved_regions {
            numerator += weights.conserved_regions * conserved_regions_score;
        }
        if presence.block_conservation {
            numerator += weights.block_conservation * block_conservation_score;
        }
        if presence.genomic {
            numerator += weights.genomic * genomic_score;
        }
        if presence.rnaseq {
            numerator += weights.rnaseq * rnaseq_score;
        }
        let score = if total_weight > 1e-6 {
            (numerator / total_weight).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let sv = structvar_map.and_then(|sv| sv.get(&m.gene_id));
        let structvar_multiplier = if homology_present {
            compute_structvar_multiplier(sv)
        } else {
            1.0
        };
        let mut score = (score * structvar_multiplier).clamp(0.0, 1.0);
        if let (Some(cfg), Some(sv)) = (caps_cfg, sv) {
            let cap = match sv.classification.as_str() {
                "FusionPossible" => cfg.structvar_fusion_max,
                "SplitPossible" => cfg.structvar_split_max,
                "InternalDuplicationPossible" => cfg.structvar_dup_max,
                _ => None,
            };
            if let Some(max_score) = cap {
                score = score.min(max_score.clamp(0.0, 1.0));
            }
        }
        staging.push(ScoreTemp {
            gene_id: m.gene_id.clone(),
            raw_score: score,
            final_score: score,
            available_weight,
            missing,
        });
    }
    if !matches!(calibration.mode, CalibrationMode::Off) {
        if calibration_has_min_samples(&staging, calibration.min_samples, calibration.min_unique) {
            match calibration.mode {
                CalibrationMode::Percentile => apply_percentile_calibration(&mut staging),
                CalibrationMode::Isotonic => apply_isotonic_calibration(&mut staging),
                CalibrationMode::Off => {}
            }
        } else {
            let mut unique: std::collections::HashSet<u64> = std::collections::HashSet::new();
            for e in &staging {
                unique.insert(e.raw_score.to_bits());
            }
            log::warn!(
                "calibration skipped: samples={} unique={} (min_samples={}, min_unique={})",
                staging.len(),
                unique.len(),
                calibration.min_samples,
                calibration.min_unique
            );
        }
    }
    let mut out = HashMap::new();
    let mut raw_map = HashMap::new();
    for entry in staging {
        let base = if entry.available_weight < 1e-6 {
            "NoData"
        } else if entry.final_score >= th_high {
            "High"
        } else if entry.final_score >= th_med {
            "Medium"
        } else {
            "Low"
        };
        let classif = format_classification(base, &entry.missing, classify_no_data);
        raw_map.insert(entry.gene_id.clone(), entry.raw_score);
        out.insert(entry.gene_id, (entry.final_score, classif));
    }
    (out, raw_map)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_component_scores(
    metrics: &[GeneMetrics],
    stats: &HashMap<String, diamond::DiamondHitStats>,
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    taxonomy_enabled: bool,
    orphan_map: Option<&HashMap<String, hmmer::OrphanAnalysis>>,
    taxonomy_map: Option<&HashMap<String, Option<TaxonomyEvidence>>>,
    alignment_map: Option<&HashMap<String, mafft::AlignmentMetrics>>,
    len_map: Option<&HashMap<String, LengthSummary>>,
    genomic_map: Option<&HashMap<String, genomic::GenomicMetrics>>,
    rnaseq_map: Option<&HashMap<String, RnaseqMetrics>>,
    rnaseq_enabled: bool,
    scoring: &Option<ScoringConfigOverride>,
) -> HashMap<String, ComponentScores> {
    let genomic_cfg = scoring.as_ref().and_then(|s| s.genomic.as_ref());
    let mut out = HashMap::new();
    for m in metrics {
        let s = stats.get(&m.gene_id);
        let default_im = metrics::IntrinsicMetrics::default();
        let im = intrinsic
            .get(&m.gene_id)
            .map(|t| &t.0)
            .unwrap_or(&default_im);
        let h_raw = compute_homology_score(s);
        let i = compute_intrinsic_score(im);
        let taxonomy_component = if taxonomy_enabled {
            let evidence = taxonomy_map
                .and_then(|tm| tm.get(&m.gene_id))
                .and_then(|opt| opt.as_ref());
            if evidence
                .map(|ev| ev.considered > 0 || ev.top_hit.is_some())
                .unwrap_or(false)
            {
                Some(compute_taxonomy_score(evidence))
            } else {
                None
            }
        } else {
            None
        };
        let orphan_component = orphan_map
            .and_then(|om| om.get(&m.gene_id))
            .map(|oa| oa.score)
            .unwrap_or(1.0);
        let subject_cov_score = compute_subject_cov_score(s);
        let termini_component = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .and_then(|a| {
                if a.mafft_enabled && a.sequences_aligned >= mafft::MIN_PANEL_FOR_TERMINI {
                    Some((a.start_concordance + a.end_concordance) / 2.0)
                } else {
                    None
                }
            });
        let conserved_regions_component =
            alignment_map
                .and_then(|am| am.get(&m.gene_id))
                .and_then(|a| {
                    if a.mafft_enabled
                        && a.sequences_aligned >= mafft::MIN_PANEL_FOR_CONSERVED_REGIONS
                    {
                        Some(compute_conserved_regions_score(Some(a)))
                    } else {
                        None
                    }
                });
        let divergence_component = alignment_map
            .and_then(|am| am.get(&m.gene_id))
            .map(|a| scoring::compute_divergence_score(Some(a)));
        let divergence_present = divergence_component.is_some();
        let length_score = len_map
            .and_then(|lm| lm.get(&m.gene_id).map(|t| t.0))
            .unwrap_or(0.0);
        let length_present = len_map
            .map(|lm| lm.contains_key(&m.gene_id))
            .unwrap_or(false);
        let h = adjust_homology_score(
            h_raw,
            divergence_component.unwrap_or(0.0),
            divergence_present,
            length_score,
            length_present,
        );
        let genomic_component = compute_genomic_score_with_cfg(
            genomic_map.and_then(|gm| gm.get(&m.gene_id)),
            genomic_cfg.and_then(|g| g.min_canonical),
            genomic_cfg.and_then(|g| g.max_noncanonical),
            genomic_cfg.and_then(|g| g.max_weird),
        );
        let rnaseq_component = if rnaseq_enabled {
            compute_rnaseq_score(rnaseq_map.and_then(|rm| rm.get(&m.gene_id)))
        } else {
            0.0
        };
        out.insert(
            m.gene_id.clone(),
            ComponentScores {
                homology: h,
                intrinsic: i,
                taxonomy: taxonomy_component,
                orphan: orphan_component,
                subject_cov: subject_cov_score,
                termini: termini_component,
                divergence: divergence_component,
                conserved_regions: conserved_regions_component,
                genomic: genomic_component,
                rnaseq: rnaseq_component,
            },
        );
    }
    out
}
