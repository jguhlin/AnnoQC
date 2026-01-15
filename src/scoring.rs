use crate::diamond::DiamondHitStats;
use crate::genomic::GenomicMetrics;
use crate::hmmer::HmmscanSummary;
use crate::mafft::AlignmentMetrics;
use crate::metrics::IntrinsicMetrics;
use crate::structvar::StructVar;
use crate::taxonomy::TaxonomyEvidence;

pub fn compute_homology_score(s: Option<&DiamondHitStats>) -> f64 {
    if let Some(s) = s {
        let hits = (s.count as f64 / 10.0).min(1.0);
        let density = if s.top_len > 0 {
            (s.top_bitscore / s.top_len as f64) / 5.0
        } else {
            0.0
        };
        let density = density.min(1.0);
        let cov_quality = ((s.top_qcov + s.top_scov) / 2.0).min(1.0);
        let cov_penalty = 1.0 - (s.coverage_delta).min(1.0);
        let raw = 0.3 * hits + 0.3 * density + 0.3 * cov_quality + 0.1 * cov_penalty;
        raw.clamp(0.0, 1.0)
    } else {
        0.0
    }
}

/// Computes intrinsic sequence quality from ambiguity, low complexity, and homopolymer metrics.
/// Returns a score in [0.0, 1.0].
pub fn compute_intrinsic_score(im: &IntrinsicMetrics) -> f64 {
    let amb_pen = (1.0 - (im.ambiguous_fraction).min(1.0)).max(0.0);
    let lc_pen = (1.0 - im.low_complexity_fraction.min(1.0)).max(0.0);
    let hp_pen = (1.0 - (im.max_homopolymer as f64 / 30.0).min(1.0)).max(0.0);
    (0.5 * amb_pen + 0.4 * lc_pen + 0.1 * hp_pen).clamp(0.0, 1.0)
}

/// Computes taxonomy congruence score from available evidence.
/// Returns a score in [0.0, 1.0], treating non-finite inputs as 0.0.
pub fn compute_taxonomy_score(evidence: Option<&TaxonomyEvidence>) -> f64 {
    evidence
        .map(|e| e.congruence_score)
        .filter(|score| score.is_finite())
        .map(|score| score.clamp(0.0, 1.0))
        .unwrap_or(0.0)
}

pub fn compute_subject_cov_score(s: Option<&DiamondHitStats>) -> f64 {
    s.map(|s| s.top_scov.clamp(0.0, 1.0)).unwrap_or(0.0)
}

pub fn compute_subject_cov_penalty(s: Option<&DiamondHitStats>) -> f64 {
    1.0 - compute_subject_cov_score(s)
}

/// Computes a domain-strength score from the best (lowest) hmmscan e-value.
///
/// Returns a score in [0.0, 1.0]. Missing or hitless hmmscan results yield 0.0.
pub fn compute_domains_strength_score(summary: Option<&HmmscanSummary>) -> f64 {
    let Some(s) = summary else {
        return 0.0;
    };
    if s.hits_count == 0 {
        return 0.0;
    }
    let Some(ev) = s.top_evalue else {
        return 0.0;
    };
    let le = if ev > 0.0 { -ev.log10() } else { 100.0 };
    (le / 20.0).clamp(0.0, 1.0)
}

pub fn compute_divergence_score(am: Option<&crate::mafft::AlignmentMetrics>) -> f64 {
    if let Some(am) = am {
        // divergence_ratio = query_id / panel_id.
        // If ratio >= 0.8, it's good (1.0).
        // If ratio < 0.5, it's bad (0.0).
        // Linear ramp in between.
        let r = am.divergence_ratio;
        if r >= 0.8 {
            1.0
        } else if r < 0.5 {
            0.0
        } else {
            (r - 0.5) / 0.3
        }
    } else {
        0.0
    }
}

/// Computes a conserved-regions / alignment-integrity score from alignment metrics.
///
/// This is intended to capture "internal" alignment quality (conserved core, gaps, large runs),
/// distinct from endpoints (termini) and divergence.
///
/// Returns a score in [0.0, 1.0]. Treats missing alignment metrics as 0.0.
pub fn compute_conserved_regions_score(am: Option<&AlignmentMetrics>) -> f64 {
    let Some(am) = am else {
        return 0.0;
    };

    let conserved = am.conserved_fraction.clamp(0.0, 1.0);

    // Penalize heavy gappiness in the query row. This is a simple, tunable mapping.
    // By default: gaps >= 20% => 0.0 contribution from this term.
    let gap_pen = (1.0 - (am.query_gap_fraction / 0.20).clamp(0.0, 1.0)).clamp(0.0, 1.0);

    // Penalize long contiguous gap runs.
    let run_pen = (-(am.max_gap_run as f64 / 25.0)).exp().clamp(0.0, 1.0);

    // Mild penalties for strong "missing exon" / "retained intron" signatures.
    let missing_pen = (-(am.missing_exon_run as f64 / 20.0)).exp().clamp(0.0, 1.0);
    let intron_pen = (-(am.retained_intron_run as f64 / 20.0))
        .exp()
        .clamp(0.0, 1.0);

    (0.55 * conserved + 0.20 * gap_pen + 0.15 * run_pen + 0.05 * missing_pen + 0.05 * intron_pen)
        .clamp(0.0, 1.0)
}

/// Computes a multiplicative penalty factor for structural-variation / "different genes" signals.
///
/// This is derived from DIAMOND HSP layouts and is meant to be a *major detriment* when a fusion
/// or split is suspected. The factor is applied as `score *= structvar_multiplier` before
/// calibration.
///
/// Returns a value in [0.0, 1.0]. Missing structvar evidence defaults to 1.0 (no penalty).
pub fn compute_structvar_multiplier(sv: Option<&StructVar>) -> f64 {
    let Some(sv) = sv else {
        return 1.0;
    };
    match sv.classification.as_str() {
        "FusionPossible" => 0.05,
        "SplitPossible" => 0.20,
        "InternalDuplicationPossible" => 0.70,
        _ => 1.0,
    }
}

pub fn compute_genomic_score(gm: Option<&GenomicMetrics>) -> f64 {
    compute_genomic_score_with_cfg(gm, None, None, None)
}

pub fn compute_genomic_score_with_cfg(
    gm: Option<&GenomicMetrics>,
    min_canonical: Option<f64>,
    max_noncanonical: Option<f64>,
    max_weird: Option<f64>,
) -> f64 {
    let Some(gm) = gm else {
        return 0.0;
    };
    if gm.introns_total == 0 {
        return 1.0;
    }
    let total = gm.introns_total as f64;
    let canonical = gm.splice_canonical as f64 / total;
    let noncanonical = (gm.splice_major_noncan + gm.splice_minor) as f64 / total;
    let weird = gm.splice_weird as f64 / total;
    let mut score = 0.7 * canonical + 0.2 * (1.0 - noncanonical) + 0.1 * (1.0 - weird);
    if let Some(min) = min_canonical {
        if min > 0.0 && canonical < min {
            score *= (canonical / min).clamp(0.0, 1.0);
        }
    }
    if let Some(max) = max_noncanonical {
        if max > 0.0 && noncanonical > max {
            score *= (max / noncanonical).clamp(0.0, 1.0);
        }
    }
    if let Some(max) = max_weird {
        if max > 0.0 && weird > max {
            score *= (max / weird).clamp(0.0, 1.0);
        }
    }
    score.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::{
        compute_conserved_regions_score, compute_domains_strength_score,
        compute_genomic_score_with_cfg, compute_structvar_multiplier, compute_taxonomy_score,
    };
    use crate::genomic::GenomicMetrics;
    use crate::hmmer::HmmscanSummary;
    use crate::mafft::AlignmentMetrics;
    use crate::structvar::StructVar;
    use crate::taxonomy::TaxonomyEvidence;

    #[test]
    fn genomic_score_respects_thresholds() {
        let gm = GenomicMetrics {
            introns_total: 4,
            splice_canonical: 1,
            splice_major_noncan: 2,
            splice_minor: 1,
            splice_weird: 0,
            ..Default::default()
        };
        let base = compute_genomic_score_with_cfg(Some(&gm), None, None, None);
        let tightened = compute_genomic_score_with_cfg(Some(&gm), Some(0.5), Some(0.2), None);
        assert!(tightened <= base);
    }

    #[test]
    fn taxonomy_score_ignores_non_finite_values() {
        let mut ev = TaxonomyEvidence::default();
        ev.congruence_score = f64::NAN;
        assert_eq!(compute_taxonomy_score(Some(&ev)), 0.0);
        ev.congruence_score = 1.2;
        assert_eq!(compute_taxonomy_score(Some(&ev)), 1.0);
    }

    #[test]
    fn conserved_regions_score_behaves_reasonably() {
        let mut m = AlignmentMetrics::default();
        m.conserved_fraction = 0.9;
        m.query_gap_fraction = 0.01;
        m.max_gap_run = 1;
        m.missing_exon_run = 0;
        m.retained_intron_run = 0;
        let good = compute_conserved_regions_score(Some(&m));
        assert!(good > 0.7);

        m.query_gap_fraction = 0.4;
        let gappy = compute_conserved_regions_score(Some(&m));
        assert!(gappy < good);

        m.query_gap_fraction = 0.01;
        m.max_gap_run = 80;
        let long_run = compute_conserved_regions_score(Some(&m));
        assert!(long_run < good);

        m.max_gap_run = 1;
        m.missing_exon_run = 50;
        let missing_exon = compute_conserved_regions_score(Some(&m));
        assert!(missing_exon < good);
    }

    #[test]
    fn structvar_multiplier_matches_classes() {
        let mut sv = StructVar::default();
        sv.classification = "None".into();
        assert_eq!(compute_structvar_multiplier(Some(&sv)), 1.0);
        sv.classification = "InternalDuplicationPossible".into();
        assert_eq!(compute_structvar_multiplier(Some(&sv)), 0.70);
        sv.classification = "SplitPossible".into();
        assert_eq!(compute_structvar_multiplier(Some(&sv)), 0.20);
        sv.classification = "FusionPossible".into();
        assert_eq!(compute_structvar_multiplier(Some(&sv)), 0.05);
    }

    #[test]
    fn domains_strength_score_maps_evalue() {
        let mut s = HmmscanSummary::default();
        s.hits_count = 0;
        s.top_evalue = None;
        assert_eq!(compute_domains_strength_score(Some(&s)), 0.0);

        s.hits_count = 1;
        s.top_evalue = Some(1e-40);
        let strong = compute_domains_strength_score(Some(&s));
        assert!(strong > 0.9);

        s.top_evalue = Some(1e-2);
        let weak = compute_domains_strength_score(Some(&s));
        assert!(weak < strong);
    }
}
