use crate::diamond::DiamondHitStats;
use crate::genomic::GenomicMetrics;
use crate::metrics::IntrinsicMetrics;
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
    use super::{compute_genomic_score_with_cfg, compute_taxonomy_score};
    use crate::genomic::GenomicMetrics;
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
}
