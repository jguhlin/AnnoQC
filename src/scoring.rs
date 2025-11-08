use crate::diamond::DiamondHitStats;
use crate::metrics::IntrinsicMetrics;

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

pub fn compute_intrinsic_score(im: &IntrinsicMetrics) -> f64 {
    let amb_pen = (1.0 - (im.ambiguous_fraction).min(1.0)).max(0.0);
    let lc_pen = (1.0 - im.low_complexity_fraction.min(1.0)).max(0.0);
    let hp_pen = (1.0 - (im.max_homopolymer as f64 / 30.0).min(1.0)).max(0.0);
    (0.5 * amb_pen + 0.4 * lc_pen + 0.1 * hp_pen).clamp(0.0, 1.0)
}

pub fn compute_taxonomy_score(resolved: bool) -> f64 {
    if resolved {
        1.0
    } else {
        0.0
    }
}
