use crate::config_types::ScoringConfigOverride;
use crate::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct WeightSet {
    pub(crate) homology: f64,
    pub(crate) intrinsic: f64,
    pub(crate) taxonomy: f64,
    pub(crate) domains: f64,
    pub(crate) domains_strength: f64,
    pub(crate) length: f64,
    pub(crate) orphan: f64,
    pub(crate) subject_cov: f64,
    pub(crate) termini: f64,
    pub(crate) divergence: f64,
    pub(crate) conserved_regions: f64,
    pub(crate) block_conservation: f64,
    pub(crate) genomic: f64,
    pub(crate) rnaseq: f64,
}

impl WeightSet {
    pub(crate) fn sum(&self) -> f64 {
        self.homology
            + self.intrinsic
            + self.taxonomy
            + self.domains
            + self.domains_strength
            + self.length
            + self.orphan
            + self.subject_cov
            + self.termini
            + self.divergence
            + self.conserved_regions
            + self.block_conservation
            + self.genomic
            + self.rnaseq
    }

    pub(crate) fn sum_available(&self, presence: &PillarPresence) -> f64 {
        let mut sum = 0.0;
        if presence.homology {
            sum += self.homology;
        }
        if presence.intrinsic {
            sum += self.intrinsic;
        }
        if presence.taxonomy {
            sum += self.taxonomy;
        }
        if presence.domains {
            sum += self.domains;
        }
        if presence.domains_strength {
            sum += self.domains_strength;
        }
        if presence.length {
            sum += self.length;
        }
        if presence.orphan {
            sum += self.orphan;
        }
        if presence.subject_cov {
            sum += self.subject_cov;
        }
        if presence.termini {
            sum += self.termini;
        }
        if presence.divergence {
            sum += self.divergence;
        }
        if presence.conserved_regions {
            sum += self.conserved_regions;
        }
        if presence.block_conservation {
            sum += self.block_conservation;
        }
        if presence.genomic {
            sum += self.genomic;
        }
        if presence.rnaseq {
            sum += self.rnaseq;
        }
        sum
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct PillarPresence {
    pub(crate) homology: bool,
    pub(crate) intrinsic: bool,
    pub(crate) taxonomy: bool,
    pub(crate) domains: bool,
    pub(crate) domains_strength: bool,
    pub(crate) length: bool,
    pub(crate) orphan: bool,
    pub(crate) subject_cov: bool,
    pub(crate) termini: bool,
    pub(crate) divergence: bool,
    pub(crate) conserved_regions: bool,
    pub(crate) block_conservation: bool,
    pub(crate) genomic: bool,
    pub(crate) rnaseq: bool,
}

impl PillarPresence {
    pub(crate) fn missing_with_weights(&self, weights: &WeightSet) -> Vec<&'static str> {
        let mut missing = Vec::new();
        if weights.homology > 0.0 && !self.homology {
            missing.push("homology");
        }
        if weights.intrinsic > 0.0 && !self.intrinsic {
            missing.push("intrinsic");
        }
        if weights.taxonomy > 0.0 && !self.taxonomy {
            missing.push("taxonomy");
        }
        if weights.domains > 0.0 && !self.domains {
            missing.push("domains");
        }
        if weights.domains_strength > 0.0 && !self.domains_strength {
            missing.push("domains_strength");
        }
        if weights.length > 0.0 && !self.length {
            missing.push("length");
        }
        if weights.orphan > 0.0 && !self.orphan {
            missing.push("orphan");
        }
        if weights.subject_cov > 0.0 && !self.subject_cov {
            missing.push("subject_cov");
        }
        if weights.termini > 0.0 && !self.termini {
            missing.push("termini");
        }
        if weights.divergence > 0.0 && !self.divergence {
            missing.push("divergence");
        }
        if weights.conserved_regions > 0.0 && !self.conserved_regions {
            missing.push("conserved_regions");
        }
        if weights.block_conservation > 0.0 && !self.block_conservation {
            missing.push("block_conservation");
        }
        if weights.genomic > 0.0 && !self.genomic {
            missing.push("genomic");
        }
        if weights.rnaseq > 0.0 && !self.rnaseq {
            missing.push("rnaseq");
        }
        missing
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ScoreTemp {
    pub(crate) gene_id: String,
    pub(crate) raw_score: f64,
    pub(crate) final_score: f64,
    pub(crate) available_weight: f64,
    pub(crate) missing: Vec<&'static str>,
}

pub(crate) fn scoring_thresholds(scoring: &Option<ScoringConfigOverride>) -> (f64, f64) {
    if let Some(cfg) = scoring {
        let h = cfg.thresholds.as_ref().and_then(|t| t.high).unwrap_or(0.8);
        let m = cfg
            .thresholds
            .as_ref()
            .and_then(|t| t.medium)
            .unwrap_or(0.5);
        (h, m)
    } else {
        (0.8, 0.5)
    }
}

pub(crate) fn scoring_weights(scoring: &Option<ScoringConfigOverride>) -> WeightSet {
    let mut ws = WeightSet {
        homology: 0.55,
        intrinsic: 0.3,
        taxonomy: 0.0,
        domains: 0.0,
        domains_strength: 0.05,
        length: 0.0,
        orphan: 0.0,
        subject_cov: 0.0,
        termini: 0.0,
        divergence: 0.0,
        conserved_regions: 0.1,
        block_conservation: 0.0,
        genomic: 0.0,
        rnaseq: 0.0,
    };
    if let Some(cfg) = scoring {
        if let Some(v) = cfg.weights.get("homology") {
            ws.homology = *v;
        }
        if let Some(v) = cfg.weights.get("intrinsic") {
            ws.intrinsic = *v;
        }
        if let Some(v) = cfg.weights.get("taxonomy") {
            ws.taxonomy = *v;
        }
        if let Some(v) = cfg.weights.get("domains") {
            ws.domains = *v;
        }
        if let Some(v) = cfg.weights.get("domains_strength") {
            ws.domains_strength = *v;
        }
        if let Some(v) = cfg.weights.get("length") {
            ws.length = *v;
        }
        if let Some(v) = cfg.weights.get("orphan") {
            ws.orphan = *v;
        }
        if let Some(v) = cfg.weights.get("subject_cov") {
            ws.subject_cov = *v;
        }
        if let Some(v) = cfg.weights.get("termini") {
            ws.termini = *v;
        }
        if let Some(v) = cfg.weights.get("divergence") {
            ws.divergence = *v;
        }
        if let Some(v) = cfg.weights.get("conserved_regions") {
            ws.conserved_regions = *v;
        }
        if let Some(v) = cfg.weights.get("block_conservation") {
            ws.block_conservation = *v;
        }
        if let Some(v) = cfg.weights.get("genomic") {
            ws.genomic = *v;
        }
        if let Some(v) = cfg.weights.get("rnaseq") {
            ws.rnaseq = *v;
        }
    }
    ws
}

pub(crate) fn format_classification(
    base: &str,
    missing: &[&'static str],
    classify_no_data: bool,
) -> String {
    if !classify_no_data || missing.is_empty() {
        return base.to_string();
    }
    let list = missing.join("|");
    format!("{} (no_data:{})", base, list)
}

pub(crate) fn classification_base(label: &str) -> &str {
    label.split([' ', '(', '[']).next().unwrap_or(label)
}

pub(crate) fn print_scoring_rubric(
    scoring: &Option<ScoringConfigOverride>,
    calibration: CalibrationSettings,
    classify_no_data: bool,
) {
    let weights = scoring_weights(scoring);
    let (th_high, th_med) = scoring_thresholds(scoring);
    let sum = weights.sum().max(1e-9);
    let norm = |v: f64| v / sum;
    println!("scoring_rubric");
    println!(
        "weights_raw: homology={:.3} intrinsic={:.3} taxonomy={:.3} domains={:.3} domains_strength={:.3} length={:.3} orphan={:.3} subject_cov={:.3} termini={:.3} divergence={:.3} conserved_regions={:.3} genomic={:.3} sum={:.3}",
        weights.homology,
        weights.intrinsic,
        weights.taxonomy,
        weights.domains,
        weights.domains_strength,
        weights.length,
        weights.orphan,
        weights.subject_cov,
        weights.termini,
        weights.divergence,
        weights.conserved_regions,
        weights.genomic,
        weights.sum(),
    );
    println!(
        "weights_normalized: homology={:.3} intrinsic={:.3} taxonomy={:.3} domains={:.3} domains_strength={:.3} length={:.3} orphan={:.3} subject_cov={:.3} termini={:.3} divergence={:.3} conserved_regions={:.3} genomic={:.3}",
        norm(weights.homology),
        norm(weights.intrinsic),
        norm(weights.taxonomy),
        norm(weights.domains),
        norm(weights.domains_strength),
        norm(weights.length),
        norm(weights.orphan),
        norm(weights.subject_cov),
        norm(weights.termini),
        norm(weights.divergence),
        norm(weights.conserved_regions),
        norm(weights.genomic),
    );
    println!("thresholds: high={:.3} medium={:.3}", th_high, th_med);
    println!(
        "calibration: mode={:?} min_samples={} min_unique={}",
        calibration.mode, calibration.min_samples, calibration.min_unique
    );
    println!(
        "no_data_handling: classify_no_data={} missing_pillars_excluded_from_weights=true",
        classify_no_data
    );
}
