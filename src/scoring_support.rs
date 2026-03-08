use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::*;

#[derive(Clone, Copy, Debug)]
struct WeightSet {
    homology: f64,
    intrinsic: f64,
    taxonomy: f64,
    domains: f64,
    domains_strength: f64,
    length: f64,
    orphan: f64,
    subject_cov: f64,
    termini: f64,
    divergence: f64,
    conserved_regions: f64,
    block_conservation: f64,
    genomic: f64,
    rnaseq: f64,
}

impl WeightSet {
    fn sum(&self) -> f64 {
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

    fn sum_available(&self, presence: &PillarPresence) -> f64 {
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
struct PillarPresence {
    homology: bool,
    intrinsic: bool,
    taxonomy: bool,
    domains: bool,
    domains_strength: bool,
    length: bool,
    orphan: bool,
    subject_cov: bool,
    termini: bool,
    divergence: bool,
    conserved_regions: bool,
    block_conservation: bool,
    genomic: bool,
    rnaseq: bool,
}

impl PillarPresence {
    fn missing_with_weights(&self, weights: &WeightSet) -> Vec<&'static str> {
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

#[derive(Clone, Debug)]
struct ScoreTemp {
    gene_id: String,
    raw_score: f64,
    final_score: f64,
    available_weight: f64,
    missing: Vec<&'static str>,
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

fn scoring_weights(scoring: &Option<ScoringConfigOverride>) -> WeightSet {
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

fn format_classification(base: &str, missing: &[&'static str], classify_no_data: bool) -> String {
    if !classify_no_data || missing.is_empty() {
        return base.to_string();
    }
    let list = missing.join("|");
    format!("{} (no_data:{})", base, list)
}

fn classification_base(label: &str) -> &str {
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

fn apply_percentile_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let len = order.len();
    for (rank, idx) in order.into_iter().enumerate() {
        let percentile = if len > 1 {
            (rank as f64 + 0.5) / len as f64
        } else {
            1.0
        };
        entries[idx].final_score = percentile.clamp(0.0, 1.0);
    }
}

fn apply_isotonic_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let n = order.len();
    #[derive(Clone)]
    struct Block {
        start: usize,
        end: usize,
        sum: f64,
        weight: usize,
    }
    let mut blocks: Vec<Block> = Vec::new();
    for (rank, _idx) in order.iter().enumerate() {
        let y = if n > 1 {
            (rank as f64 + 0.5) / n as f64
        } else {
            1.0
        };
        blocks.push(Block {
            start: rank,
            end: rank,
            sum: y,
            weight: 1,
        });
        while blocks.len() >= 2 {
            let k = blocks.len() - 1;
            let prev = &blocks[k - 1];
            let curr = &blocks[k];
            let avg_prev = prev.sum / prev.weight as f64;
            let avg_curr = curr.sum / curr.weight as f64;
            if avg_prev <= avg_curr {
                break;
            }
            let merged = Block {
                start: prev.start,
                end: curr.end,
                sum: prev.sum + curr.sum,
                weight: prev.weight + curr.weight,
            };
            blocks.pop();
            blocks.pop();
            blocks.push(merged);
        }
    }
    let mut fitted = vec![0.0; n];
    for block in blocks {
        let avg = (block.sum / block.weight as f64).clamp(0.0, 1.0);
        for value in fitted
            .iter_mut()
            .take(block.end.saturating_add(1))
            .skip(block.start)
        {
            *value = avg;
        }
    }
    for (rank, idx) in order.into_iter().enumerate() {
        entries[idx].final_score = fitted[rank];
    }
}

fn calibration_has_min_samples(
    entries: &[ScoreTemp],
    min_samples: usize,
    min_unique: usize,
) -> bool {
    if entries.len() < min_samples {
        return false;
    }
    let mut unique: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for e in entries {
        unique.insert(e.raw_score.to_bits());
        if unique.len() >= min_unique {
            return true;
        }
    }
    false
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

pub(crate) fn export_high_sequences(
    out_dir: &str,
    override_path: Option<&str>,
    metrics: &[GeneMetrics],
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scores_map: &HashMap<String, (f64, String)>,
) -> Result<Option<(PathBuf, usize)>, Box<dyn std::error::Error>> {
    let target_path = override_path
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(out_dir).join("high_scoring.faa"));
    let mut writer: Option<BufWriter<File>> = None;
    let mut count = 0usize;
    for m in metrics {
        let Some((_, classif)) = scores_map.get(&m.gene_id) else {
            continue;
        };
        if classification_base(classif) != "High" {
            continue;
        }
        let Some((_im, seq)) = intrinsic.get(&m.gene_id) else {
            continue;
        };
        if writer.is_none() {
            if let Some(parent) = target_path.parent() {
                if !parent.as_os_str().is_empty() {
                    fs::create_dir_all(parent)?;
                }
            }
            writer = Some(BufWriter::new(File::create(&target_path)?));
        }
        if let Some(buf) = writer.as_mut() {
            write_fasta_record(buf, &m.gene_id, seq)?;
            count += 1;
        }
    }
    if let Some(mut buf) = writer {
        buf.flush()?;
        if count == 0 {
            drop(buf);
            fs::remove_file(&target_path).ok();
            Ok(None)
        } else {
            Ok(Some((target_path, count)))
        }
    } else {
        Ok(None)
    }
}

fn write_fasta_record<W: Write>(writer: &mut W, gene_id: &str, seq: &[u8]) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(gene_id.as_bytes())?;
    writer.write_all(b"\n")?;
    for chunk in seq.chunks(60) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
