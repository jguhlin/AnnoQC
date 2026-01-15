use crate::diamond::{DiamondHitRow, HitSource};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ConsensusConfig {
    pub min_hits: usize,
    pub max_panel: usize,
    pub filt_qcov: f64,
    pub filt_scov: f64,
    pub filt_evalue: f64,
    pub filt_pident: f64,
    pub redundancy_pident: f64,
    pub max_high_identity: usize,
    pub len_ratio_tolerance: f64,
    pub redundancy_trigger: usize,
    pub diversity_backfill: bool,
    pub clusters: Option<Arc<HashMap<String, Vec<String>>>>,
    pub backfill_enabled: bool,
    pub backfill_min_primary_hits: usize,
    pub backfill_max_added: usize,
    pub refprot_proteome_cap: usize,
    pub diversity_rank_index: usize,
    pub diversity_rank_cap: usize,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_hits: 5,
            max_panel: 40,
            filt_qcov: 0.70,
            filt_scov: 0.50,
            filt_evalue: 1e-5,
            filt_pident: 25.0,
            redundancy_pident: 90.0,
            max_high_identity: 3,
            len_ratio_tolerance: 0.30,
            redundancy_trigger: 20,
            diversity_backfill: true,
            clusters: None,
            backfill_enabled: true,
            backfill_min_primary_hits: 1,
            backfill_max_added: 5,
            refprot_proteome_cap: 3,
            diversity_rank_index: 3,
            diversity_rank_cap: 3,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct PanelStats {
    pub total_hits: usize,
    pub filtered_hits: usize,
    pub selected: usize,
    pub phase: usize,
    pub len_ratio_window_min: f64,
    pub len_ratio_window_max: f64,
    pub median_len_ratio: f64,
    pub len_ratio_min: f64,
    pub len_ratio_max: f64,
    pub median_pident: f64,
    pub high_identity_dropped: usize,
    pub backfill_from_clusters: usize,
    pub diversity_rank_index: usize,
    pub diversity_cap: usize,
    pub diversity_keys_used: usize,
    pub diversity_skipped: usize,
}

#[derive(Debug, Clone, Default)]
pub struct PanelSelection {
    pub ids: Vec<String>,
    pub stats: PanelStats,
}

#[derive(Debug, Clone, Default)]
pub struct LenStats {
    pub median: f64,
    pub mad: f64,
    pub expected_len: usize,
    pub qlen: usize,
}

#[derive(Debug, Clone)]
pub struct PanelResult {
    pub selection: PanelSelection,
    pub len_stats: LenStats,
    pub aggregated_hits: Vec<AggregatedHit>,
}

#[derive(Debug, Clone)]
pub struct AggregatedHit {
    pub sseqid: String,
    pub qcov: f64,
    pub scov: f64,
    pub bitscore: f64,
    pub evalue: String,
    pub qlen: usize,
    #[allow(dead_code)]
    pub slen: usize,
    pub len_ratio: f64,
    pub pident: f64,
    pub hit_count: usize,
    pub source: HitSource,
    pub quality: f64,
    pub taxid: Option<u32>,
    pub lineage: Vec<String>,
}

impl AggregatedHit {
    pub fn diversity_key(&self, rank_index: usize) -> Option<String> {
        if !self.lineage.is_empty() {
            if rank_index < self.lineage.len() {
                let key = self.lineage[rank_index].clone();
                if !key.is_empty() {
                    return Some(key);
                }
            } else if let Some(last) = self.lineage.last() {
                if !last.is_empty() {
                    return Some(last.clone());
                }
            }
        }
        if let Some(tid) = self.taxid {
            return Some(format!("taxid:{}", tid));
        }
        if let HitSource::RefProt(pid) = &self.source {
            return Some(format!("refprot:{}", pid));
        }
        None
    }
}

#[derive(Clone)]
struct PhaseFilter {
    qcov: f64,
    #[allow(dead_code)]
    scov: f64,
    evalue: f64,
    pident: f64,
    len_min: f64,
    len_max: f64,
    idx: usize,
}

struct AggregatedHitBuilder {
    q_spans: Vec<(usize, usize)>,
    s_spans: Vec<(usize, usize)>,
    qlen_max: usize,
    slen_max: usize,
    hit_count: usize,
    best_bitscore: f64,
    best_pident: f64,
    best_evalue: Option<String>,
    best_source: HitSource,
    best_taxid: Option<u32>,
    best_lineage: Vec<String>,
}

impl AggregatedHitBuilder {
    fn new() -> Self {
        Self {
            q_spans: Vec::new(),
            s_spans: Vec::new(),
            qlen_max: 0,
            slen_max: 0,
            hit_count: 0,
            best_bitscore: f64::NEG_INFINITY,
            best_pident: 0.0,
            best_evalue: None,
            best_source: HitSource::SwissProt,
            best_taxid: None,
            best_lineage: Vec::new(),
        }
    }

    fn push(&mut self, row: &DiamondHitRow) {
        if row.bitscore > self.best_bitscore {
            self.best_bitscore = row.bitscore;
            self.best_evalue = Some(row.evalue.clone());
            self.best_source = row.source.clone();
        }
        if row.pident > self.best_pident {
            self.best_pident = row.pident;
        }
        if self.best_taxid.is_none() {
            self.best_taxid = row.staxid;
        }
        if self.best_lineage.is_empty() && !row.lineage.is_empty() {
            self.best_lineage = row.lineage.clone();
        }
        if row.qlen > 0 {
            self.qlen_max = self.qlen_max.max(row.qlen);
        }
        if row.slen > 0 {
            self.slen_max = self.slen_max.max(row.slen);
        }
        if row.qend >= row.qstart {
            let start = row.qstart.min(row.qend);
            let end = row.qstart.max(row.qend);
            self.q_spans.push((start, end));
        }
        if row.send >= row.sstart {
            let start = row.sstart.min(row.send);
            let end = row.sstart.max(row.send);
            self.s_spans.push((start, end));
        }
        self.hit_count += 1;
    }

    fn finalize(self, sseqid: String) -> AggregatedHit {
        let mut q_spans = self.q_spans;
        let mut s_spans = self.s_spans;
        let qlen = if self.qlen_max > 0 {
            self.qlen_max
        } else {
            q_spans.iter().map(|span| span.1).max().unwrap_or(0)
        };
        let slen = if self.slen_max > 0 {
            self.slen_max
        } else {
            s_spans.iter().map(|span| span.1).max().unwrap_or(0)
        };
        let qcov_len = union_span_len(&mut q_spans);
        let scov_len = union_span_len(&mut s_spans);
        let qcov = if qlen > 0 {
            (qcov_len as f64 / qlen as f64).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let scov = if slen > 0 {
            (scov_len as f64 / slen as f64).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let len_ratio = if qlen > 0 {
            slen as f64 / (qlen as f64)
        } else {
            0.0
        };
        let base_quality = if len_ratio > 0.0 {
            self.best_bitscore * qcov * (1.0 / (1.0 + (len_ratio - 1.0).abs())).max(0.01)
        } else {
            self.best_bitscore * qcov
        };
        let scov_factor = match self.best_source {
            HitSource::RefProt(_) => (0.5 + 0.5 * scov).clamp(0.1, 1.0),
            _ => 1.0,
        };
        let quality = base_quality * scov_factor;
        AggregatedHit {
            sseqid,
            qcov,
            scov,
            bitscore: self.best_bitscore.max(0.0),
            evalue: self.best_evalue.unwrap_or_else(|| "1".to_string()),
            qlen,
            slen,
            len_ratio,
            pident: self.best_pident,
            hit_count: self.hit_count,
            source: self.best_source,
            quality,
            taxid: self.best_taxid,
            lineage: self.best_lineage,
        }
    }
}

impl Default for AggregatedHitBuilder {
    fn default() -> Self {
        Self::new()
    }
}

fn union_span_len(spans: &mut [(usize, usize)]) -> usize {
    if spans.is_empty() {
        return 0;
    }
    spans.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut total = 0;
    let mut current = spans[0];
    for span in spans.iter().skip(1) {
        if span.0 <= current.1 {
            current.1 = current.1.max(span.1);
        } else {
            total += current.1 - current.0 + 1;
            current = *span;
        }
    }
    total + (current.1 - current.0 + 1)
}

fn compute_len_stats(hits: &[AggregatedHit]) -> LenStats {
    if hits.is_empty() {
        return LenStats::default();
    }
    let qlen = hits.iter().map(|h| h.qlen).max().unwrap_or(0);
    let mut ratios: Vec<f64> = hits
        .iter()
        .map(|h| h.len_ratio)
        .filter(|&lr| lr.is_finite() && lr > 0.0)
        .collect();
    if ratios.is_empty() {
        return LenStats {
            median: 1.0,
            mad: 0.0,
            expected_len: qlen,
            qlen,
        };
    }
    let median_val = median(&mut ratios);
    let mut deviations: Vec<f64> = ratios.iter().map(|lr| (lr - median_val).abs()).collect();
    let mad_val = if deviations.is_empty() {
        0.0
    } else {
        median(&mut deviations)
    };
    let expected_len = ((median_val * qlen as f64).round() as usize).max(1);
    LenStats {
        median: median_val,
        mad: mad_val,
        expected_len,
        qlen,
    }
}

/// Expand a phase's length window using median shift and extra slack.
///
/// `len_shift` recenters the window around observed median ratios, `extra` widens
/// it based on dispersion, and `len_tol` sets the baseline tolerance.
fn adjust_phase_window(
    phase: &PhaseFilter,
    len_shift: f64,
    extra: f64,
    len_tol: f64,
) -> PhaseFilter {
    let mut dyn_phase = phase.clone();
    let width = extra + len_tol * 0.5;
    let min = (phase.len_min + len_shift - width).max(0.05);
    let max = (phase.len_max + len_shift + width).min(3.0);
    dyn_phase.len_min = min;
    dyn_phase.len_max = if max > dyn_phase.len_min {
        max
    } else {
        dyn_phase.len_min + 0.05
    };
    dyn_phase
}

fn aggregate_hits_by_subject(hits: &[DiamondHitRow]) -> Vec<AggregatedHit> {
    let mut map: HashMap<String, AggregatedHitBuilder> = HashMap::new();
    for hit in hits.iter() {
        map.entry(hit.sseqid.clone()).or_default().push(hit);
    }
    let mut aggregated: Vec<AggregatedHit> = map
        .into_iter()
        .map(|(sseqid, builder)| builder.finalize(sseqid))
        .collect();
    aggregated.sort_by(|a, b| {
        b.quality
            .partial_cmp(&a.quality)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    aggregated
}

/// Cap RefProt hits per proteome key; `cap == 0` drops all RefProt hits.
fn cap_refprot_hits(hits: Vec<AggregatedHit>, cap: usize) -> Vec<AggregatedHit> {
    if cap == 0 {
        return hits
            .into_iter()
            .filter(|h| !matches!(h.source, HitSource::RefProt(_)))
            .collect();
    }
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut kept: Vec<AggregatedHit> = Vec::new();
    for hit in hits.into_iter() {
        if let HitSource::RefProt(prot) = &hit.source {
            let entry = counts.entry(prot.clone()).or_insert(0);
            if *entry >= cap {
                continue;
            }
            *entry += 1;
        }
        kept.push(hit);
    }
    kept
}

/// Select a consensus panel for a single query from its ranked DIAMOND hits.
///
/// Filters are applied in phases (strict to relaxed) using the config thresholds,
/// with redundancy capping (`redundancy_*`) and optional diversity constraints
/// (`diversity_*`). Length ratio windows are adjusted from the median/MAD of
/// aggregated hits and `len_ratio_tolerance`. When enabled, backfill can
/// supplement scarce panels from configured clusters and length-ratio bins,
/// bounded by `backfill_min_primary_hits`, `backfill_max_added`, and `max_panel`.
#[allow(dead_code)]
pub fn select_panel(hits: &[DiamondHitRow], cfg: &ConsensusConfig) -> PanelSelection {
    select_panel_with_result(hits, cfg).selection
}

/// Full panel selection with aggregated hits and length stats returned.
///
/// Filters are applied in phases (strict to relaxed) using the config thresholds,
/// with redundancy capping (`redundancy_*`) and optional diversity constraints
/// (`diversity_*`). Length ratio windows are adjusted from the median/MAD of
/// aggregated hits and `len_ratio_tolerance`. When enabled, backfill can
/// supplement scarce panels from configured clusters and length-ratio bins,
/// bounded by `backfill_min_primary_hits`, `backfill_max_added`, and `max_panel`.
pub fn select_panel_with_result(hits: &[DiamondHitRow], cfg: &ConsensusConfig) -> PanelResult {
    let total_hits = hits.len();
    if hits.is_empty() {
        return PanelResult {
            selection: PanelSelection {
                ids: Vec::new(),
                stats: PanelStats {
                    total_hits,
                    ..Default::default()
                },
            },
            len_stats: LenStats::default(),
            aggregated_hits: Vec::new(),
        };
    }
    debug_assert!(
        hits.windows(2).all(|w| w[0].qseqid == w[1].qseqid),
        "grouped hits contain mismatched query ids"
    );

    let aggregated_hits = aggregate_hits_by_subject(hits);
    let aggregated_hits = cap_refprot_hits(aggregated_hits, cfg.refprot_proteome_cap);
    let len_stats = compute_len_stats(&aggregated_hits);

    if aggregated_hits.is_empty() {
        return PanelResult {
            selection: PanelSelection {
                ids: Vec::new(),
                stats: PanelStats {
                    total_hits,
                    ..Default::default()
                },
            },
            len_stats,
            aggregated_hits,
        };
    }

    let len_tol = if cfg.len_ratio_tolerance.is_finite() {
        cfg.len_ratio_tolerance
    } else {
        ConsensusConfig::default().len_ratio_tolerance
    }
    .clamp(0.05, 2.0);
    let len_shift = (1.0 - len_stats.median).clamp(-0.4, 0.4);
    let dyn_extra = len_stats.mad.max(0.05);
    let base_phases = [
        PhaseFilter {
            qcov: cfg.filt_qcov,
            scov: cfg.filt_scov,
            evalue: cfg.filt_evalue,
            pident: cfg.filt_pident,
            len_min: (1.0 - len_tol).max(0.1),
            len_max: 1.0 + len_tol,
            idx: 1,
        },
        PhaseFilter {
            qcov: 0.50,
            scov: 0.30,
            evalue: 1e-3,
            pident: 20.0,
            len_min: 0.50,
            len_max: 1.50,
            idx: 2,
        },
        PhaseFilter {
            qcov: 0.30,
            scov: 0.10,
            evalue: 1e-2,
            pident: 15.0,
            len_min: 0.40,
            len_max: 2.00,
            idx: 3,
        },
    ];

    let mut best = PanelSelection::default();
    for phase in base_phases.iter() {
        let dyn_phase = adjust_phase_window(phase, len_shift, dyn_extra, len_tol);
        let selection = run_phase(&aggregated_hits, cfg, &dyn_phase, total_hits);
        if selection.ids.len() >= cfg.min_hits {
            return PanelResult {
                selection,
                len_stats: len_stats.clone(),
                aggregated_hits,
            };
        }
        if selection.ids.len() > best.ids.len() {
            best = selection;
        }
    }

    PanelResult {
        selection: best,
        len_stats,
        aggregated_hits,
    }
}

fn run_phase(
    hits: &[AggregatedHit],
    cfg: &ConsensusConfig,
    phase: &PhaseFilter,
    total_hits: usize,
) -> PanelSelection {
    let mut stats = PanelStats {
        total_hits,
        phase: phase.idx,
        len_ratio_window_min: phase.len_min,
        len_ratio_window_max: phase.len_max,
        ..Default::default()
    };
    let mut selected: Vec<String> = Vec::new();
    let mut len_ratios: Vec<f64> = Vec::new();
    let mut pidents: Vec<f64> = Vec::new();
    let mut seen_roots: HashSet<String> = HashSet::new();
    let mut high_identity_selected = 0usize;
    let max_high_identity = cfg.max_high_identity.max(1);
    let enforce_redundancy = stats.filtered_hits >= cfg.redundancy_trigger;
    let enforce_diversity = cfg.diversity_rank_cap > 0;
    let mut diversity_counts: HashMap<String, usize> = HashMap::new();
    let mut diversity_spill: Vec<usize> = Vec::new();

    for (idx, h) in hits.iter().enumerate() {
        if h.qcov < phase.qcov {
            continue;
        }
        let ev: f64 = h
            .evalue
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v >= 0.0)
            .unwrap_or(f64::INFINITY);
        if ev > phase.evalue {
            continue;
        }
        if h.pident < phase.pident {
            continue;
        }
        let len_ratio = h.len_ratio;
        if len_ratio < phase.len_min || len_ratio > phase.len_max {
            continue;
        }
        stats.filtered_hits += 1;
        let root = subject_root(&h.sseqid);
        if seen_roots.contains(&root) {
            continue;
        }
        if enforce_redundancy && h.pident >= cfg.redundancy_pident {
            if high_identity_selected >= max_high_identity {
                stats.high_identity_dropped += 1;
                continue;
            }
            high_identity_selected += 1;
        }
        let mut diversity_key: Option<String> = None;
        if enforce_diversity {
            diversity_key = h.diversity_key(cfg.diversity_rank_index);
            if let Some(key) = &diversity_key {
                let count = diversity_counts.get(key).copied().unwrap_or(0);
                if count >= cfg.diversity_rank_cap {
                    stats.diversity_skipped += 1;
                    diversity_spill.push(idx);
                    continue;
                }
            }
        }
        seen_roots.insert(root);
        selected.push(h.sseqid.clone());
        if len_ratio.is_finite() {
            len_ratios.push(len_ratio);
        }
        if h.pident.is_finite() {
            pidents.push(h.pident);
        }
        if let Some(key) = diversity_key {
            *diversity_counts.entry(key).or_insert(0) += 1;
        }
        if selected.len() >= cfg.max_panel {
            break;
        }
    }

    if selected.len() < cfg.min_hits && !diversity_spill.is_empty() {
        for idx in diversity_spill {
            if selected.len() >= cfg.min_hits {
                break;
            }
            if let Some(h) = hits.get(idx) {
                if selected.contains(&h.sseqid) {
                    continue;
                }
                let root = subject_root(&h.sseqid);
                if seen_roots.contains(&root) {
                    continue;
                }
                seen_roots.insert(root);
                let len_ratio = h.len_ratio;
                selected.push(h.sseqid.clone());
                len_ratios.push(len_ratio);
                pidents.push(h.pident);
            }
        }
    }

    stats.selected = selected.len();
    stats.diversity_rank_index = cfg.diversity_rank_index;
    stats.diversity_cap = cfg.diversity_rank_cap;
    stats.diversity_keys_used = diversity_counts.len();
    if !len_ratios.is_empty() {
        stats.len_ratio_min = len_ratios
            .iter()
            .fold(f64::INFINITY, |acc, &x| if x < acc { x } else { acc });
        stats.len_ratio_max =
            len_ratios
                .iter()
                .fold(f64::NEG_INFINITY, |acc, &x| if x > acc { x } else { acc });
        let mut tmp = len_ratios.clone();
        stats.median_len_ratio = median(&mut tmp);
    }
    if !pidents.is_empty() {
        let mut tmp = pidents.clone();
        stats.median_pident = median(&mut tmp);
    }

    PanelSelection {
        ids: finalize_with_backfill(hits, cfg, selected, &mut stats),
        stats,
    }
}

fn finalize_with_backfill(
    hits: &[AggregatedHit],
    cfg: &ConsensusConfig,
    mut selected: Vec<String>,
    stats: &mut PanelStats,
) -> Vec<String> {
    // Backfill from clusters if scarce
    if cfg.backfill_enabled
        && selected.len() < cfg.min_hits
        && selected.len() >= cfg.backfill_min_primary_hits
    {
        if let Some(clmap) = &cfg.clusters {
            let mut added = 0usize;
            let mut seen: HashSet<String> = selected.iter().cloned().collect();
            for h in hits.iter().take(2) {
                let root = subject_root(&h.sseqid);
                if let Some(members) = clmap.get(&root).or_else(|| clmap.get(&h.sseqid)) {
                    for m in members {
                        if seen.insert(m.clone()) {
                            selected.push(m.clone());
                            added += 1;
                            if added >= cfg.backfill_max_added {
                                break;
                            }
                            if selected.len() >= cfg.min_hits || selected.len() >= cfg.max_panel {
                                break;
                            }
                        }
                    }
                }
                if selected.len() >= cfg.min_hits || selected.len() >= cfg.max_panel {
                    break;
                }
            }
            stats.backfill_from_clusters = added;
        }
    }

    // Diversity backfill by length ratio bins
    if cfg.diversity_backfill && selected.len() < cfg.max_panel {
        let mut bins_seen: HashSet<i32> = HashSet::new();
        let mut seen_ids: HashSet<String> = selected.iter().cloned().collect();
        // mark bins of already selected
        for h in hits {
            if seen_ids.contains(&h.sseqid) {
                let lr = h.len_ratio;
                bins_seen.insert((lr * 10.0).floor() as i32);
            }
        }
        // prefer hits in unseen bins
        for h in hits {
            if selected.len() >= cfg.max_panel {
                break;
            }
            if seen_ids.contains(&h.sseqid) {
                continue;
            }
            let lr = h.len_ratio;
            let bin = (lr * 10.0).floor() as i32;
            if !bins_seen.contains(&bin) {
                selected.push(h.sseqid.clone());
                seen_ids.insert(h.sseqid.clone());
                bins_seen.insert(bin);
            }
        }
        // fill remaining with best leftover
        for h in hits {
            if selected.len() >= cfg.max_panel {
                break;
            }
            if !seen_ids.contains(&h.sseqid) {
                selected.push(h.sseqid.clone());
                seen_ids.insert(h.sseqid.clone());
            }
        }
    }

    selected.truncate(cfg.max_panel);
    selected
}

fn median(values: &mut Vec<f64>) -> f64 {
    values.retain(|v| v.is_finite());
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

fn subject_root(sseqid: &str) -> String {
    // Try to canonicalize: take up to first whitespace, then strip trailing version after '.'
    let head = sseqid.split_whitespace().next().unwrap_or("");
    let mut parts = head.split('.');
    parts.next().unwrap_or("").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(gene: &str, subj: &str, len_ratio: f64, source: HitSource) -> DiamondHitRow {
        let qlen = 1000usize;
        let slen = (len_ratio * qlen as f64).round().max(1.0) as usize;
        DiamondHitRow {
            qseqid: gene.to_string(),
            sseqid: subj.to_string(),
            bitscore: 200.0,
            evalue: "1e-30".into(),
            length: 150,
            qcov: 0.85,
            scov: 0.70,
            pident: 45.0,
            qstart: 1,
            qend: 850,
            sstart: 10,
            send: 860,
            qlen,
            slen,
            source,
            staxid: Some(1),
            lineage: vec!["root".into(), "Test".into()],
        }
    }

    #[test]
    fn basic_selection() {
        let q = "Q".to_string();
        let mk = |i: usize, p: f64| DiamondHitRow {
            qseqid: q.clone(),
            sseqid: format!("S{}", i),
            bitscore: 100.0 - i as f64,
            evalue: "1e-20".to_string(),
            length: 100,
            qcov: 0.8,
            scov: 0.7,
            pident: p,
            qstart: 1,
            qend: 100,
            sstart: 1,
            send: 100,
            qlen: 100,
            slen: 100,
            source: HitSource::SwissProt,
            staxid: Some(1),
            lineage: vec!["root".into(), format!("Class{}", i % 2)],
        };
        let hits = vec![
            mk(1, 96.0),
            mk(2, 92.0),
            mk(3, 70.0),
            mk(4, 60.0),
            mk(5, 55.0),
            mk(6, 50.0),
            mk(7, 45.0),
        ];
        let cfg = ConsensusConfig {
            max_panel: 6,
            ..Default::default()
        };
        let sel = select_panel(&hits, &cfg);
        assert!(sel.ids.len() >= cfg.min_hits.min(cfg.max_panel));
        assert_eq!(sel.stats.selected, sel.ids.len());
        assert!(sel.stats.median_pident > 40.0);
    }

    #[test]
    fn truncation_expands_length_window() {
        let gene = "GeneX";
        // All hits have len_ratio 1.4 (> default 1.3 max), so only dynamic window should rescue them
        let hits: Vec<_> = (0..6)
            .map(|i| make_hit(gene, &format!("Ref{i}"), 1.4, HitSource::SwissProt))
            .collect();
        let sel = select_panel(&hits, &ConsensusConfig::default());
        assert!(
            sel.ids.len() >= 5,
            "dynamic window failed to keep truncated hits"
        );
    }

    #[test]
    fn cap_refprot_limits_per_proteome() {
        let mut hits: Vec<AggregatedHit> = Vec::new();
        for i in 0..5 {
            hits.push(AggregatedHit {
                sseqid: format!("RefA_{i}"),
                qcov: 0.9,
                scov: 0.8,
                bitscore: 300.0 - i as f64,
                evalue: "1e-40".into(),
                qlen: 1000,
                slen: 950,
                len_ratio: 0.95,
                pident: 40.0,
                hit_count: 1,
                source: HitSource::RefProt("P0AAA".into()),
                quality: 10.0,
                taxid: Some(1),
                lineage: vec![],
            });
        }
        hits.push(AggregatedHit {
            sseqid: "Swiss1".into(),
            qcov: 0.9,
            scov: 0.9,
            bitscore: 320.0,
            evalue: "1e-50".into(),
            qlen: 1000,
            slen: 1020,
            len_ratio: 1.02,
            pident: 60.0,
            hit_count: 1,
            source: HitSource::SwissProt,
            quality: 11.0,
            taxid: Some(1),
            lineage: vec![],
        });
        let kept = super::cap_refprot_hits(hits.clone(), 2);
        let refprot_kept = kept
            .iter()
            .filter(|h| matches!(h.source, HitSource::RefProt(_)))
            .count();
        assert_eq!(refprot_kept, 2, "proteome cap should enforce limit");
        assert!(
            kept.iter()
                .any(|h| matches!(h.source, HitSource::SwissProt)),
            "non-refprot hits must remain"
        );
    }

    #[test]
    fn refprot_quality_prefers_high_coverage() {
        let q = "Q".to_string();
        let row_low = DiamondHitRow {
            qseqid: q.clone(),
            sseqid: "R0".into(),
            bitscore: 200.0,
            evalue: "1e-20".into(),
            length: 100,
            qcov: 1.0,
            scov: 0.25,
            pident: 40.0,
            qstart: 1,
            qend: 100,
            sstart: 1,
            send: 100,
            qlen: 100,
            slen: 400,
            source: HitSource::RefProt("PX".into()),
            staxid: Some(1),
            lineage: vec![],
        };
        let mut row_high = row_low.clone();
        row_high.sseqid = "R1".into();
        row_high.slen = 100;
        let hits = vec![row_low, row_high];
        let cfg = ConsensusConfig::default();
        let res = select_panel_with_result(&hits, &cfg);
        assert_eq!(res.aggregated_hits.len(), 2);
        let top = &res.aggregated_hits[0];
        assert_eq!(top.sseqid, "R1");
        assert!(top.quality > res.aggregated_hits[1].quality);
    }
}
