use crate::diamond::DiamondHitRow;
use std::cmp::Ordering;
use std::collections::HashSet;

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
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_hits: 5,
            max_panel: 20,
            filt_qcov: 0.70,
            filt_scov: 0.50,
            filt_evalue: 1e-5,
            filt_pident: 25.0,
            redundancy_pident: 90.0,
            max_high_identity: 3,
            len_ratio_tolerance: 0.30,
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
}

#[derive(Debug, Clone, Default)]
pub struct PanelSelection {
    pub ids: Vec<String>,
    pub stats: PanelStats,
}

struct PhaseFilter {
    qcov: f64,
    scov: f64,
    evalue: f64,
    pident: f64,
    len_min: f64,
    len_max: f64,
    idx: usize,
}

/// Select a consensus panel for a single query from its ranked DIAMOND hits,
/// applying filters with phased relaxation, redundancy trimming, and length heuristics.
pub fn select_panel(hits: &[DiamondHitRow], cfg: &ConsensusConfig) -> PanelSelection {
    if hits.is_empty() {
        return PanelSelection {
            ids: Vec::new(),
            stats: PanelStats {
                total_hits: 0,
                ..Default::default()
            },
        };
    }
    debug_assert!(
        hits.windows(2).all(|w| w[0].qseqid == w[1].qseqid),
        "grouped hits contain mismatched query ids"
    );

    let len_tol = cfg.len_ratio_tolerance.max(0.05);
    let tight_min = (1.0 - len_tol).max(0.1);
    let tight_max = 1.0 + len_tol;
    let phases = [
        PhaseFilter {
            qcov: cfg.filt_qcov,
            scov: cfg.filt_scov,
            evalue: cfg.filt_evalue,
            pident: cfg.filt_pident,
            len_min: tight_min,
            len_max: tight_max,
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
    for phase in phases.iter() {
        let selection = run_phase(hits, cfg, phase);
        if selection.ids.len() >= cfg.min_hits {
            return selection;
        }
        if selection.ids.len() > best.ids.len() {
            best = selection;
        }
    }

    if best.ids.is_empty() {
        // Fallback: keep first unique subjects up to max_panel without filtering
        let mut seen: HashSet<String> = HashSet::new();
        let mut ids: Vec<String> = Vec::new();
        for h in hits.iter() {
            let root = subject_root(&h.sseqid);
            if seen.insert(root) {
                ids.push(h.sseqid.clone());
            }
            if ids.len() >= cfg.max_panel {
                break;
            }
        }
        let fallback_stats = PanelStats {
            total_hits: hits.len(),
            selected: ids.len(),
            phase: 0,
            ..Default::default()
        };
        PanelSelection {
            ids,
            stats: fallback_stats,
        }
    } else {
        best
    }
}

fn run_phase(hits: &[DiamondHitRow], cfg: &ConsensusConfig, phase: &PhaseFilter) -> PanelSelection {
    let mut stats = PanelStats {
        total_hits: hits.len(),
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

    for h in hits.iter() {
        if h.qcov < phase.qcov || h.scov < phase.scov {
            continue;
        }
        let ev: f64 = h.evalue.parse::<f64>().unwrap_or(1.0);
        if ev > phase.evalue {
            continue;
        }
        if h.pident < phase.pident {
            continue;
        }
        let len_ratio = if h.qlen > 0 {
            h.slen as f64 / h.qlen as f64
        } else {
            0.0
        };
        if len_ratio < phase.len_min || len_ratio > phase.len_max {
            continue;
        }
        stats.filtered_hits += 1;
        let root = subject_root(&h.sseqid);
        if !seen_roots.insert(root) {
            continue;
        }
        if h.pident >= cfg.redundancy_pident {
            if high_identity_selected >= max_high_identity {
                stats.high_identity_dropped += 1;
                continue;
            }
            high_identity_selected += 1;
        }
        selected.push(h.sseqid.clone());
        len_ratios.push(len_ratio);
        pidents.push(h.pident);
        if selected.len() >= cfg.max_panel {
            break;
        }
    }

    stats.selected = selected.len();
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
        ids: selected,
        stats,
    }
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
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
}
