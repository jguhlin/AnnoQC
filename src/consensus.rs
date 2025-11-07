use crate::diamond::DiamondHitRow;

#[derive(Debug, Clone)]
pub struct ConsensusConfig {
    pub min_hits: usize,
    pub max_panel: usize,
    pub filt_qcov: f64,
    pub filt_scov: f64,
    pub filt_evalue: f64,
    pub filt_pident: f64,
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            min_hits: 5,
            max_panel: 10,
            filt_qcov: 0.70,
            filt_scov: 0.50,
            filt_evalue: 1e-5,
            filt_pident: 25.0,
        }
    }
}

/// Select a consensus panel for a single query from its ranked DIAMOND hits,
/// applying filters with phased relaxation and a simple redundancy heuristic.
pub fn select_panel(
    hits: &[DiamondHitRow],
    cfg: &ConsensusConfig,
) -> Vec<String> {
    if hits.is_empty() { return Vec::new(); }

    // Phased filters
    let phases = [
        (cfg.filt_qcov, cfg.filt_scov, cfg.filt_evalue, cfg.filt_pident),
        (0.50, 0.30, 1e-3, 20.0),
        (0.30, 0.10, 1e-2, 15.0),
    ];

    for (qcov_min, scov_min, eval_max, pid_min) in phases {
        let mut selected: Vec<String> = Vec::new();
        let mut seen_subject_roots: std::collections::HashSet<String> = Default::default();
        for h in hits.iter() {
            if h.qcov < qcov_min || h.scov < scov_min { continue; }
            // parse evalue string to f64
            let ev: f64 = match h.evalue.parse::<f64>() { Ok(x) => x, Err(_) => 1.0 };
            if ev > eval_max { continue; }
            if h.pident < pid_min { continue; }
            // Redundancy heuristic: keep 1 per subject root (strip after first space and version)
            let root = subject_root(&h.sseqid);
            if seen_subject_roots.contains(&root) { continue; }
            seen_subject_roots.insert(root);
            selected.push(h.sseqid.clone());
            if selected.len() >= cfg.max_panel { break; }
        }
        if selected.len() >= cfg.min_hits { return selected; }
        // else relax and retry
    }
    // If still not enough, return whatever best we have up to max_panel using last relaxed phase
    let mut selected: Vec<String> = Vec::new();
    let mut seen_subject_roots: std::collections::HashSet<String> = Default::default();
    for h in hits.iter() {
        let root = subject_root(&h.sseqid);
        if seen_subject_roots.contains(&root) { continue; }
        seen_subject_roots.insert(root);
        selected.push(h.sseqid.clone());
        if selected.len() >= cfg.max_panel { break; }
    }
    selected
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
            qseqid: q.clone(), sseqid: format!("S{}", i), bitscore: 100.0 - i as f64, evalue: "1e-20".to_string(),
            length: 100, qcov: 0.8, scov: 0.7, pident: p,
            qstart: 1, qend: 100, sstart: 1, send: 100,
            qlen: 100, slen: 100
        };
        let hits = vec![mk(1, 60.0), mk(2, 55.0), mk(3, 40.0), mk(4, 30.0), mk(5, 28.0), mk(6, 26.0), mk(7, 24.0)];
        let cfg = ConsensusConfig::default();
        let sel = select_panel(&hits, &cfg);
        assert!(sel.len() >= 5 && sel.len() <= cfg.max_panel);
    }
}
