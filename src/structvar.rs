use crate::diamond::DiamondHitRow;

#[derive(Debug, Clone, Default)]
pub struct StructVar {
    pub fusion_possible: bool,
    pub split_possible: bool,
    pub duplication_possible: bool,
    pub spans: Vec<(usize, usize)>,
    pub classification: String,
}

/// Heuristic structural variation analysis using DIAMOND top hits.
///
/// - fusion_possible: two or more high-quality hits align to disjoint query regions
///   separated by a gap >= 30 aa.
/// - split_possible: large |qcov - scov| as an indicator of split/partial mapping.
/// - duplication_possible: overlapping query spans among multiple hits (same or different
///   subjects) suggesting internal duplication.
pub fn analyze(hits: &[DiamondHitRow]) -> StructVar {
    let mut sv = StructVar::default();
    if hits.is_empty() { return sv; }
    let mut spans: Vec<(usize, usize)> = hits
        .iter()
        .map(|h| {
            let a = h.qstart.min(h.qend);
            let b = h.qstart.max(h.qend);
            (a, b)
        })
        .collect();
    spans.sort_by_key(|t| t.0);
    sv.spans = spans.clone();
    // fusion: look for two non-overlapping spans with a sizable gap
    for i in 0..spans.len().saturating_sub(1) {
        let (a1, b1) = spans[i];
        let (a2, b2) = spans[i+1];
        if a2 > b1 {
            let gap = a2.saturating_sub(b1);
            if gap >= 30 { sv.fusion_possible = true; break; }
        }
    }
    // split: coverage delta from the best hit if available
    if let Some(best) = hits.first() {
        let qcov = best.qcov;
        let scov = best.scov;
        if (qcov - scov).abs() >= 0.30 { sv.split_possible = true; }
    }
    // duplication: any substantial overlap (>50 aa) between two spans
    'dup: for i in 0..spans.len().saturating_sub(1) {
        for j in i+1..spans.len() {
            let (a1, b1) = spans[i];
            let (a2, b2) = spans[j];
            let overlap = (b1.min(b2)).saturating_sub(a1.max(a2));
            if overlap >= 50 {
                sv.duplication_possible = true;
                break 'dup;
            }
        }
    }
    sv.classification = if sv.fusion_possible { "FusionPossible".into() }
        else if sv.split_possible { "SplitPossible".into() }
        else if sv.duplication_possible { "InternalDuplicationPossible".into() }
        else { "None".into() };
    sv
}

