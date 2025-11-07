use crate::diamond::DiamondHitRow;
use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct StructVar {
    pub fusion_possible: bool,
    pub split_possible: bool,
    pub duplication_possible: bool,
    pub spans: Vec<(usize, usize)>,
    pub classification: String,
    // Diagnostics
    pub fusion_gap: Option<usize>,
    pub fusion_left_len: Option<usize>,
    pub fusion_right_len: Option<usize>,
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
    // Heuristic thresholds
    const MIN_STRONG_LEN: usize = 50;      // aa
    const MIN_STRONG_FRAC: f64 = 0.20;     // of query
    const FUSION_MIN_GAP: usize = 50;      // aa gap between disjoint spans from different subjects
    const DUP_MAX_GAP: usize = 20;         // small gap between repeated spans within same subject
    const SPLIT_DELTA: f64 = 0.30;         // |qcov - scov| threshold

    // Group HSPs by subject (canonical sseqid string as-is; upstream may canonicalize as needed)
    let mut by_subject: HashMap<&str, Vec<(usize, usize)>> = HashMap::new();
    let mut by_subject_strong: HashMap<&str, Vec<(usize, usize)>> = HashMap::new();
    for h in hits {
        let a = h.qstart.min(h.qend);
        let b = h.qstart.max(h.qend);
        by_subject.entry(&h.sseqid).or_default().push((a, b));
        // strong HSP filter
        let hsp_len = b.saturating_sub(a) + 1;
        let frac = if h.qlen > 0 { (hsp_len as f64) / (h.qlen as f64) } else { 0.0 };
        if hsp_len >= MIN_STRONG_LEN && frac >= MIN_STRONG_FRAC { by_subject_strong.entry(&h.sseqid).or_default().push((a,b)); }
    }
    // Merge overlapping HSPs per subject to get coarse coverage spans
    fn merge(mut v: Vec<(usize,usize)>) -> Vec<(usize,usize)> {
        if v.is_empty() { return v; }
        v.sort_by_key(|t| t.0);
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for (a, b) in v {
            if let Some(last) = merged.last_mut() {
                if a <= last.1 + 5 { last.1 = last.1.max(b); } else { merged.push((a, b)); }
            } else { merged.push((a, b)); }
        }
        merged
    }
    let mut subj_spans: Vec<(String, Vec<(usize, usize)>)> = Vec::new();
    let mut subj_spans_strong: Vec<(String, Vec<(usize, usize)>)> = Vec::new();
    for (sid, v) in by_subject.into_iter() {
        let mut v2 = merge(v);
        v2.sort_by_key(|t| t.0);
        subj_spans.push((sid.to_string(), v2));
    }
    for (sid, v) in by_subject_strong.into_iter() {
        subj_spans_strong.push((sid.to_string(), merge(v)));
    }
    // Flatten all spans for summary and diagnostics
    let mut spans: Vec<(usize, usize)> = subj_spans.iter().flat_map(|(_, v)| v.clone()).collect();
    spans.sort_by_key(|t| t.0);
    sv.spans = spans.clone();
    // fusion: two different subjects whose (merged, strong) spans are disjoint, separated by large gap
    'fusion: for i in 0..subj_spans_strong.len().saturating_sub(1) {
        for j in i+1..subj_spans_strong.len() {
            let (_, ref a) = subj_spans[i];
            let (_, ref b) = subj_spans[j];
            // use coarse first and last positions
            let (a_min, a_max) = (a.first().map(|x| x.0).unwrap_or(0), a.last().map(|x| x.1).unwrap_or(0));
            let (b_min, b_max) = (b.first().map(|x| x.0).unwrap_or(0), b.last().map(|x| x.1).unwrap_or(0));
            if a_max > 0 && b_max > 0 {
                let (left_max, right_min) = if a_min <= b_min { (a_max, b_min) } else { (b_max, a_min) };
                if right_min > left_max {
                    let gap = right_min - left_max;
                    if gap >= FUSION_MIN_GAP {
                        sv.fusion_possible = true;
                        let left_len = (left_max - if a_min <= b_min { a_min } else { b_min }) + 1;
                        let right_len = (if a_min <= b_min { b_max } else { a_max }) - right_min + 1;
                        sv.fusion_gap = Some(gap);
                        sv.fusion_left_len = Some(left_len);
                        sv.fusion_right_len = Some(right_len);
                        break 'fusion;
                    }
                }
            }
        }
    }
    // split: coverage delta from the best hit if available
    if let Some(best) = hits.first() {
        let qcov = best.qcov;
        let scov = best.scov;
        if (qcov - scov).abs() >= SPLIT_DELTA { sv.split_possible = true; }
    }
    // duplication: within a subject, two strong disjoint spans with small gap suggest internal duplication
    'dup: for (_, v) in &subj_spans_strong {
        for i in 0..v.len().saturating_sub(1) {
            let (a1, b1) = v[i];
            let (a2, b2) = v[i+1];
            if a2 > b1 {
                let gap = a2 - b1;
                if gap <= DUP_MAX_GAP && (b1 - a1 + 1) >= MIN_STRONG_LEN && (b2 - a2 + 1) >= MIN_STRONG_LEN {
                    sv.duplication_possible = true; break 'dup;
                }
            }
        }
    }
    sv.classification = if sv.fusion_possible { "FusionPossible".into() }
        else if sv.split_possible { "SplitPossible".into() }
        else if sv.duplication_possible { "InternalDuplicationPossible".into() }
        else { "None".into() };
    sv
}
