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
    pub fusion_subjects: Option<(String, String)>,
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
    let mut by_subject: HashMap<&str, Vec<(usize, usize, usize, usize)>> = HashMap::new();
    let mut by_subject_strong: HashMap<&str, Vec<(usize, usize, usize, usize)>> = HashMap::new();
    let mut subj_slen: HashMap<&str, usize> = HashMap::new();
    for h in hits {
        let a = h.qstart.min(h.qend);
        let b = h.qstart.max(h.qend);
        by_subject.entry(&h.sseqid).or_default().push((a, b, h.sstart.min(h.send), h.sstart.max(h.send)));
        subj_slen.entry(&h.sseqid).or_insert(h.slen);
        // strong HSP filter
        let hsp_len = b.saturating_sub(a) + 1;
        let frac = if h.qlen > 0 { (hsp_len as f64) / (h.qlen as f64) } else { 0.0 };
        if hsp_len >= MIN_STRONG_LEN && frac >= MIN_STRONG_FRAC {
            by_subject_strong.entry(&h.sseqid).or_default().push((a,b, h.sstart.min(h.send), h.sstart.max(h.send)));
        }
    }
    // Merge overlapping HSPs per subject to get coarse coverage spans
    fn merge(mut v: Vec<(usize,usize,usize,usize)>) -> Vec<(usize,usize,usize,usize)> {
        if v.is_empty() { return v; }
        v.sort_by_key(|t| t.0);
        let mut merged: Vec<(usize, usize, usize, usize)> = Vec::new();
        for (a, b, sa, sb) in v {
            if let Some(last) = merged.last_mut() {
                if a <= last.1 + 5 {
                    // extend both query and subject bounds
                    last.1 = last.1.max(b);
                    last.3 = last.3.max(sb);
                    last.2 = last.2.min(sa);
                } else {
                    merged.push((a, b, sa, sb));
                }
            } else { merged.push((a, b, sa, sb)); }
        }
        merged
    }
    let mut subj_spans: Vec<(String, Vec<(usize, usize, usize, usize)>)> = Vec::new();
    let mut subj_spans_strong: Vec<(String, Vec<(usize, usize, usize, usize)>)> = Vec::new();
    for (sid, v) in by_subject.into_iter() {
        let mut v2 = merge(v);
        v2.sort_by_key(|t| t.0);
        subj_spans.push((sid.to_string(), v2));
    }
    for (sid, v) in by_subject_strong.into_iter() {
        subj_spans_strong.push((sid.to_string(), merge(v)));
    }
    // Flatten all spans for summary and diagnostics
    let mut spans: Vec<(usize, usize)> = subj_spans
        .iter()
        .flat_map(|(_, v)| v.iter().map(|(a,b,_,_)| (*a,*b)))
        .collect();
    spans.sort_by_key(|t| t.0);
    sv.spans = spans.clone();
    // fusion: two different subjects whose (merged, strong) spans are disjoint, separated by large gap
    'fusion: for i in 0..subj_spans_strong.len().saturating_sub(1) {
        for j in i+1..subj_spans_strong.len() {
            let (ref ida, ref a) = subj_spans[i];
            let (ref idb, ref b) = subj_spans[j];
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
                        sv.fusion_subjects = Some((ida.clone(), idb.clone()));
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
            let (a1, b1, _sa1, _sb1) = v[i];
            let (a2, b2, _sa2, _sb2) = v[i+1];
            if a2 > b1 {
                let gap = a2 - b1;
                if gap <= DUP_MAX_GAP && (b1 - a1 + 1) >= MIN_STRONG_LEN && (b2 - a2 + 1) >= MIN_STRONG_LEN {
                    sv.duplication_possible = true; break 'dup;
                }
            }
        }
    }

    // Final classification preference order (fusion > split > duplication)
    sv.classification = if sv.fusion_possible { "FusionPossible".into() }
        else if sv.split_possible { "SplitPossible".into() }
        else if sv.duplication_possible { "InternalDuplicationPossible".into() }
        else { "None".into() };
    sv
}

#[cfg(test)]
mod tests {
    use super::*;
    fn mk(q: &str, s: &str, qs: usize, qe: usize, qlen: usize, sstart: usize, send: usize) -> DiamondHitRow {
        DiamondHitRow { qseqid: q.into(), sseqid: s.into(), bitscore: 200.0, evalue: "1e-20".into(), length: (qe.max(qs)-qe.min(qs)+1), qcov: ((qe.max(qs)-qe.min(qs)+1) as f64)/(qlen as f64), scov: 0.8, pident: 60.0, qstart: qs, qend: qe, sstart, send, qlen, slen: 200 }
    }
    #[test]
    fn fusion_two_subjects_far_apart() {
        // Two strong HSPs on different subjects, well-separated in query
        let hits = vec![ mk("q","A",  50, 300, 1000, 10, 100), mk("q","B", 700, 950, 1000, 20, 180) ];
        let sv = analyze(&hits);
        assert!(sv.fusion_possible);
        assert_eq!(sv.classification, "FusionPossible");
        assert!(sv.fusion_gap.unwrap() >= 50);
    }
    #[test]
    fn duplication_adjacent_within_subject() {
        // Two strong blocks separated by a small gap (<=20 aa) within the same subject
        let hits = vec![ mk("q","A", 100, 320, 1000, 10, 200), mk("q","A", 340, 540, 1000, 220, 420) ];
        let sv = analyze(&hits);
        assert!(sv.duplication_possible);
    }
}
