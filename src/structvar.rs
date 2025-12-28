use crate::diamond::DiamondHitRow;
use serde::Serialize;
use std::collections::{HashMap, HashSet};

type HspBlock = (usize, usize, usize, usize, bool);
type MergedSpan = (usize, usize, usize, usize);
type SubjectHspMap = HashMap<String, Vec<HspBlock>>;

#[derive(Debug, Clone, Default, Serialize)]
pub struct StructVarSubjectSummary {
    pub subject_id: String,
    pub query_coverage: f64,
    pub subject_coverage: f64,
    pub span_count: usize,
    pub orientation: String,
    pub strand_confidence: f64,
    pub order_conflict: bool,
}

#[derive(Debug, Clone)]
pub struct StructVarThresholds {
    pub min_strong_len: usize,
    pub min_strong_frac: f64,
    pub fusion_min_gap: usize,
    pub dup_max_gap: usize,
    pub split_delta: f64,
    pub min_subject_cov: f64,
    pub orient_majority: f64,
}

impl Default for StructVarThresholds {
    fn default() -> Self {
        Self {
            min_strong_len: 50,
            min_strong_frac: 0.20,
            fusion_min_gap: 160,
            dup_max_gap: 20,
            split_delta: 0.30,
            min_subject_cov: 0.27,
            orient_majority: 0.70,
        }
    }
}

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
    pub fusion_cover_fracs: Option<(f64, f64)>,
    pub subjects: Vec<StructVarSubjectSummary>,
    pub warnings: Vec<String>,
}

/// Heuristic structural variation analysis using DIAMOND top hits.
///
/// - fusion_possible: two or more high-quality hits align to disjoint query regions
///   separated by a gap >= 30 aa.
/// - split_possible: large |qcov - scov| as an indicator of split/partial mapping.
/// - duplication_possible: overlapping query spans among multiple hits (same or different
///   subjects) suggesting internal duplication.
pub fn analyze(hits: &[DiamondHitRow], th: &StructVarThresholds) -> StructVar {
    let mut sv = StructVar::default();
    if hits.is_empty() {
        return sv;
    }
    let qlen = hits.first().map(|h| h.qlen.max(1)).unwrap_or(1);
    // Heuristic thresholds (configurable)
    let min_strong_len = th.min_strong_len;
    let min_strong_frac = th.min_strong_frac;
    let fusion_min_gap = th.fusion_min_gap;
    let dup_max_gap = th.dup_max_gap;
    let split_delta = th.split_delta;
    let min_subj_cov = th.min_subject_cov;
    let orient_maj = th.orient_majority;

    // Group HSPs by subject (canonical sseqid string as-is; upstream may canonicalize as needed)
    let mut by_subject: SubjectHspMap = HashMap::new();
    let mut by_subject_strong: SubjectHspMap = HashMap::new();
    let mut subj_slen: HashMap<String, usize> = HashMap::new();
    let mut subj_orient: HashMap<String, (usize, usize)> = HashMap::new(); // (fwd, rev)
    for h in hits {
        let a = h.qstart.min(h.qend);
        let b = h.qstart.max(h.qend);
        let is_rev = h.sstart > h.send;
        let key = h.sseqid.clone();
        by_subject.entry(key.clone()).or_default().push((
            a,
            b,
            h.sstart.min(h.send),
            h.sstart.max(h.send),
            is_rev,
        ));
        subj_slen.entry(key.clone()).or_insert(h.slen);
        let e = subj_orient.entry(key.clone()).or_insert((0, 0));
        if is_rev {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
        // strong HSP filter
        let hsp_len = h.length.max(1);
        let frac = if h.qlen > 0 {
            (hsp_len as f64) / (h.qlen as f64)
        } else {
            0.0
        };
        if hsp_len >= min_strong_len && frac >= min_strong_frac {
            by_subject_strong.entry(key).or_default().push((
                a,
                b,
                h.sstart.min(h.send),
                h.sstart.max(h.send),
                is_rev,
            ));
        }
    }
    // Merge overlapping HSPs per subject to get coarse coverage spans
    fn merge(mut v: Vec<HspBlock>) -> Vec<MergedSpan> {
        if v.is_empty() {
            return Vec::new();
        }
        v.sort_by_key(|t| t.0);
        let mut merged: Vec<(usize, usize, usize, usize)> = Vec::new();
        for (a, b, sa, sb, _rev) in v {
            if let Some(last) = merged.last_mut() {
                if a <= last.1 + 5 {
                    // extend both query and subject bounds
                    last.1 = last.1.max(b);
                    last.3 = last.3.max(sb);
                    last.2 = last.2.min(sa);
                } else {
                    merged.push((a, b, sa, sb));
                }
            } else {
                merged.push((a, b, sa, sb));
            }
        }
        merged
    }
    #[derive(Clone)]
    struct SubjectWorking {
        id: String,
        spans_all: Vec<MergedSpan>,
        spans_strong: Vec<MergedSpan>,
        qcov: f64,
        scov: f64,
        orientation: String,
        strand_conf: f64,
        orientation_conflict: bool,
        order_conflict: bool,
    }

    let mut warning_set: HashSet<String> = HashSet::new();
    let mut subjects_working: Vec<SubjectWorking> = Vec::new();
    for (sid, spans) in by_subject.into_iter() {
        let mut merged_all = merge(spans);
        merged_all.sort_by_key(|t| t.0);
        let strong_raw = by_subject_strong.get(&sid).cloned().unwrap_or_default();
        let mut merged_strong = merge(strong_raw);
        merged_strong.sort_by_key(|t| t.0);
        let spans_for_cov = if merged_strong.is_empty() {
            &merged_all
        } else {
            &merged_strong
        };
        let qcovered: usize = spans_for_cov
            .iter()
            .map(|(qa, qb, _, _)| qb.saturating_sub(*qa) + 1)
            .sum();
        let qcov = (qcovered as f64) / (qlen as f64);
        let slen = *subj_slen.get(&sid).unwrap_or(&1);
        let scov = if slen > 0 {
            let scov_cov: usize = spans_for_cov
                .iter()
                .map(|(_, _, sa, sb)| sb.saturating_sub(*sa) + 1)
                .sum();
            (scov_cov as f64) / (slen as f64)
        } else {
            0.0
        };
        let (fwd, rev) = subj_orient.get(&sid).copied().unwrap_or((0, 0));
        let (orientation, strand_conf, orientation_conflict) =
            determine_orientation(fwd, rev, orient_maj);
        let order_conflict = detect_order_conflict(&merged_all, &orientation);
        if orientation_conflict {
            warning_set.insert(format!("StructVarOrientationConflict:{}", sid));
        }
        if order_conflict {
            warning_set.insert(format!("StructVarOrderConflict:{}", sid));
        }
        subjects_working.push(SubjectWorking {
            id: sid.clone(),
            spans_all: merged_all,
            spans_strong: merged_strong,
            qcov,
            scov,
            orientation,
            strand_conf,
            orientation_conflict,
            order_conflict,
        });
    }
    let mut spans: Vec<(usize, usize)> = subjects_working
        .iter()
        .flat_map(|s| s.spans_all.iter().map(|(a, b, _, _)| (*a, *b)))
        .collect();
    spans.sort_by_key(|t| t.0);
    sv.spans = spans;
    // fusion: disjoint query spans with a large gap, using strong spans when available
    #[derive(Clone)]
    struct FusionCandidate {
        idx: usize,
        min: usize,
        max: usize,
    }
    let mut fusion_candidates: Vec<FusionCandidate> = subjects_working
        .iter()
        .enumerate()
        .filter_map(|(idx, subj)| {
            if subj.qcov < min_subj_cov || subj.orientation_conflict || subj.scov < 0.20 {
                return None;
            }
            let spans = if subj.spans_strong.is_empty() {
                &subj.spans_all
            } else {
                &subj.spans_strong
            };
            let (min, max) = span_bounds(spans);
            if min == 0 && max == 0 {
                return None;
            }
            Some(FusionCandidate { idx, min, max })
        })
        .collect();
    fusion_candidates.sort_by_key(|c| (c.min, c.max, c.idx));
    let mut by_max = fusion_candidates.clone();
    by_max.sort_by_key(|c| (c.max, c.min, c.idx));
    let mut max_idx = 0usize;
    let mut best_left: Option<FusionCandidate> = None;
    'fusion: for right in &fusion_candidates {
        let threshold = right.min.saturating_sub(fusion_min_gap);
        while max_idx < by_max.len() && by_max[max_idx].max <= threshold {
            let cand = by_max[max_idx].clone();
            if best_left
                .as_ref()
                .map_or(true, |best| cand.max > best.max)
            {
                best_left = Some(cand);
            }
            max_idx += 1;
        }
        if let Some(left) = best_left.as_ref() {
            let gap = right.min - left.max;
            if gap >= fusion_min_gap {
                let left_subj = &subjects_working[left.idx];
                let right_subj = &subjects_working[right.idx];
                sv.fusion_possible = true;
                sv.fusion_gap = Some(gap);
                sv.fusion_left_len = Some(left.max - left.min + 1);
                sv.fusion_right_len = Some(right.max - right.min + 1);
                sv.fusion_subjects = Some((left_subj.id.clone(), right_subj.id.clone()));
                sv.fusion_cover_fracs =
                    Some((left_subj.scov.min(1.0), right_subj.scov.min(1.0)));
                break 'fusion;
            }
        }
    }
    // split: coverage delta from the best hit if available
    if let Some(best_subject) = subjects_working.iter().max_by(|a, b| {
        a.qcov
            .partial_cmp(&b.qcov)
            .unwrap_or(std::cmp::Ordering::Equal)
    }) {
        if (best_subject.qcov - best_subject.scov).abs() >= split_delta {
            sv.split_possible = true;
        }
    }
    // duplication: within a subject, two strong disjoint spans with small gap suggest internal duplication
    'dup: for subj in &subjects_working {
        if subj.qcov < min_subj_cov {
            continue;
        }
        let v = if subj.spans_strong.is_empty() {
            &subj.spans_all
        } else {
            &subj.spans_strong
        };
        for i in 0..v.len().saturating_sub(1) {
            let (a1, b1, _sa1, _sb1) = v[i];
            let (a2, b2, _sa2, _sb2) = v[i + 1];
            if a2 > b1 {
                let gap = a2 - b1;
                if gap <= dup_max_gap
                    && (b1 - a1 + 1) >= min_strong_len
                    && (b2 - a2 + 1) >= min_strong_len
                {
                    sv.duplication_possible = true;
                    break 'dup;
                }
            }
        }
    }

    // Final classification preference order (fusion > split > duplication)
    sv.classification = if sv.fusion_possible {
        "FusionPossible".into()
    } else if sv.split_possible {
        "SplitPossible".into()
    } else if sv.duplication_possible {
        "InternalDuplicationPossible".into()
    } else {
        "None".into()
    };
    sv.subjects = subjects_working
        .iter()
        .map(|s| StructVarSubjectSummary {
            subject_id: s.id.clone(),
            query_coverage: s.qcov,
            subject_coverage: s.scov,
            span_count: s.spans_all.len(),
            orientation: s.orientation.clone(),
            strand_confidence: s.strand_conf,
            order_conflict: s.order_conflict,
        })
        .collect();
    let mut warnings: Vec<String> = warning_set.into_iter().collect();
    warnings.sort();
    sv.warnings = warnings;
    sv
}

fn span_bounds(spans: &[MergedSpan]) -> (usize, usize) {
    if spans.is_empty() {
        return (0, 0);
    }
    let min = spans.first().map(|t| t.0).unwrap_or(0);
    let max = spans.last().map(|t| t.1).unwrap_or(0);
    (min, max)
}

fn determine_orientation(fwd: usize, rev: usize, thresh: f64) -> (String, f64, bool) {
    let total = fwd + rev;
    if total == 0 {
        return ("Unknown".into(), 0.0, true);
    }
    if rev == 0 {
        let conf = (fwd as f64) / (total as f64);
        ("Forward".into(), conf, conf < thresh)
    } else if fwd == 0 {
        let conf = (rev as f64) / (total as f64);
        ("Reverse".into(), conf, conf < thresh)
    } else {
        let conf = (fwd.max(rev) as f64) / (total as f64);
        ("Mixed".into(), conf, true)
    }
}

fn detect_order_conflict(spans: &[MergedSpan], orientation: &str) -> bool {
    if spans.len() <= 1 {
        return false;
    }
    match orientation {
        "Forward" => {
            let mut last = spans[0].2;
            for span in spans.iter().skip(1) {
                if span.2 < last {
                    return true;
                }
                last = span.2;
            }
            false
        }
        "Reverse" => {
            let mut last = spans[0].2;
            for span in spans.iter().skip(1) {
                if span.2 > last {
                    return true;
                }
                last = span.2;
            }
            false
        }
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diamond::HitSource;
    fn mk(
        q: &str,
        s: &str,
        qs: usize,
        qe: usize,
        qlen: usize,
        sstart: usize,
        send: usize,
    ) -> DiamondHitRow {
        DiamondHitRow {
            qseqid: q.into(),
            sseqid: s.into(),
            bitscore: 200.0,
            evalue: "1e-20".into(),
            length: (qe.max(qs) - qe.min(qs) + 1),
            qcov: ((qe.max(qs) - qe.min(qs) + 1) as f64) / (qlen as f64),
            scov: 0.8,
            pident: 60.0,
            qstart: qs,
            qend: qe,
            sstart,
            send,
            qlen,
            slen: 200,
            source: HitSource::SwissProt,
            staxid: Some(1),
            lineage: vec!["root".into(), "TestClass".into()],
        }
    }
    #[test]
    fn fusion_two_subjects_far_apart() {
        // Two strong HSPs on different subjects, well-separated in query
        let hits = vec![
            mk("q", "A", 50, 300, 1000, 10, 100),
            mk("q", "B", 700, 950, 1000, 20, 180),
        ];
        let mut th = StructVarThresholds::default();
        th.min_subject_cov = 0.20;
        let sv = analyze(&hits, &th);
        assert!(sv.fusion_possible);
        assert_eq!(sv.classification, "FusionPossible");
        assert!(sv.fusion_gap.unwrap() >= 50);
    }
    #[test]
    fn duplication_adjacent_within_subject() {
        // Two strong blocks separated by a small gap (<=20 aa) within the same subject
        let hits = vec![
            mk("q", "A", 100, 320, 1000, 10, 200),
            mk("q", "A", 340, 540, 1000, 220, 420),
        ];
        let sv = analyze(&hits, &StructVarThresholds::default());
        assert!(sv.duplication_possible);
    }

    #[test]
    fn split_possible_when_subject_is_partial() {
        let mut hit = mk("q", "A", 5, 600, 1000, 10, 700);
        hit.scov = 0.3;
        hit.qcov = 0.9;
        let sv = analyze(&[hit], &StructVarThresholds::default());
        assert!(sv.split_possible);
        assert_eq!(sv.classification, "SplitPossible");
    }

    #[test]
    fn orientation_conflict_triggers_warning() {
        let hits = vec![
            mk("q", "A", 50, 150, 800, 20, 120),
            mk("q", "A", 200, 300, 800, 300, 200),
        ];
        let sv = analyze(&hits, &StructVarThresholds::default());
        assert!(!sv.subjects.is_empty());
        assert!(sv
            .warnings
            .iter()
            .any(|w| w.contains("StructVarOrientationConflict")));
    }
}
