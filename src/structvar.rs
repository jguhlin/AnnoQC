use crate::diamond::DiamondHitRow;
use std::collections::HashMap;

type HspBlock = (usize, usize, usize, usize, bool);
type MergedSpan = (usize, usize, usize, usize);
type SubjectHspMap<'a> = HashMap<&'a str, Vec<HspBlock>>;
type SubjectSpanVec = Vec<(String, Vec<MergedSpan>)>;

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
    let mut by_subject: SubjectHspMap<'_> = HashMap::new();
    let mut by_subject_strong: SubjectHspMap<'_> = HashMap::new();
    let mut subj_slen: HashMap<&str, usize> = HashMap::new();
    let mut subj_orient: HashMap<&str, (usize, usize)> = HashMap::new(); // (fwd, rev)
    for h in hits {
        let a = h.qstart.min(h.qend);
        let b = h.qstart.max(h.qend);
        let is_rev = h.sstart > h.send;
        by_subject.entry(&h.sseqid).or_default().push((
            a,
            b,
            h.sstart.min(h.send),
            h.sstart.max(h.send),
            is_rev,
        ));
        subj_slen.entry(&h.sseqid).or_insert(h.slen);
        let e = subj_orient.entry(&h.sseqid).or_insert((0, 0));
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
            by_subject_strong.entry(&h.sseqid).or_default().push((
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
    let mut subj_spans: SubjectSpanVec = Vec::new();
    let mut subj_spans_strong: SubjectSpanVec = Vec::new();
    for (sid, v) in by_subject.into_iter() {
        let mut v2 = merge(v);
        v2.sort_by_key(|t| t.0);
        subj_spans.push((sid.to_string(), v2));
    }
    for (sid, v) in by_subject_strong.into_iter() {
        subj_spans_strong.push((sid.to_string(), merge(v)));
    }
    let mut subj_qcov: HashMap<String, f64> = HashMap::new();
    for (sid, spans) in &subj_spans_strong {
        let covered: usize = spans
            .iter()
            .map(|(qa, qb, _, _)| qb.saturating_sub(*qa) + 1)
            .sum();
        subj_qcov.insert(sid.clone(), (covered as f64) / (qlen as f64));
    }
    // Flatten all spans for summary and diagnostics
    let mut spans: Vec<(usize, usize)> = subj_spans
        .iter()
        .flat_map(|(_, v)| v.iter().map(|(a, b, _, _)| (*a, *b)))
        .collect();
    spans.sort_by_key(|t| t.0);
    sv.spans = spans.clone();
    // fusion: two different subjects whose (merged, strong) spans are disjoint, separated by large gap
    'fusion: for i in 0..subj_spans_strong.len().saturating_sub(1) {
        for j in i + 1..subj_spans_strong.len() {
            let (ref ida, ref a) = subj_spans[i];
            let (ref idb, ref b) = subj_spans[j];
            let qcov_a = *subj_qcov.get(ida).unwrap_or(&0.0);
            let qcov_b = *subj_qcov.get(idb).unwrap_or(&0.0);
            if qcov_a < min_subj_cov || qcov_b < min_subj_cov {
                continue;
            }
            // use coarse first and last positions
            let (a_min, a_max) = (
                a.first().map(|x| x.0).unwrap_or(0),
                a.last().map(|x| x.1).unwrap_or(0),
            );
            let (b_min, b_max) = (
                b.first().map(|x| x.0).unwrap_or(0),
                b.last().map(|x| x.1).unwrap_or(0),
            );
            if a_max > 0 && b_max > 0 {
                let (left_max, right_min) = if a_min <= b_min {
                    (a_max, b_min)
                } else {
                    (b_max, a_min)
                };
                if right_min > left_max {
                    let gap = right_min - left_max;
                    if gap >= fusion_min_gap {
                        // per-subject coverage fractions (strong merged spans)
                        // compute dominant orientation (>= th.orient_majority votes) for each subject
                        let dominant = |id: &str| -> bool {
                            let (fwd, rev) = *subj_orient.get(id).unwrap_or(&(1, 0));
                            let tot = fwd + rev;
                            if tot == 0 {
                                return true;
                            }
                            let (maj, _min) = if fwd >= rev { (fwd, rev) } else { (rev, fwd) };
                            (maj as f64) / (tot as f64) >= orient_maj
                        };
                        if !dominant(ida) || !dominant(idb) {
                            continue;
                        }
                        let (cov_a, cov_b) = {
                            let sa_cov: usize = a
                                .iter()
                                .map(|(_qa, _qb, sa, sb)| sb.saturating_sub(*sa) + 1)
                                .sum();
                            let sb_cov: usize = b
                                .iter()
                                .map(|(_qa, _qb, sa, sb)| sb.saturating_sub(*sa) + 1)
                                .sum();
                            let la = *subj_slen.get(ida.as_str()).unwrap_or(&1);
                            let lb = *subj_slen.get(idb.as_str()).unwrap_or(&1);
                            ((sa_cov as f64) / (la as f64), (sb_cov as f64) / (lb as f64))
                        };
                        if cov_a < 0.20 || cov_b < 0.20 {
                            continue;
                        }
                        sv.fusion_possible = true;
                        let left_len = (left_max - if a_min <= b_min { a_min } else { b_min }) + 1;
                        let right_len =
                            (if a_min <= b_min { b_max } else { a_max }) - right_min + 1;
                        sv.fusion_gap = Some(gap);
                        sv.fusion_left_len = Some(left_len);
                        sv.fusion_right_len = Some(right_len);
                        sv.fusion_subjects = Some((ida.clone(), idb.clone()));
                        sv.fusion_cover_fracs = Some((cov_a.min(1.0), cov_b.min(1.0)));
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
        if (qcov - scov).abs() >= split_delta {
            sv.split_possible = true;
        }
    }
    // duplication: within a subject, two strong disjoint spans with small gap suggest internal duplication
    'dup: for (sid, v) in &subj_spans_strong {
        if *subj_qcov.get(sid).unwrap_or(&0.0) < min_subj_cov {
            continue;
        }
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
    sv
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
