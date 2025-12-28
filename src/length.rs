#[derive(Debug, Clone, PartialEq)]
pub struct LengthConsistency {
    pub score: f64,
    pub z: f64,
    pub ratio: f64,
    pub class_: String,
}

const RATIO_TRUNCATED: f64 = 0.8;
const RATIO_EXTENDED: f64 = 1.2;

fn median_unstable(values: &mut [f64]) -> f64 {
    let n = values.len();
    if n == 0 {
        return 0.0;
    }
    let mid = n / 2;
    if n % 2 == 1 {
        values.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
        values[mid]
    } else {
        values.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
        let upper = values[mid];
        let lower = values[..mid]
            .iter()
            .copied()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap_or(upper);
        (lower + upper) / 2.0
    }
}

pub fn compute_length_consistency(
    query_len: usize,
    subject_lens: &[usize],
) -> Option<LengthConsistency> {
    if subject_lens.is_empty() {
        return None;
    }
    let mut vals: Vec<f64> = subject_lens
        .iter()
        .copied()
        .filter(|&x| x > 0)
        .map(|x| x as f64)
        .collect();
    if vals.is_empty() {
        return None;
    }
    let med = median_unstable(&mut vals);
    let mut devs = Vec::with_capacity(vals.len());
    devs.extend(vals.iter().map(|v| (v - med).abs()));
    let mad = median_unstable(&mut devs);
    let madn = (1.4826 * mad).max(1.0);
    let lq = query_len as f64;
    let z = if med > 0.0 { (lq - med) / madn } else { 0.0 };
    let ratio = if med > 0.0 { lq / med } else { 0.0 };
    let score = (-(z.abs()) / 2.0).exp().clamp(0.0, 1.0);
    let class_ = if ratio < RATIO_TRUNCATED {
        "LikelyNTruncated"
    } else if ratio > RATIO_EXTENDED {
        "LikelyNExtended"
    } else {
        "InRange"
    }
    .to_string();
    Some(LengthConsistency {
        score,
        z,
        ratio,
        class_,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn length_consistency_basic() {
        let q = 1000usize;
        let sl = vec![990, 995, 1005, 1010, 1000];
        let lc = compute_length_consistency(q, &sl).unwrap();
        assert!(lc.score > 0.95);
        assert_eq!(lc.class_, "InRange");
        // Clearly truncated
        let q2 = 700usize;
        let lc2 = compute_length_consistency(q2, &sl).unwrap();
        assert_eq!(lc2.class_, "LikelyNTruncated");
        // Clearly extended
        let q3 = 1300usize;
        let lc3 = compute_length_consistency(q3, &sl).unwrap();
        assert_eq!(lc3.class_, "LikelyNExtended");
    }
}
