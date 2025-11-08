#[derive(Debug, Clone, PartialEq)]
pub struct LengthConsistency {
    pub score: f64,
    pub z: f64,
    pub ratio: f64,
    pub class_: String,
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
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = vals.len();
    let med = if n % 2 == 1 {
        vals[n / 2]
    } else {
        (vals[n / 2 - 1] + vals[n / 2]) / 2.0
    };
    let mut devs: Vec<f64> = vals.iter().map(|v| (v - med).abs()).collect();
    devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = devs.len();
    let mad = if m == 0 {
        0.0
    } else if m % 2 == 1 {
        devs[m / 2]
    } else {
        (devs[m / 2 - 1] + devs[m / 2]) / 2.0
    };
    let madn = (1.4826 * mad).max(1.0);
    let lq = query_len as f64;
    let z = if med > 0.0 { (lq - med) / madn } else { 0.0 };
    let ratio = if med > 0.0 { lq / med } else { 0.0 };
    let score = (-(z.abs()) / 2.0).exp().clamp(0.0, 1.0);
    let class_ = if ratio < 0.8 {
        "LikelyNTruncated"
    } else if ratio > 1.2 {
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
