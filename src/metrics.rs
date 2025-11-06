#[derive(Debug, Clone, Default)]
pub struct IntrinsicMetrics {
    pub ambiguous_fraction: f64,
    pub max_homopolymer: usize,
    pub low_complexity_fraction: f64,
    pub low_complexity_windows: usize,
}

pub fn compute_intrinsic(seq: &[u8]) -> IntrinsicMetrics {
    let len = seq.len().max(1);
    let amb = seq
        .iter()
        .filter(|c| {
            matches!(
                c.to_ascii_uppercase(),
                b'X' | b'B' | b'Z' | b'J' | b'U' | b'O'
            )
        })
        .count();
    let ambiguous_fraction = amb as f64 / len as f64;

    // max homopolymer length
    let mut max_run = 1usize;
    let mut run = 1usize;
    for i in 1..seq.len() {
        if seq[i] == seq[i - 1] {
            run += 1;
            max_run = max_run.max(run);
        } else {
            run = 1;
        }
    }

    // naive low-complexity: sliding window uniqueness threshold
    let w = 25usize;
    let mut low_windows = 0usize;
    let mut total_windows = 0usize;
    if seq.len() >= w {
        for i in 0..=seq.len() - w {
            total_windows += 1;
            let window = &seq[i..i + w];
            let mut mask = [false; 26];
            let mut uniq = 0usize;
            for &b in window {
                let u = (b.to_ascii_uppercase() as i32) - ('A' as i32);
                if u >= 0 && u < 26 {
                    let ui = u as usize;
                    if !mask[ui] {
                        mask[ui] = true;
                        uniq += 1;
                    }
                }
            }
            if uniq <= 6 {
                low_windows += 1;
            }
        }
    }
    let low_complexity_fraction = if total_windows > 0 {
        low_windows as f64 / total_windows as f64
    } else {
        0.0
    };

    IntrinsicMetrics {
        ambiguous_fraction,
        max_homopolymer: max_run,
        low_complexity_fraction,
        low_complexity_windows: low_windows,
    }
}
