#[derive(Debug, Clone, Default)]
pub struct IntrinsicMetrics {
    pub ambiguous_fraction: f64,
    pub max_homopolymer: usize,
    pub low_complexity_fraction: f64,
    pub low_complexity_windows: usize,
    #[allow(dead_code)] pub start_methionine: bool,
    #[allow(dead_code)] pub alt_start_pos: Option<usize>,
    #[allow(dead_code)] pub internal_stop_count: usize,
    #[allow(dead_code)] pub terminal_stop: bool,
    pub orf_start_score: f64,
}

#[allow(clippy::needless_range_loop)]
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
                if (0..26).contains(&u) {
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

    // ORF start/stop heuristics on protein sequence
    let mut internal_stop_count = 0usize;
    for &b in seq {
        if b == b'*' { internal_stop_count += 1; }
    }
    let terminal_stop = seq.last().copied() == Some(b'*');
    if terminal_stop && internal_stop_count>0 { internal_stop_count -= 1; }
    let start_methionine = seq.first().copied() == Some(b'M');
    let k = 10usize.min(seq.len());
    let mut alt_start_pos = None;
    if !start_methionine {
        for i in 1..k {
            if seq[i] == b'M' { alt_start_pos = Some(i); break; }
        }
    }
    let mut orf_start_score = if start_methionine { 1.0 } else if alt_start_pos.is_some() { 0.8 } else { 0.5 };
    if terminal_stop { orf_start_score = (orf_start_score - 0.2f64).max(0.0); }
    if internal_stop_count > 0 { orf_start_score = (orf_start_score - 0.3f64).max(0.0); }

    IntrinsicMetrics {
        ambiguous_fraction,
        max_homopolymer: max_run,
        low_complexity_fraction,
        low_complexity_windows: low_windows,
        start_methionine,
        alt_start_pos,
        internal_stop_count,
        terminal_stop,
        orf_start_score,
    }
}
