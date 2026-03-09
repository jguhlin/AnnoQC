use std::cmp::Ordering;

use crate::scoring_thresholds::ScoreTemp;

pub(crate) fn apply_percentile_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let len = order.len();
    for (rank, idx) in order.into_iter().enumerate() {
        let percentile = if len > 1 {
            (rank as f64 + 0.5) / len as f64
        } else {
            1.0
        };
        entries[idx].final_score = percentile.clamp(0.0, 1.0);
    }
}

pub(crate) fn apply_isotonic_calibration(entries: &mut [ScoreTemp]) {
    if entries.is_empty() {
        return;
    }
    let mut order: Vec<usize> = (0..entries.len()).collect();
    order.sort_by(|&a, &b| {
        entries[a]
            .raw_score
            .partial_cmp(&entries[b].raw_score)
            .unwrap_or(Ordering::Equal)
    });
    let n = order.len();
    #[derive(Clone)]
    struct Block {
        start: usize,
        end: usize,
        sum: f64,
        weight: usize,
    }
    let mut blocks: Vec<Block> = Vec::new();
    for (rank, _idx) in order.iter().enumerate() {
        let y = if n > 1 {
            (rank as f64 + 0.5) / n as f64
        } else {
            1.0
        };
        blocks.push(Block {
            start: rank,
            end: rank,
            sum: y,
            weight: 1,
        });
        while blocks.len() >= 2 {
            let k = blocks.len() - 1;
            let prev = &blocks[k - 1];
            let curr = &blocks[k];
            let avg_prev = prev.sum / prev.weight as f64;
            let avg_curr = curr.sum / curr.weight as f64;
            if avg_prev <= avg_curr {
                break;
            }
            let merged = Block {
                start: prev.start,
                end: curr.end,
                sum: prev.sum + curr.sum,
                weight: prev.weight + curr.weight,
            };
            blocks.pop();
            blocks.pop();
            blocks.push(merged);
        }
    }
    let mut fitted = vec![0.0; n];
    for block in blocks {
        let avg = (block.sum / block.weight as f64).clamp(0.0, 1.0);
        for value in fitted
            .iter_mut()
            .take(block.end.saturating_add(1))
            .skip(block.start)
        {
            *value = avg;
        }
    }
    for (rank, idx) in order.into_iter().enumerate() {
        entries[idx].final_score = fitted[rank];
    }
}

pub(crate) fn calibration_has_min_samples(
    entries: &[ScoreTemp],
    min_samples: usize,
    min_unique: usize,
) -> bool {
    if entries.len() < min_samples {
        return false;
    }
    let mut unique: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for e in entries {
        unique.insert(e.raw_score.to_bits());
        if unique.len() >= min_unique {
            return true;
        }
    }
    false
}
