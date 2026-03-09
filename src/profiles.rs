use crate::config_types::{ScoringConfigOverride, ScoringThresholds};
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Serialize, Deserialize, Default)]
pub enum Profile {
    #[default]
    Standard,
    Strict,
    Discovery,
    Bacteria,
}

#[allow(dead_code)]
pub struct ProfileConfig {
    pub weights: HashMap<String, f64>,
    pub thresholds: ScoringThresholds,
    // Future: alignment settings, etc.
}

const STANDARD_WEIGHTS: [(&str, f64); 9] = [
    ("homology", 0.6),
    ("intrinsic", 0.4),
    ("taxonomy", 0.0),
    ("domains", 0.0),
    ("length", 0.0),
    ("orphan", 0.0),
    ("subject_cov", 0.0),
    ("termini", 0.0),
    ("divergence", 0.0),
];

const STRICT_WEIGHTS: [(&str, f64); 9] = [
    ("homology", 0.5),
    ("intrinsic", 0.2),
    ("taxonomy", 0.1),
    ("domains", 0.1),
    ("length", 0.1),
    ("orphan", 0.0),
    ("subject_cov", 0.0),
    ("termini", 0.0),
    ("divergence", 0.0),
];

const DISCOVERY_WEIGHTS: [(&str, f64); 9] = [
    ("homology", 0.3),
    ("intrinsic", 0.6),
    ("taxonomy", 0.0),
    ("domains", 0.1),
    ("length", 0.0),
    ("orphan", 0.0),
    ("subject_cov", 0.0),
    ("termini", 0.0),
    ("divergence", 0.0),
];

const BACTERIA_WEIGHTS: [(&str, f64); 9] = [
    ("homology", 0.7),
    ("intrinsic", 0.3),
    ("taxonomy", 0.0),
    ("domains", 0.0),
    ("length", 0.0),
    ("orphan", 0.0),
    ("subject_cov", 0.0),
    ("termini", 0.0),
    ("divergence", 0.0),
];

#[allow(dead_code)]
fn build_weights(weights: &[(&str, f64)]) -> HashMap<String, f64> {
    let mut map = HashMap::with_capacity(weights.len());
    for (name, value) in weights {
        map.insert((*name).to_string(), *value);
    }
    map
}

fn validate_weights_sum(weights: &[(&str, f64)]) {
    let sum: f64 = weights.iter().map(|(_, v)| *v).sum();
    debug_assert!(
        (sum - 1.0).abs() < 1e-6,
        "profile weights must sum to 1.0, got {}",
        sum
    );
}

fn validate_thresholds(thresholds: &ScoringThresholds) {
    if let (Some(high), Some(medium)) = (thresholds.high, thresholds.medium) {
        debug_assert!(
            high >= medium,
            "thresholds must be monotonic: high {} < medium {}",
            high,
            medium
        );
    }
}

impl Profile {
    fn weights(&self) -> &'static [(&'static str, f64)] {
        match self {
            Profile::Standard => &STANDARD_WEIGHTS,
            Profile::Strict => &STRICT_WEIGHTS,
            Profile::Discovery => &DISCOVERY_WEIGHTS,
            Profile::Bacteria => &BACTERIA_WEIGHTS,
        }
    }

    fn thresholds(&self) -> ScoringThresholds {
        match self {
            Profile::Standard => ScoringThresholds {
                high: Some(0.8),
                medium: Some(0.5),
            },
            Profile::Strict => ScoringThresholds {
                high: Some(0.9),
                medium: Some(0.7),
            },
            Profile::Discovery => ScoringThresholds {
                high: Some(0.8),
                medium: Some(0.5),
            },
            Profile::Bacteria => ScoringThresholds {
                high: Some(0.85),
                medium: Some(0.6),
            },
        }
    }

    #[allow(dead_code)]
    pub fn config(&self) -> ProfileConfig {
        let weights = self.weights();
        let thresholds = self.thresholds();
        validate_weights_sum(weights);
        validate_thresholds(&thresholds);

        ProfileConfig {
            weights: build_weights(weights),
            thresholds,
        }
    }

    /// Applies profile defaults without overriding explicit config values.
    /// Thresholds only fill missing fields; weights only fill missing keys.
    pub fn apply_to(&self, cfg: &mut ScoringConfigOverride) {
        let weights = self.weights();
        let thresholds = self.thresholds();
        validate_weights_sum(weights);
        validate_thresholds(&thresholds);

        // Only apply if not already set in config
        if cfg.thresholds.is_none() {
            cfg.thresholds = Some(thresholds);
        } else if let Some(ref mut th) = cfg.thresholds {
            if th.high.is_none() {
                th.high = thresholds.high;
            }
            if th.medium.is_none() {
                th.medium = thresholds.medium;
            }
        }

        // Merge weights: Profile acts as "Smart Defaults"
        // If config has weight, keep it. Else use profile.
        for (name, value) in weights {
            cfg.weights.entry((*name).to_string()).or_insert(*value);
        }
    }
}
