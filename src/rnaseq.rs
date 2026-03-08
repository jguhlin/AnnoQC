use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// RNA-seq expression metrics for a single gene
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaseqMetrics {
    /// Transcripts Per Million (normalized expression)
    pub tpm: Option<f64>,
    /// Raw read count
    pub num_reads: Option<u64>,
    /// Effective transcript length
    pub effective_length: Option<f64>,
    /// Whether this gene has any expression support
    pub has_support: bool,
    /// Normalized expression score [0, 1]
    pub expression_score: f64,
}

impl Default for RnaseqMetrics {
    fn default() -> Self {
        Self {
            tpm: None,
            num_reads: None,
            effective_length: None,
            has_support: false,
            expression_score: 0.0,
        }
    }
}

/// Parse Salmon/Kallisto quant table (TSV format)
/// Expected columns: Name, Length, EffectiveLength, TPM, NumReads
pub fn parse_quant_table(path: &str) -> Result<HashMap<String, RnaseqMetrics>, String> {
    let mut map = HashMap::new();
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;

    for (idx, line) in content.lines().enumerate() {
        if idx == 0 {
            continue; // Skip header
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }

        let name = parts[0].to_string();
        let tpm: f64 = parts[3].parse().unwrap_or(0.0);
        let num_reads: u64 = parts[4].parse().unwrap_or(0);
        let eff_length: f64 = parts[2].parse().unwrap_or(0.0);

        let expression_score = normalize_tpm(tpm);

        map.insert(
            name,
            RnaseqMetrics {
                tpm: Some(tpm),
                num_reads: Some(num_reads),
                effective_length: Some(eff_length),
                has_support: num_reads > 0,
                expression_score,
            },
        );
    }

    Ok(map)
}

/// Parse precomputed expression support JSON/TSV
pub fn parse_expression_file(path: &str) -> Result<HashMap<String, RnaseqMetrics>, String> {
    let ext = std::path::Path::new(path)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    match ext {
        "json" => parse_expression_json(path),
        "tsv" | "txt" => parse_expression_tsv(path),
        _ => Err("Unsupported file format".to_string()),
    }
}

fn parse_expression_json(path: &str) -> Result<HashMap<String, RnaseqMetrics>, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let json: serde_json::Value = serde_json::from_str(&content).map_err(|e| e.to_string())?;

    let mut map = HashMap::new();

    if let Some(obj) = json.as_object() {
        for (gene_id, data) in obj {
            let data_obj = data.as_object().unwrap();
            let tpm = data_obj.get("tpm").and_then(|v| v.as_f64());
            let num_reads = data_obj.get("num_reads").and_then(|v| v.as_u64());
            let support_score = data_obj.get("support_score").and_then(|v| v.as_f64());

            let expression_score = support_score.unwrap_or_else(|| tpm.map_or(0.0, normalize_tpm));

            map.insert(
                gene_id.clone(),
                RnaseqMetrics {
                    tpm,
                    num_reads,
                    effective_length: None,
                    has_support: num_reads.map_or(false, |n| n > 0)
                        || support_score.map_or(false, |s| s > 0.0),
                    expression_score,
                },
            );
        }
    }

    Ok(map)
}

fn parse_expression_tsv(path: &str) -> Result<HashMap<String, RnaseqMetrics>, String> {
    let mut map = HashMap::new();
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;

    for (idx, line) in content.lines().enumerate() {
        if idx == 0 && line.starts_with("gene") {
            continue; // Skip header if present
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        let gene_id = parts[0].to_string();
        let tpm: f64 = parts[1].parse().unwrap_or(0.0);
        let num_reads: u64 = if parts.len() > 2 {
            parts[2].parse().unwrap_or(0)
        } else {
            0
        };
        let support_score: f64 = if parts.len() > 3 {
            parts[3].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        let expression_score = if support_score > 0.0 {
            support_score
        } else {
            normalize_tpm(tpm)
        };

        map.insert(
            gene_id,
            RnaseqMetrics {
                tpm: Some(tpm),
                num_reads: Some(num_reads),
                effective_length: None,
                has_support: num_reads > 0 || support_score > 0.0,
                expression_score,
            },
        );
    }

    Ok(map)
}

/// Normalize TPM to [0, 1] score
/// TPM < 1 → score 0, TPM >= 100 → score 1
fn normalize_tpm(tpm: f64) -> f64 {
    const MIN_TPM: f64 = 1.0;
    const MAX_TPM: f64 = 100.0;

    if tpm < MIN_TPM {
        return 0.0;
    }
    if tpm >= MAX_TPM {
        return 1.0;
    }

    // Log scale normalization
    let log_min = MIN_TPM.log10();
    let log_max = MAX_TPM.log10();
    let log_val = tpm.log10();

    (log_val - log_min) / (log_max - log_min).max(0.01)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_tpm() {
        assert_eq!(normalize_tpm(0.0), 0.0);
        assert_eq!(normalize_tpm(0.5), 0.0);
        assert_eq!(normalize_tpm(1.0), 0.0);
        assert!(normalize_tpm(10.0) > 0.0 && normalize_tpm(10.0) < 1.0);
        assert_eq!(normalize_tpm(100.0), 1.0);
        assert_eq!(normalize_tpm(1000.0), 1.0);
    }
}
