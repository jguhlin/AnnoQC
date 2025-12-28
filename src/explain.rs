use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use serde_json::Value;

/// Find and print the breakdown for `gene_id` from `qc_report.jsonl` in `out_dir`.
///
/// Returns `Ok(())` after printing the first matching record; returns an error if no
/// matching `gene_id` is found. The output includes the final score, score components,
/// homology/taxonomy summaries, warnings, and plugin penalty (if present).
pub fn explain_gene(out_dir: &str, gene_id: &str) -> Result<(), Box<dyn Error>> {
    let path = Path::new(out_dir).join("qc_report.jsonl");
    let file =
        File::open(&path).map_err(|e| format!("failed to open {}: {}", path.display(), e))?;
    let reader = BufReader::new(file);
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let val: Value = match serde_json::from_str(line) {
            Ok(v) => v,
            Err(_) => continue,
        };
        if val
            .get("type")
            .and_then(|v| v.as_str())
            .is_some_and(|v| v == "metadata")
        {
            continue;
        }
        if val
            .get("gene_id")
            .and_then(|v| v.as_str())
            .is_some_and(|v| v == gene_id)
        {
            print_breakdown(&val);
            return Ok(());
        }
    }
    Err(format!("gene_id '{}' not found in {}", gene_id, path.display()).into())
}

fn fmt_f64(val: Option<f64>) -> String {
    val.map(|v| format!("{:.4}", v))
        .unwrap_or_else(|| "NA".to_string())
}

fn print_breakdown(val: &Value) {
    let gene_id = val
        .get("gene_id")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown");
    let final_score = val.get("final_score").and_then(|v| v.as_f64());
    let raw_score = val.get("final_score_raw").and_then(|v| v.as_f64());
    println!("Gene: {}", gene_id);
    println!(
        "Final score: {} (raw {})",
        fmt_f64(final_score),
        fmt_f64(raw_score)
    );

    if let Some(comp) = val.get("score_components").and_then(|v| v.as_object()) {
        println!("Score components:");
        for (k, v) in comp {
            let s = v
                .as_f64()
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string());
            println!("  {}: {}", k, s);
        }
    }

    if let Some(homology) = val.get("homology").and_then(|v| v.as_object()) {
        let hits = homology
            .get("hits_count")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        let top_hit = homology
            .get("top_hit")
            .and_then(|v| v.as_str())
            .unwrap_or("NA");
        let top_qcov = homology.get("top_qcov").and_then(|v| v.as_f64());
        let top_scov = homology.get("top_scov").and_then(|v| v.as_f64());
        println!(
            "Homology: hits={} top_hit={} qcov={} scov={}",
            hits,
            top_hit,
            fmt_f64(top_qcov),
            fmt_f64(top_scov)
        );
    }

    if let Some(tax) = val.get("taxonomy").and_then(|v| v.as_object()) {
        let status = tax
            .get("status")
            .and_then(|v| v.as_str())
            .unwrap_or("unknown");
        let detail = tax.get("detail").and_then(|v| v.as_str()).unwrap_or("NA");
        println!("Taxonomy: status={} detail={}", status, detail);
    }

    if let Some(warnings) = val.get("warnings").and_then(|v| v.as_array()) {
        if !warnings.is_empty() {
            let mut joined = String::new();
            for warning in warnings.iter().filter_map(|v| v.as_str()) {
                if !joined.is_empty() {
                    joined.push_str(", ");
                }
                joined.push_str(warning);
            }
            if !joined.is_empty() {
                println!("Warnings: {}", joined);
            }
        }
    }

    if let Some(plugins) = val.get("plugins").and_then(|v| v.as_object()) {
        if let Some(pen) = plugins.get("total_penalty").and_then(|v| v.as_f64()) {
            println!("Plugin penalty: {:.4}", pen);
        }
    }
}
