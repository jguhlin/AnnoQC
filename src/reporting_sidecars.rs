use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::hmmer_support::DomainsArchDebugRow;
use crate::structvar;

pub fn write_domains_arch_debug(
    out_dir: &str,
    rows: &[DomainsArchDebugRow],
) -> Result<(), Box<dyn std::error::Error>> {
    if rows.is_empty() {
        return Ok(());
    }

    let dbg_path = Path::new(out_dir).join("domains_arch_debug.csv");
    let mut fdbg = File::create(dbg_path)?;
    writeln!(
        fdbg,
        "gene_id,panel_size,refs_with_domains,query_domains,core_count,accessory_count,overlap_core,overlap_accessory,recall_core,precision_acc,extras_pen,score"
    )?;
    for (gid, ps, rwd, qd, cc, ac, oc, oa, rc, pa, ep, sc) in rows {
        writeln!(
            fdbg,
            "{},{},{},{},{},{},{},{},{:.4},{:.4},{:.4},{:.4}",
            gid, ps, rwd, qd, cc, ac, oc, oa, rc, pa, ep, sc
        )?;
    }
    Ok(())
}

pub fn write_structvar_summary(
    out_dir: &str,
    structvar_map: &HashMap<String, structvar::StructVar>,
) -> Result<(), Box<dyn std::error::Error>> {
    if structvar_map.is_empty() {
        return Ok(());
    }

    let mut counts: HashMap<&str, usize> = Default::default();
    let mut gaps: Vec<usize> = Vec::new();
    let mut covs: Vec<f64> = Vec::new();
    for sv in structvar_map.values() {
        let k = sv.classification.as_str();
        *counts.entry(k).or_insert(0) += 1;
        if let Some(g) = sv.fusion_gap {
            gaps.push(g);
        }
        if let Some((a, b)) = sv.fusion_cover_fracs {
            covs.push(a);
            covs.push(b);
        }
    }
    covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    gaps.sort_unstable();
    let pct = |v: &[usize], p: f64| -> usize {
        if v.is_empty() {
            0
        } else {
            let i = ((p * (v.len() as f64)).clamp(0.0, (v.len() - 1) as f64)) as usize;
            v[i]
        }
    };
    let pctf = |v: &[f64], p: f64| -> f64 {
        if v.is_empty() {
            0.0
        } else {
            let i = ((p * (v.len() as f64)).clamp(0.0, (v.len() - 1) as f64)) as usize;
            v[i]
        }
    };
    let summ = serde_json::json!({
        "counts": counts,
        "gap_percentiles": {"p10": pct(&gaps,0.10), "p50": pct(&gaps,0.50), "p90": pct(&gaps,0.90)},
        "cov_percentiles": {"p10": format!("{:.3}", pctf(&covs,0.10)), "p50": format!("{:.3}", pctf(&covs,0.50)), "p90": format!("{:.3}", pctf(&covs,0.90))},
        "suggested_cutoffs": {"fusion_min_gap": pct(&gaps,0.10).max(50), "min_subject_cov": pctf(&covs,0.10)},
    });
    std::fs::write(
        Path::new(out_dir).join("structvar_summary.json"),
        serde_json::to_string_pretty(&summ)?,
    )?;
    Ok(())
}
