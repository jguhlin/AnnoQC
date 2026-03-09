use std::collections::HashMap;
use std::path::Path;

use crate::ecs::{GeneMetrics, RenderSummary};
use crate::orchestration::{
    RunCounts, RunFeatures, RunMetrics, RunSummary, RunThroughput, StepDuration,
};
use crate::refprot;
use crate::taxonomy::{TaxonomyEvidence, TaxonomyResolver};

pub fn build_steps(
    diamond_secs: f64,
    ecs_secs: f64,
    intrinsic_secs: f64,
    consensus_secs: f64,
    emit_secs: f64,
) -> Vec<StepDuration<'static>> {
    vec![
        StepDuration {
            name: "diamond",
            seconds: diamond_secs,
        },
        StepDuration {
            name: "ecs",
            seconds: ecs_secs,
        },
        StepDuration {
            name: "intrinsic",
            seconds: intrinsic_secs,
        },
        StepDuration {
            name: "consensus",
            seconds: consensus_secs,
        },
        StepDuration {
            name: "emit_outputs",
            seconds: emit_secs,
        },
    ]
}

pub fn build_base_totals(
    total_genes: usize,
    hmmer_genes: u64,
    hmmer_hits_total: u64,
    mafft_alignments: u64,
    diamond_mode: &str,
    backfill_used: u64,
    backfill_added_total: u64,
    refprot_used_panels: u64,
) -> HashMap<&'static str, serde_json::Value> {
    let mut totals: HashMap<&'static str, serde_json::Value> = Default::default();
    totals.insert("total_genes", serde_json::json!(total_genes));
    totals.insert("hmmer_genes", serde_json::json!(hmmer_genes));
    totals.insert("hmmer_hits_total", serde_json::json!(hmmer_hits_total));
    totals.insert("mafft_alignments", serde_json::json!(mafft_alignments));
    totals.insert("diamond_mode", serde_json::json!(diamond_mode));
    totals.insert(
        "cluster_backfill_used_panels",
        serde_json::json!(backfill_used),
    );
    totals.insert(
        "cluster_backfill_added_total",
        serde_json::json!(backfill_added_total),
    );
    totals.insert(
        "refprot_panels_with_hits",
        serde_json::json!(refprot_used_panels),
    );
    totals
}

pub struct RefprotSelectionConfig<'a> {
    pub enabled: bool,
    pub readme_path: Option<&'a str>,
    pub max_scopes: usize,
    pub max_proteomes: usize,
}

pub fn augment_totals_from_taxonomy(
    out_dir: &str,
    totals: &mut HashMap<&'static str, serde_json::Value>,
    metrics: &[GeneMetrics],
    taxsum_map: &HashMap<String, Option<TaxonomyEvidence>>,
    resolver: Option<&TaxonomyResolver>,
    refprot_cfg: Option<RefprotSelectionConfig<'_>>,
) {
    let mut status_counts: HashMap<String, u64> = HashMap::new();
    let mut superkingdom: HashMap<String, u64> = HashMap::new();
    let mut class_counts: HashMap<String, u64> = HashMap::new();
    let mut species_counts: HashMap<String, u64> = HashMap::new();
    let mut class_scope_counts: HashMap<String, u64> = HashMap::new();
    for m in metrics {
        if let Some(Some(ev)) = taxsum_map.get(&m.gene_id) {
            *status_counts.entry(format!("{}", ev.detail)).or_default() += 1;
            if let Some(cons) = &ev.consensus {
                if cons.lineage.len() > 1 {
                    *superkingdom.entry(cons.lineage[1].clone()).or_default() += 1;
                }
                let mut class_label: Option<String> = None;
                for (i, tid) in cons.lineage_ids.iter().enumerate() {
                    if let Some(r) = resolver {
                        if r.rank_of(*tid) == Some("class") {
                            class_label = cons.lineage.get(i).cloned();
                            break;
                        }
                    }
                }
                if let Some(lbl) = class_label {
                    *class_counts.entry(lbl.clone()).or_default() += 1;
                    *class_scope_counts.entry(lbl).or_default() += 1;
                }
                if let Some(last) = cons.lineage.last() {
                    *species_counts.entry(last.clone()).or_default() += 1;
                }
            }
        }
    }
    totals.insert(
        "taxonomy_status_counts",
        serde_json::to_value(status_counts).unwrap(),
    );
    totals.insert(
        "taxonomy_superkingdom_top",
        serde_json::to_value(superkingdom).unwrap(),
    );
    totals.insert(
        "taxonomy_class_top",
        serde_json::to_value(class_counts).unwrap(),
    );
    totals.insert(
        "taxonomy_species_top",
        serde_json::to_value(species_counts).unwrap(),
    );

    if let (Some(cfg), Some(resolver)) = (refprot_cfg, resolver) {
        if cfg.enabled {
            if let Some(readme) = cfg.readme_path {
                if let Ok(entries) = refprot::parse_readme(readme) {
                    let mut v: Vec<(String, u64)> = class_scope_counts.into_iter().collect();
                    v.sort_by(|a, b| b.1.cmp(&a.1));
                    let top_classes: Vec<String> = v
                        .into_iter()
                        .map(|x| x.0)
                        .take(cfg.max_scopes.max(1))
                        .collect();
                    let mut class_taxids: Vec<u32> = Vec::new();
                    for e in &entries {
                        let (lineage, lids, _) = resolver.reconstruct_lineage_public(e.taxid);
                        for (i, name) in lineage.iter().enumerate() {
                            if top_classes.iter().any(|c| c == name) {
                                if let Some(tid) = lids.get(i) {
                                    class_taxids.push(*tid);
                                }
                                break;
                            }
                        }
                    }
                    class_taxids.sort();
                    class_taxids.dedup();
                    let selected = refprot::select_by_taxon(
                        &entries,
                        &class_taxids,
                        resolver,
                        cfg.max_proteomes.max(1),
                    );
                    let mut lines = Vec::new();
                    for e in &selected {
                        lines.push(format!("{}\t{}\t{}", e.proteome_id, e.taxid, e.organism));
                    }
                    let _ = std::fs::write(
                        Path::new(out_dir).join("refprot_selected.txt"),
                        lines.join("\n"),
                    );
                    totals.insert("refprot_selected_count", serde_json::json!(selected.len()));
                }
            }
        }
    }
}

pub struct RunSummaryInputs<'a> {
    pub schema_version: &'a str,
    pub total_genes: usize,
    pub render_summary: &'a RenderSummary,
    pub aligner: String,
    pub report_format: String,
    pub resume: bool,
    pub taxonomy_enabled: bool,
    pub hmmer_enabled: bool,
    pub alignment_enabled: bool,
    pub genomic_enabled: bool,
    pub plugins: usize,
    pub rules: usize,
    pub nucleotide: bool,
    pub calibration: String,
    pub timings: Vec<StepDuration<'a>>,
}

pub fn build_run_summary(inputs: RunSummaryInputs<'_>) -> RunSummary<'_> {
    let total_secs: f64 = inputs.timings.iter().map(|s| s.seconds).sum();
    let genes_per_sec = if total_secs > 0.0 {
        inputs.total_genes as f64 / total_secs
    } else {
        0.0
    };
    RunSummary {
        schema_version: inputs.schema_version,
        counts: RunCounts {
            total_genes: inputs.total_genes,
            rendered_genes: inputs.render_summary.total,
            high: inputs.render_summary.high,
            low: inputs.render_summary.low,
            x: inputs.render_summary.x,
            high_complete: inputs.render_summary.high_complete,
            high_fragmented: inputs.render_summary.high_fragmented,
            low_novel: inputs.render_summary.low_novel,
            low_artifact: inputs.render_summary.low_artifact,
        },
        throughput: RunThroughput {
            total_seconds: total_secs,
            genes_per_second: genes_per_sec,
        },
        features: RunFeatures {
            aligner: inputs.aligner,
            report_format: inputs.report_format,
            resume: inputs.resume,
            taxonomy_enabled: inputs.taxonomy_enabled,
            hmmer_enabled: inputs.hmmer_enabled,
            alignment_enabled: inputs.alignment_enabled,
            genomic_enabled: inputs.genomic_enabled,
            plugins: inputs.plugins,
            rules: inputs.rules,
            nucleotide: inputs.nucleotide,
            calibration: inputs.calibration,
        },
        timings: inputs.timings,
        errors: Vec::new(),
    }
}

pub fn build_run_metrics<'a>(
    schema_version: &'a str,
    steps: Vec<StepDuration<'a>>,
    totals: HashMap<&'a str, serde_json::Value>,
) -> RunMetrics<'a> {
    RunMetrics {
        schema_version,
        steps,
        totals,
    }
}
