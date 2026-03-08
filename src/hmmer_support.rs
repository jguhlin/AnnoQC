use std::collections::HashMap;

use crate::ecs::{run_heavy_pipelines, GeneMetrics, HeavyPipelineConfig, HmmerPipelineConfig};
use crate::hmmer;
use crate::hmmer::HmmscanSummary;
use crate::mafft::load_sequences_by_ids;
use crate::taxonomy;

#[derive(Clone, Debug, Default)]
pub struct HmmerOverrideInputs {
    pub top_n: Option<usize>,
    pub threads: Option<usize>,
    pub query_ievalue: Option<f64>,
    pub ref_ievalue: Option<f64>,
    pub orphan_analysis: Option<bool>,
}

pub struct HmmerSetup {
    pub config: Option<HmmerPipelineConfig>,
    pub bin_owned: Option<String>,
    pub db_owned: Option<String>,
    pub top_n: usize,
    pub thread_cap: usize,
    pub ref_ievalue: Option<f64>,
    pub orphan_analysis_enabled: bool,
    pub job_count: usize,
}

pub type DomainsArchDebugRow = (
    String,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
);

pub struct HmmerPostProcessOutput {
    pub domains_arch_map: HashMap<String, f64>,
    pub domains_arch_dbg: Vec<DomainsArchDebugRow>,
    pub orphan_map: HashMap<String, hmmer::OrphanAnalysis>,
}

pub fn build_hmmer_setup(
    hmmscan_bin: Option<&str>,
    pfam_db: Option<&str>,
    items: Vec<(String, Vec<u8>)>,
    threads_default: usize,
    cli: &HmmerOverrideInputs,
    config: &HmmerOverrideInputs,
    disable_orphan_analysis: bool,
) -> HmmerSetup {
    let top_n = cli.top_n.or(config.top_n).unwrap_or(5);
    let thread_cap = cli.threads.or(config.threads).unwrap_or(threads_default);
    let query_ievalue = cli.query_ievalue.or(config.query_ievalue);
    let ref_ievalue = cli.ref_ievalue.or(config.ref_ievalue);
    let orphan_cfg = config.orphan_analysis.unwrap_or(true);
    let orphan_analysis_enabled = orphan_cfg && !disable_orphan_analysis;

    let mut cfg = None;
    let mut bin_owned = None;
    let mut db_owned = None;
    let job_count = items.len();

    if let (Some(hmm), Some(db)) = (hmmscan_bin, pfam_db) {
        bin_owned = Some(hmm.to_string());
        db_owned = Some(db.to_string());
        if !items.is_empty() {
            cfg = Some(HmmerPipelineConfig {
                hmmscan_bin: hmm.to_string(),
                db_path: db.to_string(),
                items,
                threads_per_job: 1,
                max_jobs: thread_cap.max(1),
                top_n,
                max_ievalue: query_ievalue,
            });
        }
    }

    HmmerSetup {
        config: cfg,
        bin_owned,
        db_owned,
        top_n,
        thread_cap,
        ref_ievalue,
        orphan_analysis_enabled,
        job_count,
    }
}

pub struct HmmerPostProcessInputs<'a> {
    pub hmm_bin: Option<&'a str>,
    pub db_path: Option<&'a str>,
    pub reference_fasta: Option<&'a str>,
    pub panel_map: &'a HashMap<String, Vec<String>>,
    pub hmmsum_map: &'a HashMap<String, HmmscanSummary>,
    pub metrics: &'a [GeneMetrics],
    pub cfg_threads: usize,
    pub log_json: bool,
    pub hmmer_thread_cap: usize,
    pub hmmer_top_n: usize,
    pub hmmer_ref_ievalue: Option<f64>,
    pub orphan_analysis_enabled: bool,
    pub pfam_clans_path: Option<&'a str>,
    pub domain_order_weight: f64,
}

pub fn postprocess_hmmer(inputs: HmmerPostProcessInputs<'_>) -> HmmerPostProcessOutput {
    let mut ref_hmmsum_map: HashMap<String, HmmscanSummary> = HashMap::new();
    let mut domains_arch_map: HashMap<String, f64> = HashMap::new();
    let mut domains_arch_dbg: Vec<DomainsArchDebugRow> = Vec::new();
    let mut orphan_map: HashMap<String, hmmer::OrphanAnalysis> = HashMap::new();

    let Some(hmm_bin) = inputs.hmm_bin else {
        return HmmerPostProcessOutput {
            domains_arch_map,
            domains_arch_dbg,
            orphan_map,
        };
    };

    if inputs.orphan_analysis_enabled {
        orphan_map = inputs
            .hmmsum_map
            .iter()
            .filter(|(_, summary)| !summary.hits.is_empty())
            .map(|(gid, summary)| (gid.clone(), hmmer::analyze_orphan_domains(summary)))
            .collect();
    }

    if let (Some(ref_fasta), Some(db_path)) = (inputs.reference_fasta, inputs.db_path) {
        let all_ids: Vec<String> = inputs.panel_map.values().flat_map(|v| v.clone()).collect();
        let ref_seqs = load_sequences_by_ids(ref_fasta, &all_ids).unwrap_or_default();
        if !ref_seqs.is_empty() {
            let ref_items: Vec<(String, Vec<u8>)> = ref_seqs.into_iter().collect();
            if !ref_items.is_empty() {
                log::info!("hmmscan reference jobs queued: {}", ref_items.len());
                let reserve_threads = inputs.cfg_threads.min(2);
                let ref_cfg = HmmerPipelineConfig {
                    hmmscan_bin: hmm_bin.to_string(),
                    db_path: db_path.to_string(),
                    items: ref_items,
                    threads_per_job: 1,
                    max_jobs: inputs.hmmer_thread_cap.max(1),
                    top_n: inputs.hmmer_top_n,
                    max_ievalue: inputs.hmmer_ref_ievalue,
                };
                let ref_results = run_heavy_pipelines(HeavyPipelineConfig {
                    cpu_threads: inputs.cfg_threads,
                    reserve_threads,
                    log_json: inputs.log_json,
                    alignment: None,
                    hmmer: Some(ref_cfg),
                });
                ref_hmmsum_map = ref_results.hmmer_map;
            }
        }
    }

    let clan_map = inputs
        .pfam_clans_path
        .and_then(|p| hmmer::load_pfam_clans(p).ok());
    for g in inputs.metrics {
        if let Some(qsum) = inputs.hmmsum_map.get(&g.gene_id) {
            let qsum_c = if let Some(ref clans) = clan_map {
                hmmer::collapse_by_clan(qsum, clans)
            } else {
                qsum.clone()
            };
            let ref_ids = inputs
                .panel_map
                .get(&g.gene_id)
                .cloned()
                .unwrap_or_default();
            let ref_ids_canonical: Vec<String> = ref_ids
                .iter()
                .map(|rid| taxonomy::canonical_accession(rid))
                .collect();
            let mut ref_map_c: HashMap<String, HmmscanSummary> = HashMap::new();
            for key in &ref_ids_canonical {
                if let Some(s) = ref_hmmsum_map.get(key) {
                    let val = if let Some(ref clans) = clan_map {
                        hmmer::collapse_by_clan(s, clans)
                    } else {
                        s.clone()
                    };
                    ref_map_c.insert(key.clone(), val);
                }
            }
            if ref_map_c.is_empty() && !ref_ids_canonical.is_empty() {
                log::debug!(
                    "ref_join_empty: gene={} panel_n={} keys_example={:?} ref_map_total={}",
                    g.gene_id,
                    ref_ids_canonical.len(),
                    &ref_ids_canonical.iter().take(5).collect::<Vec<_>>(),
                    ref_hmmsum_map.len()
                );
            }
            let dbg = hmmer::domains_architecture_diagnostics(
                &qsum_c,
                &ref_ids_canonical,
                &ref_map_c,
                inputs.domain_order_weight,
            );
            domains_arch_map.insert(g.gene_id.clone(), dbg.score);
            domains_arch_dbg.push((
                g.gene_id.clone(),
                dbg.panel_size,
                dbg.refs_with_domains,
                dbg.query_domains,
                dbg.core_count,
                dbg.accessory_count,
                dbg.overlap_core,
                dbg.overlap_accessory,
                dbg.recall_core,
                dbg.precision_acc,
                dbg.extras_pen,
                dbg.score,
            ));
        }
    }

    HmmerPostProcessOutput {
        domains_arch_map,
        domains_arch_dbg,
        orphan_map,
    }
}
