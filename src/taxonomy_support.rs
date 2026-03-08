use std::collections::HashMap;
use std::sync::Arc;

use crate::diamond;
use crate::ecs::GeneMetrics;
use crate::taxonomy::{
    TaxonomyConsensusConfig, TaxonomyDetail, TaxonomyEvidence, TaxonomyResolver,
};

#[derive(Clone)]
pub struct TaxonomyWarnings {
    pub expected_domain: Option<String>,
    pub warn_non_target_min_frac: f64,
    pub warn_non_target_min_hits: usize,
    pub warn_non_target_strong_frac: f64,
    pub warn_non_target_strong_hits: usize,
    pub warn_genus_min_frac: f64,
    pub warn_genus_min_hits: usize,
    pub low_coverage_frac: f64,
}

pub struct TaxonomySetup {
    pub enabled: bool,
    pub resolver: Option<Arc<TaxonomyResolver>>,
    pub taxsum_map: Arc<HashMap<String, Option<TaxonomyEvidence>>>,
    pub warnings: TaxonomyWarnings,
}

pub fn build_taxonomy_setup(
    metrics: &[GeneMetrics],
    taxonomy_hits_map: &HashMap<String, Vec<diamond::DiamondHitRow>>,
    stats: &HashMap<String, diamond::DiamondHitStats>,
    enable_requested: bool,
    cache_path: Option<&str>,
    reference_fasta: Option<&str>,
    taxdump_dir: Option<&str>,
    consensus: TaxonomyConsensusConfig,
    warnings: TaxonomyWarnings,
) -> Result<TaxonomySetup, String> {
    let resolver = TaxonomyResolver::from_sources(cache_path, reference_fasta, taxdump_dir)
        .map_err(|e| format!("taxonomy setup failed: {}", e))?
        .map(Arc::new);
    let enabled = enable_requested || resolver.is_some();

    let mut taxsum_local: HashMap<String, Option<TaxonomyEvidence>> = HashMap::new();
    if let Some(ref resolver) = resolver {
        for m in metrics {
            let hit_ids: Vec<String> = taxonomy_hits_map
                .get(&m.gene_id)
                .map(|rows| {
                    rows.iter()
                        .take(consensus.top_hits)
                        .map(|r| r.sseqid.clone())
                        .collect()
                })
                .unwrap_or_default();
            let mut evidence = resolver.summarize_panel(&hit_ids, &consensus);
            if evidence.top_hit.is_none() {
                if let Some(acc) = stats.get(&m.gene_id).and_then(|s| s.top_sseqid.clone()) {
                    evidence.top_hit = resolver.lookup(&acc);
                }
            }
            taxsum_local.insert(m.gene_id.clone(), Some(evidence));
        }
        propagate_transcript_taxonomy(&mut taxsum_local, metrics);
    } else {
        for m in metrics {
            taxsum_local.insert(m.gene_id.clone(), None);
        }
    }

    Ok(TaxonomySetup {
        enabled,
        resolver,
        taxsum_map: Arc::new(taxsum_local),
        warnings,
    })
}

fn transcript_root(id: &str) -> Option<String> {
    if let Some((prefix, suffix)) = id.rsplit_once('.') {
        if suffix.len() >= 2
            && suffix.starts_with('t')
            && suffix[1..].chars().all(|c| c.is_ascii_digit())
        {
            return Some(prefix.to_string());
        }
    }
    None
}

pub(crate) fn propagate_transcript_taxonomy(
    map: &mut HashMap<String, Option<TaxonomyEvidence>>,
    metrics: &[GeneMetrics],
) {
    let mut groups: HashMap<String, Vec<String>> = HashMap::new();
    for m in metrics {
        if let Some(root) = transcript_root(&m.gene_id) {
            groups.entry(root).or_default().push(m.gene_id.clone());
        }
    }
    for (_, members) in groups {
        let mut best: Option<TaxonomyEvidence> = None;
        for gene in &members {
            if let Some(Some(ev)) = map.get(gene) {
                if matches!(
                    ev.detail,
                    TaxonomyDetail::Consensus | TaxonomyDetail::CoarseConsensus
                ) && best
                    .as_ref()
                    .is_none_or(|b| ev.congruence_score > b.congruence_score)
                {
                    best = Some(ev.clone());
                }
            }
        }
        let Some(best_ev) = best else { continue };
        for gene in &members {
            let entry = map.entry(gene.clone()).or_insert(None);
            let should_copy = match entry {
                Some(ev) => matches!(
                    ev.detail,
                    TaxonomyDetail::NoHits | TaxonomyDetail::InsufficientHits
                ),
                None => true,
            };
            if should_copy {
                let mut clone = best_ev.clone();
                clone.detail = TaxonomyDetail::Borrowed;
                clone.support = 0;
                clone.support_fraction = 0.0;
                clone.congruence_score = (clone.congruence_score * 0.9).clamp(0.0, 1.0);
                clone.contamination_score = (1.0 - clone.congruence_score).clamp(0.0, 1.0);
                *entry = Some(clone);
            }
        }
    }
}
