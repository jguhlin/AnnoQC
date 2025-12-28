use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use needletail::parse_fastx_file;

#[derive(Debug, Clone)]
pub struct TaxonomyResolver {
    accession_to_taxid: HashMap<String, u32>,
    nodes: HashMap<u32, TaxonomyNode>,
    scientific_names: HashMap<u32, String>,
}

#[derive(Debug, Clone)]
pub struct TaxonomyConsensusConfig {
    /// Minimum resolved hits required before consensus is considered reliable.
    pub min_hits: usize,
    /// Maximum number of accessions to consider from the hit list.
    pub top_hits: usize,
    /// Fraction of supporting hits required for consensus (0.0-1.0).
    pub min_support: f64,
    /// Lineage index used for coarse consensus when fine consensus fails.
    pub coarse_rank_index: usize,
    /// Minimum support fraction for coarse consensus (0.0-1.0).
    pub coarse_min_support: f64,
}

impl Default for TaxonomyConsensusConfig {
    fn default() -> Self {
        Self {
            min_hits: 5,
            top_hits: 20,
            min_support: 0.75,
            coarse_rank_index: 1,
            coarse_min_support: 0.6,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaxonomyDetail {
    NoHits,
    InsufficientHits,
    Consensus,
    CoarseConsensus,
    Borrowed,
}

impl fmt::Display for TaxonomyDetail {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            TaxonomyDetail::NoHits => "NoHits",
            TaxonomyDetail::InsufficientHits => "InsufficientHits",
            TaxonomyDetail::Consensus => "Consensus",
            TaxonomyDetail::CoarseConsensus => "CoarseConsensus",
            TaxonomyDetail::Borrowed => "Borrowed",
        };
        write!(f, "{}", text)
    }
}

#[derive(Debug, Clone)]
pub struct TaxonomyEvidence {
    pub top_hit: Option<TaxonSummary>,
    pub consensus: Option<TaxonSummary>,
    pub consensus_rank: Option<String>,
    pub consensus_depth: usize,
    pub support: usize,
    pub considered: usize,
    pub support_fraction: f64,
    pub congruence_score: f64,
    pub contamination_score: f64,
    pub detail: TaxonomyDetail,
}

impl Default for TaxonomyEvidence {
    fn default() -> Self {
        Self {
            top_hit: None,
            consensus: None,
            consensus_rank: None,
            consensus_depth: 0,
            support: 0,
            considered: 0,
            support_fraction: 0.0,
            congruence_score: 0.0,
            contamination_score: 0.0,
            detail: TaxonomyDetail::NoHits,
        }
    }
}

#[derive(Debug, Clone)]
struct FallbackConsensus {
    summary: TaxonSummary,
    support: usize,
    support_fraction: f64,
    depth: usize,
    rank: Option<String>,
    congruence: f64,
}

#[derive(Debug, Clone)]
struct TaxonomyNode {
    parent: u32,
    #[allow(dead_code)]
    rank: String,
}

#[derive(Debug, Clone)]
pub struct TaxonSummary {
    pub taxid: u32,
    pub name: Option<String>,
    pub lineage: Vec<String>,
    pub lineage_ids: Vec<u32>,
}

impl TaxonomyResolver {
    pub fn from_sources(
        cache_path: Option<&str>,
        reference_fasta: Option<&str>,
        taxdump_dir: Option<&str>,
    ) -> Result<Option<Self>, String> {
        let mut nodes = HashMap::new();
        let mut names = HashMap::new();
        if let Some(dir) = taxdump_dir {
            if Path::new(dir).exists() {
                load_taxdump(dir, &mut nodes, &mut names)?;
            } else {
                log::warn!(
                    "Taxdump dir '{}' missing; taxonomy lineage will be limited",
                    dir
                );
            }
        }

        let mut accession_to_taxid = HashMap::new();
        if let Some(cache) = cache_path {
            if Path::new(cache).exists() {
                accession_to_taxid = load_cache(cache)?;
            } else {
                log::warn!("taxonomy cache '{}' not found", cache);
            }
        }

        if accession_to_taxid.is_empty() {
            if let Some(fasta) = reference_fasta {
                if Path::new(fasta).exists() {
                    log::info!("Building taxonomy cache from FASTA {}", fasta);
                    accession_to_taxid = build_cache_from_fasta(fasta)?;
                } else {
                    log::warn!(
                        "reference fasta '{}' missing; taxonomy pillar limited",
                        fasta
                    );
                }
            }
        }

        if accession_to_taxid.is_empty() && nodes.is_empty() {
            return Ok(None);
        }

        Ok(Some(TaxonomyResolver::new(
            accession_to_taxid,
            nodes,
            names,
        )))
    }

    fn new(
        accession_to_taxid: HashMap<String, u32>,
        nodes: HashMap<u32, TaxonomyNode>,
        scientific_names: HashMap<u32, String>,
    ) -> Self {
        Self {
            accession_to_taxid,
            nodes,
            scientific_names,
        }
    }

    pub fn lookup(&self, accession: &str) -> Option<TaxonSummary> {
        let key = canonical_accession(accession);
        let taxid = self.accession_to_taxid.get(&key).copied()?;
        Some(self.summary_for_taxid(taxid))
    }

    pub fn summarize_panel(
        &self,
        accessions: &[String],
        cfg: &TaxonomyConsensusConfig,
    ) -> TaxonomyEvidence {
        let mut evidence = TaxonomyEvidence::default();
        if accessions.is_empty() {
            return evidence;
        }
        let mut resolved: Vec<TaxonSummary> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for acc in accessions.iter().take(cfg.top_hits) {
            let canon = canonical_accession(acc);
            if !seen.insert(canon.clone()) {
                continue;
            }
            if let Some(summary) = self.lookup(acc) {
                if evidence.top_hit.is_none() {
                    evidence.top_hit = Some(summary.clone());
                }
                resolved.push(summary);
            }
        }
        evidence.considered = resolved.len();
        if resolved.is_empty() {
            evidence.detail = TaxonomyDetail::NoHits;
            return evidence;
        }
        if resolved.len() < cfg.min_hits {
            evidence.detail = TaxonomyDetail::InsufficientHits;
        }

        let mut stats: HashMap<u32, (usize, usize)> = HashMap::new();
        let mut max_depth = 0usize;
        for summary in &resolved {
            let depth = summary.lineage_ids.len();
            if depth > max_depth {
                max_depth = depth;
            }
            for (idx, taxid) in summary.lineage_ids.iter().enumerate() {
                stats
                    .entry(*taxid)
                    .and_modify(|entry| {
                        entry.0 += 1;
                        if idx > entry.1 {
                            entry.1 = idx;
                        }
                    })
                    .or_insert((1, idx));
            }
        }
        let total = resolved.len();
        let min_support = cfg.min_support.clamp(0.0, 1.0);
        let quorum = ((min_support * total as f64).ceil() as usize).max(1);
        let mut best: Option<(usize, usize, u32)> = None;
        let mut stats_vec: Vec<(u32, usize, usize)> = stats
            .into_iter()
            .map(|(taxid, (count, depth))| (taxid, count, depth))
            .collect();
        stats_vec.sort_by(|a, b| a.0.cmp(&b.0));
        for (taxid, count, depth) in stats_vec {
            if count < quorum {
                continue;
            }
            match best {
                Some((best_depth, best_count, best_taxid)) => {
                    if depth > best_depth
                        || (depth == best_depth && count > best_count)
                        || (depth == best_depth && count == best_count && taxid < best_taxid)
                    {
                        best = Some((depth, count, taxid));
                    }
                }
                None => best = Some((depth, count, taxid)),
            }
        }
        if let Some((depth, support, taxid)) = best {
            let support_fraction = support as f64 / total as f64;
            let max_depth = max_depth.max(1);
            let depth_norm = (depth as f64 / (max_depth as f64 - 1.0).max(1.0)).clamp(0.0, 1.0);
            let congruence = (support_fraction * depth_norm).clamp(0.0, 1.0);
            let contamination = (1.0 - support_fraction).clamp(0.0, 1.0);
            let consensus = resolved
                .iter()
                .find(|ts| ts.taxid == taxid)
                .cloned()
                .or_else(|| Some(self.summary_for_taxid(taxid)));
            evidence.consensus = consensus;
            evidence.consensus_rank = self
                .nodes
                .get(&taxid)
                .map(|node| node.rank.clone())
                .filter(|s| !s.is_empty());
            evidence.consensus_depth = depth;
            evidence.support = support;
            evidence.support_fraction = support_fraction;
            evidence.congruence_score = congruence;
            evidence.contamination_score = contamination;
            evidence.detail = TaxonomyDetail::Consensus;
        } else if let Some(coarse) = self.fallback_coarse_consensus(&resolved, cfg) {
            evidence.support = coarse.support;
            evidence.considered = total;
            evidence.support_fraction = coarse.support_fraction;
            evidence.congruence_score = coarse.congruence;
            evidence.contamination_score = (1.0 - coarse.support_fraction).clamp(0.0, 1.0);
            evidence.consensus = Some(coarse.summary);
            evidence.consensus_rank = coarse.rank;
            evidence.consensus_depth = coarse.depth;
            evidence.detail = TaxonomyDetail::CoarseConsensus;
        }
        evidence
    }

    fn fallback_coarse_consensus(
        &self,
        resolved: &[TaxonSummary],
        cfg: &TaxonomyConsensusConfig,
    ) -> Option<FallbackConsensus> {
        let idx = cfg.coarse_rank_index;
        let mut counts: HashMap<u32, usize> = HashMap::new();
        let mut label_for: HashMap<u32, String> = HashMap::new();
        for summary in resolved {
            if summary.lineage_ids.len() <= idx {
                continue;
            }
            let taxid = summary.lineage_ids[idx];
            *counts.entry(taxid).or_insert(0) += 1;
            if let Some(name) = summary.lineage.get(idx) {
                label_for.entry(taxid).or_insert_with(|| name.clone());
            }
        }
        if counts.is_empty() {
            return None;
        }
        let total = resolved.len();
        let mut best_tax: Option<(usize, u32)> = None;
        for (taxid, count) in counts {
            if (count as f64 / total as f64) < cfg.coarse_min_support {
                continue;
            }
            match best_tax {
                Some((best_count, best_taxid)) => {
                    if count > best_count || (count == best_count && taxid < best_taxid) {
                        best_tax = Some((count, taxid));
                    }
                }
                None => best_tax = Some((count, taxid)),
            }
        }
        let (support, taxid) = best_tax?;
        let summary = self.summary_for_taxid(taxid);
        let support_fraction = support as f64 / total as f64;
        Some(FallbackConsensus {
            summary,
            support,
            support_fraction,
            depth: idx,
            rank: label_for.get(&taxid).cloned(),
            congruence: (support_fraction * 0.5).clamp(0.0, 1.0),
        })
    }

    fn summary_for_taxid(&self, taxid: u32) -> TaxonSummary {
        let (lineage, lineage_ids, name) = self.reconstruct_lineage(taxid);
        TaxonSummary {
            taxid,
            name,
            lineage,
            lineage_ids,
        }
    }

    fn reconstruct_lineage(&self, taxid: u32) -> (Vec<String>, Vec<u32>, Option<String>) {
        if self.nodes.is_empty() {
            let name = self.scientific_names.get(&taxid).cloned();
            return (Vec::new(), vec![taxid], name);
        }
        let mut lineage = Vec::new();
        let mut lineage_ids = Vec::new();
        let mut current = taxid;
        let mut depth = 0;
        let mut seen = HashSet::new();
        while let Some(node) = self.nodes.get(&current) {
            if let Some(name) = self.scientific_names.get(&current) {
                lineage.push(name.clone());
            }
            lineage_ids.push(current);
            if seen.contains(&current) {
                break;
            }
            seen.insert(current);
            if node.parent == current {
                break;
            }
            current = node.parent;
            depth += 1;
            if depth > 1024 {
                break;
            }
        }
        if lineage_ids.last().copied() != Some(current) {
            lineage_ids.push(current);
        }
        if let Some(root_name) = self.scientific_names.get(&current) {
            if lineage.last().map(|s| s != root_name).unwrap_or(true) {
                lineage.push(root_name.clone());
            }
        }
        let name = self.scientific_names.get(&taxid).cloned();
        lineage.reverse();
        lineage_ids.reverse();
        (lineage, lineage_ids, name)
    }

    pub fn rank_of(&self, taxid: u32) -> Option<&str> {
        self.nodes.get(&taxid).map(|n| n.rank.as_str())
    }

    pub fn parent_of(&self, taxid: u32) -> Option<u32> {
        self.nodes.get(&taxid).map(|n| n.parent)
    }

    // Public helper to expose lineage for external selection code (name/ids).
    pub fn reconstruct_lineage_public(
        &self,
        taxid: u32,
    ) -> (Vec<String>, Vec<u32>, Option<String>) {
        self.reconstruct_lineage(taxid)
    }

    pub fn find_taxid_by_name_exact(&self, name: &str) -> Option<u32> {
        for (tid, n) in &self.scientific_names {
            if n == name {
                return Some(*tid);
            }
        }
        None
    }
}

pub fn canonical_accession(raw: &str) -> String {
    if let Some(idx) = raw.find(' ') {
        return canonical_accession(&raw[..idx]);
    }
    if raw.contains('|') {
        let parts: Vec<&str> = raw.split('|').collect();
        if parts.len() >= 2 && !parts[1].is_empty() {
            return parts[1].to_string();
        }
        if let Some(first) = parts.first() {
            return first.to_string();
        }
    }
    raw.to_string()
}

fn load_cache(path: &str) -> Result<HashMap<String, u32>, String> {
    let file = File::open(path).map_err(|e| format!("unable to open cache {}: {}", path, e))?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {}: {}", path, e))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split('\t');
        if let (Some(acc), Some(tax)) = (parts.next(), parts.next()) {
            if let Ok(taxid) = tax.parse::<u32>() {
                map.insert(acc.to_string(), taxid);
            }
        }
    }
    Ok(map)
}

fn load_taxdump(
    dir: &str,
    nodes: &mut HashMap<u32, TaxonomyNode>,
    scientific: &mut HashMap<u32, String>,
) -> Result<(), String> {
    let nodes_path = Path::new(dir).join("nodes.dmp");
    let names_path = Path::new(dir).join("names.dmp");
    if !nodes_path.exists() {
        return Err(format!("nodes.dmp missing in {}", dir));
    }
    if !names_path.exists() {
        return Err(format!("names.dmp missing in {}", dir));
    }
    let nodes_file = File::open(&nodes_path)
        .map_err(|e| format!("unable to open {}: {}", nodes_path.display(), e))?;
    let reader = BufReader::new(nodes_file);
    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {}: {}", nodes_path.display(), e))?;
        let fields: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if fields.len() < 3 {
            continue;
        }
        let taxid = match fields[0].parse::<u32>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let parent = fields[1].parse::<u32>().unwrap_or(taxid);
        let rank = fields[2].to_string();
        nodes.insert(taxid, TaxonomyNode { parent, rank });
    }

    let names_file = File::open(&names_path)
        .map_err(|e| format!("unable to open {}: {}", names_path.display(), e))?;
    let reader = BufReader::new(names_file);
    for line in reader.lines() {
        let line = line.map_err(|e| format!("error reading {}: {}", names_path.display(), e))?;
        let fields: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if fields.len() < 4 {
            continue;
        }
        if fields[3] == "scientific name" {
            if let Ok(taxid) = fields[0].parse::<u32>() {
                scientific.insert(taxid, fields[1].to_string());
            }
        }
    }

    Ok(())
}

fn build_cache_from_fasta(fasta_path: &str) -> Result<HashMap<String, u32>, String> {
    let mut map = HashMap::new();
    let mut reader = parse_fastx_file(fasta_path).map_err(|e| e.to_string())?;
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| e.to_string())?;
        let defline = String::from_utf8_lossy(record.id());
        let mut tokens = defline.split_whitespace();
        let accession = tokens.next().unwrap_or("");
        if let Some(taxid) = parse_taxid(defline.as_ref()) {
            let key = canonical_accession(accession);
            if !key.is_empty() {
                map.entry(key).or_insert(taxid);
            }
        }
    }
    Ok(map)
}

fn parse_taxid(description: &str) -> Option<u32> {
    for token in description.split_whitespace() {
        if let Some(rest) = token.strip_prefix("OX=") {
            let candidate = rest
                .split(';')
                .next()
                .unwrap_or(rest)
                .trim()
                .trim_end_matches(';');
            if candidate.is_empty() {
                continue;
            }
            if let Ok(val) = candidate.parse::<u32>() {
                return Some(val);
            }
        }
    }
    None
}

#[allow(dead_code)]
pub fn write_cache_from_fasta(fasta_path: &str, out_path: &str) -> Result<usize, String> {
    let map = build_cache_from_fasta(fasta_path)?;
    let mut entries: Vec<_> = map.into_iter().collect();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let mut file = std::fs::File::create(out_path)
        .map_err(|e| format!("unable to create {}: {}", out_path, e))?;
    use std::io::Write;
    for (acc, taxid) in &entries {
        writeln!(file, "{}\t{}", acc, taxid)
            .map_err(|e| format!("error writing {}: {}", out_path, e))?;
    }
    Ok(entries.len())
}

#[allow(dead_code)]
pub fn infer_taxdump_dir(path: Option<&str>) -> Option<PathBuf> {
    path.map(PathBuf::from).filter(|p| p.exists())
}

#[cfg(test)]
mod tests {
    use super::TaxonomyResolver;
    use std::collections::HashMap;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn resolver_from_fasta_and_taxdump() -> Result<(), Box<dyn std::error::Error>> {
        let tmp = TempDir::new()?;
        let fasta = tmp.path().join("ref.faa");
        fs::write(
            &fasta,
            b">sp|P12345|SOME_PROT OX=562; GN=abc\nMKTIIALSYIFCLVFADYKDDD\n",
        )?;
        let taxdump = tmp.path().join("taxdump");
        fs::create_dir_all(&taxdump)?;
        fs::write(
            taxdump.join("nodes.dmp"),
            b"1 | 1 | no rank |\n562 | 1 | species |\n",
        )?;
        fs::write(
            taxdump.join("names.dmp"),
            b"1 | root | | scientific name |\n562 | Escherichia coli | | scientific name |\n",
        )?;

        let resolver = TaxonomyResolver::from_sources(
            None,
            Some(fasta.to_str().unwrap()),
            Some(taxdump.to_str().unwrap()),
        )?
        .expect("resolver constructed");

        let ts = resolver
            .lookup("sp|P12345|SOME_PROT")
            .expect("lookup returns");
        assert_eq!(ts.taxid, 562);
        assert_eq!(ts.name.as_deref(), Some("Escherichia coli"));
        assert_eq!(
            ts.lineage,
            vec!["root".to_string(), "Escherichia coli".to_string()]
        );

        Ok(())
    }

    #[test]
    fn consensus_scores_and_support() {
        use super::{TaxonomyConsensusConfig, TaxonomyDetail, TaxonomyNode};
        let mut accession_to_taxid = HashMap::new();
        accession_to_taxid.insert("A".to_string(), 4);
        accession_to_taxid.insert("B".to_string(), 4);
        accession_to_taxid.insert("C".to_string(), 3);
        accession_to_taxid.insert("D".to_string(), 5);
        accession_to_taxid.insert("E".to_string(), 6);

        let mut nodes = HashMap::new();
        nodes.insert(
            1,
            TaxonomyNode {
                parent: 1,
                rank: "root".into(),
            },
        );
        nodes.insert(
            2,
            TaxonomyNode {
                parent: 1,
                rank: "kingdom".into(),
            },
        );
        nodes.insert(
            3,
            TaxonomyNode {
                parent: 2,
                rank: "genus".into(),
            },
        );
        nodes.insert(
            4,
            TaxonomyNode {
                parent: 3,
                rank: "species".into(),
            },
        );
        nodes.insert(
            5,
            TaxonomyNode {
                parent: 2,
                rank: "genus".into(),
            },
        );
        nodes.insert(
            6,
            TaxonomyNode {
                parent: 2,
                rank: "genus".into(),
            },
        );

        let mut names = HashMap::new();
        names.insert(1, "root".into());
        names.insert(2, "Bacteria".into());
        names.insert(3, "Escherichia".into());
        names.insert(4, "Escherichia coli".into());
        names.insert(5, "Outlierus".into());
        names.insert(6, "Driftus".into());

        let resolver = TaxonomyResolver::new(accession_to_taxid, nodes, names);
        let cfg = TaxonomyConsensusConfig {
            min_hits: 3,
            top_hits: 10,
            min_support: 0.9,
            coarse_rank_index: 1,
            coarse_min_support: 0.6,
        };
        let hits = vec![
            "sp|A|".to_string(),
            "sp|B|".to_string(),
            "tr|C|".to_string(),
            "sp|D|".to_string(),
        ];
        let evidence = resolver.summarize_panel(&hits, &cfg);
        assert!(matches!(
            evidence.detail,
            TaxonomyDetail::Consensus | TaxonomyDetail::CoarseConsensus
        ));
        let consensus = evidence.consensus.expect("consensus taxon present").taxid;
        assert!(consensus == 2 || consensus == 3);
        assert_eq!(evidence.considered, 4);
        assert!(evidence.support >= 2);
        assert!(evidence.support_fraction >= 0.5);
        assert!(evidence.congruence_score > 0.1);
        assert!((evidence.contamination_score - 0.25).abs() < 0.26);

        let scarce_hits = vec!["sp|A|".into(), "sp|B|".into()];
        let scarce = resolver.summarize_panel(&scarce_hits, &cfg);
        assert!(matches!(
            scarce.detail,
            TaxonomyDetail::InsufficientHits
                | TaxonomyDetail::CoarseConsensus
                | TaxonomyDetail::Consensus
        ));

        let cfg_strict = TaxonomyConsensusConfig {
            min_hits: 3,
            top_hits: 10,
            min_support: 0.8,
            coarse_rank_index: 1,
            coarse_min_support: 0.5,
        };
        let mixed_hits = vec!["sp|A|".into(), "sp|D|".into(), "sp|E|".into()];
        let mixed = resolver.summarize_panel(&mixed_hits, &cfg_strict);
        assert_eq!(mixed.considered, 3);
        assert!(matches!(
            mixed.detail,
            TaxonomyDetail::Consensus | TaxonomyDetail::CoarseConsensus
        ));
    }
}
