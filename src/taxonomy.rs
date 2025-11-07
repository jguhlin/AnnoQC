use std::collections::HashMap;
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

        if accession_to_taxid.is_empty() {
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
        let (lineage, name) = self.reconstruct_lineage(taxid);
        Some(TaxonSummary {
            taxid,
            name,
            lineage,
        })
    }

    fn reconstruct_lineage(&self, taxid: u32) -> (Vec<String>, Option<String>) {
        if self.nodes.is_empty() {
            let name = self.scientific_names.get(&taxid).cloned();
            return (Vec::new(), name);
        }
        let mut lineage = Vec::new();
        let mut current = taxid;
        let mut depth = 0;
        let mut seen = std::collections::HashSet::new();
        while let Some(node) = self.nodes.get(&current) {
            if let Some(name) = self.scientific_names.get(&current) {
                lineage.push(name.clone());
            }
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
        if let Some(root_name) = self.scientific_names.get(&current) {
            if lineage.last().map(|s| s != root_name).unwrap_or(true) {
                lineage.push(root_name.clone());
            }
        }
        let name = self.scientific_names.get(&taxid).cloned();
        lineage.reverse();
        (lineage, name)
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
            if let Ok(val) = rest.trim_end_matches(';').parse::<u32>() {
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
}
