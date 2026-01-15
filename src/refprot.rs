use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteomeEntry {
    pub proteome_id: String, // e.g. UP000005640
    pub taxid: u32,
    pub organism: String,
    pub division: String, // eukaryota/bacteria/viruses/archaea
}

/// Parse UniProt Reference Proteomes README to extract proteome_id → (taxid, organism).
/// The README format may change; we use regexes tolerant to whitespace.
pub fn parse_readme(path: &str) -> Result<Vec<ProteomeEntry>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("unable to read {}: {}", path, e))?;
    let mut entries = Vec::new();
    for line in text.lines() {
        // Format observed: UP000216524\t463024\t...\torganism name
        let line = line.trim_start();
        if !line.starts_with("UP") {
            continue;
        }
        let fields: Vec<&str> = if line.contains('\t') {
            line.split('\t')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .collect()
        } else {
            line.split_whitespace().collect()
        };
        let pid = match fields.first() {
            Some(s) if s.len() >= 11 && s.starts_with("UP") => (*s).to_string(),
            _ => continue,
        };
        let taxid = match fields.get(1).and_then(|s| s.parse::<u32>().ok()) {
            Some(t) => t,
            None => continue,
        };
        let division = match fields.get(3).map(|s| s.to_ascii_lowercase()) {
            Some(div)
                if matches!(
                    div.as_str(),
                    "bacteria" | "archaea" | "eukaryota" | "viruses"
                ) =>
            {
                div
            }
            _ => "unknown".to_string(),
        };
        let org = if fields.len() >= 8 {
            fields[7..].join(" ")
        } else {
            fields.last().copied().unwrap_or_default().to_string()
        };
        entries.push(ProteomeEntry {
            proteome_id: pid,
            taxid,
            organism: org,
            division,
        });
    }
    Ok(entries)
}

/// Select proteomes whose taxid is a descendant of any target taxid according to the resolver hierarchy.
/// Stops after collecting `max_proteomes` entries.
pub fn select_by_taxon(
    entries: &[ProteomeEntry],
    targets: &[u32],
    resolver: &crate::taxonomy::TaxonomyResolver,
    max_proteomes: usize,
) -> Vec<ProteomeEntry> {
    let target_set: HashSet<u32> = targets.iter().copied().collect();
    let mut memo: HashMap<u32, bool> = HashMap::new();
    let mut out = Vec::new();
    'outer: for e in entries {
        if out.len() >= max_proteomes {
            break;
        }
        // climb parents until root; if any equals a target, accept
        let mut cur = e.taxid;
        let mut chain: Vec<u32> = Vec::new();
        let mut hops = 0usize;
        let result: Option<bool>;
        loop {
            if let Some(cached) = memo.get(&cur).copied() {
                result = Some(cached);
                break;
            }
            if target_set.contains(&cur) {
                result = Some(true);
                break;
            }
            chain.push(cur);
            let parent = match resolver.parent_of(cur) {
                Some(parent) => parent,
                None => {
                    result = Some(false);
                    break;
                }
            };
            if parent == cur || chain.contains(&parent) {
                result = Some(false);
                break;
            }
            cur = parent;
            hops += 1;
            if hops > 2000 {
                result = Some(false);
                break;
            }
        }
        let accepted = result.unwrap_or(false);
        for taxid in chain {
            memo.insert(taxid, accepted);
        }
        if accepted {
            memo.entry(cur).or_insert(true);
            out.push(e.clone());
            continue 'outer;
        }
    }
    out
}
