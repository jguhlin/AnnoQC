use std::collections::HashMap;

use crate::{refprot, taxonomy};

pub(crate) fn run_taxonomy_count(args: crate::TaxonomyCountArgs) -> Result<(), String> {
    let entries = refprot::parse_readme(&args.readme)
        .map_err(|e| format!("refprot README parse failed: {}", e))?;
    if entries.is_empty() {
        return Err(format!("no proteomes found in {}", args.readme));
    }
    let resolver = taxonomy::TaxonomyResolver::from_sources(None, None, Some(&args.taxdump_dir))
        .map_err(|e| format!("taxonomy resolver setup failed: {}", e))?
        .ok_or_else(|| {
            format!(
                "taxonomy data missing; ensure {} contains nodes.dmp/names.dmp",
                args.taxdump_dir
            )
        })?;
    let target_taxid = if let Some(tid) = args.taxid {
        tid
    } else if let Some(name) = args.name.as_deref() {
        resolver
            .find_taxid_by_name_exact(name)
            .ok_or_else(|| format!("taxonomy name '{}' not found in taxdump", name))?
    } else {
        return Err("taxonomy-count requires --taxid or --name".into());
    };
    let (_, _, target_name) = resolver.reconstruct_lineage_public(target_taxid);
    let selected = refprot::select_by_taxon(&entries, &[target_taxid], &resolver, usize::MAX);
    println!(
        "Taxon: {} ({})",
        target_taxid,
        target_name.unwrap_or_else(|| "unknown".to_string())
    );
    println!("Reference proteomes total: {}", entries.len());
    println!("Proteomes at/under target: {}", selected.len());
    if let Some(rank) = args.rank.as_deref() {
        let buckets = group_proteomes_by_rank(&resolver, &selected, rank);
        if buckets.is_empty() {
            println!(
                "No descendants expose rank '{}' under taxid {}",
                rank, target_taxid
            );
        } else {
            println!("Top {} {} descendants (by proteome count):", args.top, rank);
            for (idx, (taxid, name, count)) in buckets.iter().enumerate() {
                if idx >= args.top {
                    break;
                }
                println!("  {} (taxid {})\t{}", name, taxid, count);
            }
        }
    }
    Ok(())
}

fn group_proteomes_by_rank(
    resolver: &taxonomy::TaxonomyResolver,
    entries: &[refprot::ProteomeEntry],
    rank: &str,
) -> Vec<(u32, String, usize)> {
    let mut counts: HashMap<u32, (usize, String)> = HashMap::new();
    let rank_lower = rank.to_ascii_lowercase();
    for entry in entries {
        let (lineage_names, lineage_ids, _) = resolver.reconstruct_lineage_public(entry.taxid);
        for (tid, name) in lineage_ids.iter().zip(lineage_names.iter()) {
            if let Some(r) = resolver.rank_of(*tid) {
                if r.eq_ignore_ascii_case(&rank_lower) {
                    let entry = counts.entry(*tid).or_insert((0, name.clone()));
                    entry.0 += 1;
                    break;
                }
            }
        }
    }
    let mut out: Vec<(u32, String, usize)> = counts
        .into_iter()
        .map(|(tid, (count, name))| (tid, name, count))
        .collect();
    out.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| a.1.cmp(&b.1)));
    out
}
