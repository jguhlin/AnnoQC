use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::consensus_support::PanelProvenanceCounts;
use crate::ecs::GeneMetrics;
use crate::mafft::load_sequences_by_ids;
use crate::metrics;
use crate::taxonomy;

pub fn dump_matches_fasta(
    out_dir: &str,
    ref_fasta: Option<&str>,
    dump_matches_gene: Option<&str>,
    dump_matches_best: bool,
    metrics: &[GeneMetrics],
    panel_map: &HashMap<String, Vec<String>>,
    intrinsic_map: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
) -> Result<(), Box<dyn std::error::Error>> {
    let Some(ref_fasta) = ref_fasta else {
        return Ok(());
    };

    let mut target_gene: Option<String> = dump_matches_gene.map(|s| s.to_string());
    if target_gene.is_none() && dump_matches_best {
        let mut best: Option<(String, usize)> = None;
        for m in metrics {
            let n = panel_map.get(&m.gene_id).map(|v| v.len()).unwrap_or(0);
            if best.as_ref().map(|b| n > b.1).unwrap_or(true) {
                best = Some((m.gene_id.clone(), n));
            }
        }
        if let Some((gid, _)) = best {
            target_gene = Some(gid);
        }
    }

    if let Some(gid) = target_gene {
        let ids = panel_map.get(&gid).cloned().unwrap_or_default();
        if !ids.is_empty() {
            let seqs = load_sequences_by_ids(ref_fasta, &ids).unwrap_or_default();
            let path = Path::new(out_dir).join("matches.fasta");
            let mut f = File::create(path)?;
            if let Some((_im, qseq)) = intrinsic_map.get(&gid) {
                writeln!(f, ">{}", gid)?;
                writeln!(f, "{}", String::from_utf8_lossy(qseq))?;
            }
            for id in &ids {
                let key = taxonomy::canonical_accession(id);
                if let Some(s) = seqs.get(&key) {
                    writeln!(f, ">{}", key)?;
                    writeln!(f, "{}", String::from_utf8_lossy(s))?;
                }
            }
        }
    }

    Ok(())
}

pub fn write_panel_sources_csv(
    out_dir: &str,
    panel_prov_rows: &[(String, PanelProvenanceCounts)],
) -> Result<(), Box<dyn std::error::Error>> {
    if panel_prov_rows.is_empty() {
        return Ok(());
    }

    let prov_path = Path::new(out_dir).join("panel_sources.csv");
    let mut f = File::create(prov_path)?;
    writeln!(f, "gene_id,swissprot_count,refprot_count,cluster_count")?;
    for (gid, prov) in panel_prov_rows {
        writeln!(
            f,
            "{},{},{},{}",
            gid, prov.swissprot, prov.refprot, prov.cluster
        )?;
    }

    Ok(())
}
