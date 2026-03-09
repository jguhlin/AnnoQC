use std::path::Path;

use super::support::enforce_taxonomy_metadata;
use crate::*;

pub(super) fn run_primary_prepare_steps(
    p: &PrepareArgs,
    diamond_bin: &str,
    prep_json: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let fasta_path = Path::new(&p.fasta);
    if !fasta_path.exists() {
        eprintln!(
            "prepare: FASTA {} not found. Download it per book/prepare.md or run scripts/fetch_reference_data.sh",
            fasta_path.display()
        );
        return Ok(());
    }
    let db_out = Path::new(&p.db_out);
    let db_done_buf = format!("{}.done", p.db_out);
    let db_done = Path::new(&db_done_buf);
    enforce_taxonomy_metadata(diamond_bin, db_out, db_done, "diamond makedb");
    checkpoint::run_step(db_done, p.resume, "diamond makedb", prep_json, || {
        if db_out.exists() {
            log::info!(
                "makedb: {} exists; rebuilding due to resume=false",
                db_out.display()
            );
        }
        let map_path = format!("{}.acc_taxid.tsv", &p.db_out);
        match taxonomy::write_cache_from_fasta(&p.fasta, &map_path) {
            Ok(n) => {
                log::info!("makedb: wrote {} accessions to {}", n, map_path);
                let ncbi_map = format!("{}.ncbi.tsv", &p.db_out);
                if let Ok(text) = std::fs::read_to_string(&map_path) {
                    use std::io::Write as _;
                    if let Ok(mut out) = std::fs::File::create(&ncbi_map) {
                        let _ = writeln!(out, "accession.version\ttaxid");
                        for line in text.lines() {
                            let mut it = line.split('\t');
                            if let (Some(acc), Some(tax)) = (it.next(), it.next()) {
                                let _ = writeln!(out, "{}\t{}", acc, tax);
                            }
                        }
                    }
                }
            }
            Err(e) => log::warn!("makedb: unable to build taxon map from FASTA: {}", e),
        }
        let taxdump_dir = std::path::Path::new("share/taxonomy/new_taxdump");
        let has_taxdump =
            taxdump_dir.join("nodes.dmp").exists() && taxdump_dir.join("names.dmp").exists();

        let mut cmd = std::process::Command::new(diamond_bin);
        cmd.arg("makedb")
            .arg("--in")
            .arg(&p.fasta)
            .arg("--db")
            .arg(&p.db_out);
        let ncbi_map = format!("{}.ncbi.tsv", &p.db_out);
        if std::path::Path::new(&ncbi_map).exists() {
            cmd.arg("--taxonmap").arg(&ncbi_map);
        }
        if has_taxdump {
            cmd.arg("--taxonnodes")
                .arg(taxdump_dir.join("nodes.dmp"))
                .arg("--taxonnames")
                .arg(taxdump_dir.join("names.dmp"));
        }
        let status = cmd.status().map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("diamond makedb failed with status {}", status));
        }
        Ok(())
    })?;

    let fasta_for_cluster = if p.fasta.ends_with(".gz") {
        let uncompressed = Path::new(&p.fasta)
            .with_extension("")
            .to_string_lossy()
            .to_string();
        if !Path::new(&uncompressed).exists() {
            log::info!(
                "prepare: decompressing {} -> {} for linclust/cluster",
                p.fasta,
                uncompressed
            );
            let status = std::process::Command::new("sh")
                .arg("-c")
                .arg(format!("gunzip -c '{}' > '{}'", p.fasta, uncompressed))
                .status()
                .map_err(|e| e.to_string())?;
            if !status.success() {
                return Err(format!("gunzip failed for {} with status {}", p.fasta, status).into());
            }
        }
        uncompressed
    } else {
        p.fasta.clone()
    };
    let clusters = Path::new("clusters");
    let clusters_done = Path::new("clusters.done");
    checkpoint::run_step(
        clusters_done,
        p.resume,
        "diamond linclust",
        prep_json,
        || {
            diamond_linclust(
                diamond_bin,
                &fasta_for_cluster,
                clusters,
                p.approx_id,
                p.threads,
            )
            .map_err(|e| format!("linclust failed: {}", e))
        },
    )?;

    let realign = Path::new("clusters.realign");
    let realign_done = Path::new("clusters.realign.done");
    checkpoint::run_step(realign_done, p.resume, "diamond cluster", prep_json, || {
        diamond_cluster(
            diamond_bin,
            &fasta_for_cluster,
            realign,
            p.approx_id,
            p.threads,
        )
        .map_err(|e| format!("cluster failed: {}", e))
    })?;

    let recluster = Path::new("clusters.recluster");
    let recluster_done = Path::new("clusters.recluster.done");
    checkpoint::run_step(
        recluster_done,
        p.resume,
        "diamond recluster",
        prep_json,
        || {
            diamond_recluster(
                diamond_bin,
                &fasta_for_cluster,
                realign,
                recluster,
                p.approx_id,
                p.member_cover,
                p.threads,
            )
            .map_err(|e| format!("recluster failed: {}", e))
        },
    )?;
    if recluster.exists() {
        if let Ok(meta) = recluster.metadata() {
            log::info!(
                "prepare: recluster output {} bytes at {}",
                meta.len(),
                recluster.display()
            );
        }
    } else {
        log::warn!(
            "prepare: recluster output missing at {}",
            recluster.display()
        );
    }

    Ok(())
}
