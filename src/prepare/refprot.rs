use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

use super::support::{
    build_refprot_proteome_map, dbinfo_extract_hash, diamond_db_has_taxonomy,
    diamond_db_integrity_ok, enforce_taxonomy_metadata, ensure_refprot_integrity, http_download,
    http_get_string, rebuild_combined_gzip, run_diamond_dbinfo, validate_gzip_file, DownloadStatus,
};
use crate::*;

pub(super) fn run_aves_refprot_step(
    p: &PrepareArgs,
    diamond_bin: &str,
    prep_json: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let refprot_dir = Path::new("share/uniprot/reference_proteomes");
    let aves_done = Path::new("share/refprot/aves/refprot_aves.done");
    let aves_root = Path::new("share/refprot/aves");
    let aves_db = aves_root.join("aves_refprot.dmnd");
    let aves_fasta_gz = aves_root.join("aves_refprot.fasta.gz");
    let aves_proteome_map = aves_root.join("aves_refprot.proteome_map.tsv");
    enforce_taxonomy_metadata(diamond_bin, &aves_db, aves_done, "refprot_aves");
    ensure_refprot_integrity(diamond_bin, &aves_db, &aves_fasta_gz, aves_done);
    checkpoint::run_step(aves_done, p.resume, "refprot_aves", prep_json, || {
        let readme_path = refprot_dir.join("README");
        let mut readme_refreshed = false;
        if readme_path.exists() {
            if let Ok(meta) = std::fs::metadata(&readme_path) {
                if let Ok(modified) = meta.modified() {
                    if let Ok(age) = modified.elapsed() {
                        if age.as_secs() > 14 * 24 * 60 * 60 {
                            let url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README";
                            log::info!(
                                "refprot: README older than 14 days; refreshing from {}",
                                url
                            );
                            let dest_path = readme_path.to_string_lossy().to_string();
                            match http_download(url, &dest_path) {
                                Ok(DownloadStatus::Downloaded(_)) => {
                                    readme_refreshed = true;
                                }
                                Ok(DownloadStatus::NotModified) => {}
                                Err(e) => {
                                    log::warn!("refprot: README refresh failed: {}", e);
                                }
                            }
                        }
                    }
                }
            }
        }
        if !readme_path.exists() {
            log::info!("refprot: README not found; skipping Aves build");
            return Ok(());
        }
        std::fs::create_dir_all(aves_root).map_err(|e| e.to_string())?;
        let mut cached_count = 0usize;
        for entry in std::fs::read_dir(aves_root).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if !path.is_file() {
                continue;
            }
            let name = path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or_default();
            if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
                continue;
            }
            cached_count += 1;
        }
        let mut new_count = 0usize;
        let need_download = cached_count == 0 || readme_refreshed;
        if need_download {
            let resolver_opt = taxonomy::TaxonomyResolver::from_sources(
                None,
                None,
                Some("share/taxonomy/new_taxdump"),
            )
            .map_err(|e| e.to_string())?;
            let resolver = if let Some(r) = resolver_opt {
                r
            } else {
                return Ok(());
            };
            let entries =
                refprot::parse_readme(readme_path.to_str().unwrap()).map_err(|e| e.to_string())?;
            let aves_list = refprot::select_by_taxon(&entries, &[8782u32], &resolver, usize::MAX);
            log::info!(
                "refprot: selected {} Aves proteomes for download",
                aves_list.len()
            );
            let base = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes";
            let (tx, rx) = std::sync::mpsc::channel::<(std::path::PathBuf, bool)>();
            let jobs = aves_list
                .into_iter()
                .map(|e| (e.proteome_id, e.division, base.to_string()))
                .collect::<Vec<_>>();
            let mut handles = Vec::new();
            let workers = p.refprot_workers.max(1);
            let chunk_size = jobs.len().div_ceil(workers);
            let max_retries = p.download_retries;
            for chunk in jobs.chunks(chunk_size.max(1)) {
                let chunk = chunk.to_vec();
                let txc = tx.clone();
                let root = aves_root.to_path_buf();
                handles.push(std::thread::spawn(move || {
                    for (pid, division, base_url) in chunk {
                        let div = {
                            let mut d = division.clone();
                            if d.is_empty() {
                                d = "eukaryota".into();
                            }
                            let mut ch = d.chars();
                            match ch.next() {
                                Some(c) => format!("{}{}", c.to_ascii_uppercase(), ch.as_str()),
                                None => "Eukaryota".into(),
                            }
                        };
                        let url_dir = format!("{}/{}/{}/", base_url, div, pid);
                        let html = http_get_string(&url_dir).unwrap_or_default();
                        let mut best: Option<String> = None;
                        for tok in html.split(|c: char| c == '"' || c.is_whitespace()) {
                            if tok.starts_with(&pid)
                                && tok.ends_with(".fasta.gz")
                                && !tok.contains("_DNA")
                            {
                                if !tok.contains("additional") {
                                    best = Some(tok.to_string());
                                    break;
                                }
                                if best.is_none() {
                                    best = Some(tok.to_string());
                                }
                            }
                        }
                        if let Some(fname) = best {
                            let dest = root.join(format!("{}.fasta.gz", pid));
                            let dest_string = dest.to_string_lossy().to_string();
                            let file_url = format!("{}{}", url_dir, fname);
                            let mut downloaded_now = false;
                            let mut attempts = 0usize;
                            loop {
                                match http_download(&file_url, &dest_string) {
                                    Ok(DownloadStatus::Downloaded(bytes)) => {
                                        log::info!("refprot: fetched {} ({} bytes)", pid, bytes);
                                        downloaded_now = true;
                                        break;
                                    }
                                    Ok(DownloadStatus::NotModified) => break,
                                    Err(err) => {
                                        if attempts >= max_retries {
                                            log::warn!(
                                                "refprot: failed to fetch {} after {} retries: {}",
                                                pid,
                                                max_retries,
                                                err
                                            );
                                            break;
                                        }
                                        let backoff =
                                            Duration::from_secs(2u64.pow(attempts.min(4) as u32));
                                        log::warn!(
                                            "refprot: retry {} for {} after {}s ({})",
                                            attempts + 1,
                                            pid,
                                            backoff.as_secs(),
                                            err
                                        );
                                        std::thread::sleep(backoff);
                                        attempts += 1;
                                    }
                                }
                            }
                            if dest.exists() {
                                let _ = txc.send((dest.clone(), downloaded_now));
                            }
                        }
                    }
                }));
            }
            drop(tx);
            for h in handles {
                let _ = h.join();
            }
            let mut downloaded: Vec<(std::path::PathBuf, bool)> = Vec::new();
            while let Ok(p) = rx.recv_timeout(std::time::Duration::from_millis(10)) {
                downloaded.push(p);
            }
            new_count = downloaded.iter().filter(|(_, fresh)| *fresh).count();
        } else {
            log::info!(
                "refprot: cached {} proteomes present; skipping download",
                cached_count
            );
        }
        let mut concat_invalid = false;
        if aves_fasta_gz.exists() {
            if let Err(e) = validate_gzip_file(&aves_fasta_gz) {
                log::warn!(
                    "refprot: combined FASTA {} failed validation ({}); rebuilding",
                    aves_fasta_gz.display(),
                    e
                );
                let _ = fs::remove_file(&aves_fasta_gz);
                concat_invalid = true;
            }
        }
        let need_concat = concat_invalid || !aves_fasta_gz.exists() || new_count > 0;
        if need_concat {
            log::info!(
                "refprot: rebuilding combined FASTA ({} new proteomes)",
                new_count
            );
            let mut files: Vec<PathBuf> = Vec::new();
            for entry in std::fs::read_dir(aves_root).map_err(|e| e.to_string())? {
                let path = entry.map_err(|e| e.to_string())?.path();
                if !path.is_file() {
                    continue;
                }
                let name = path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or_default();
                if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
                    continue;
                }
                files.push(path);
            }
            files.sort();
            if files.is_empty() {
                log::warn!(
                    "refprot: no proteome FASTAs available under {}",
                    aves_root.display()
                );
            } else {
                rebuild_combined_gzip(&files, &aves_fasta_gz)?;
            }
        }

        if !aves_proteome_map.exists() || new_count > 0 {
            build_refprot_proteome_map(aves_root, &aves_proteome_map)
                .map_err(|e| format!("refprot: proteome map build failed: {}", e))?;
        }

        let missing_taxonomy = aves_db.exists() && !diamond_db_has_taxonomy(diamond_bin, &aves_db);
        let db_failed_validation =
            aves_db.exists() && !diamond_db_integrity_ok(diamond_bin, &aves_db);
        let need_makedb = missing_taxonomy
            || db_failed_validation
            || !aves_db.exists()
            || new_count > 0
            || !aves_fasta_gz.exists();
        if need_makedb {
            let aves_acc = aves_root.join("aves_refprot.acc_taxid.tsv");
            let aves_ncbi = aves_root.join("aves_refprot.ncbi.tsv");
            match taxonomy::write_cache_from_fasta(
                &aves_fasta_gz.to_string_lossy(),
                &aves_acc.to_string_lossy(),
            ) {
                Ok(n) => {
                    log::info!("refprot: wrote {} accessions to {}", n, aves_acc.display());
                    if let Ok(text) = std::fs::read_to_string(&aves_acc) {
                        if let Ok(mut f) = std::fs::File::create(&aves_ncbi) {
                            use std::io::Write as _;
                            let _ = writeln!(f, "accession.version\ttaxid");
                            for line in text.lines() {
                                let mut it = line.split('\t');
                                if let (Some(acc), Some(tid)) = (it.next(), it.next()) {
                                    let _ = writeln!(f, "{}\t{}", acc, tid);
                                }
                            }
                        }
                    }
                }
                Err(e) => log::warn!("refprot: unable to build taxon map: {}", e),
            }
            let taxdump_dir = std::path::Path::new("share/taxonomy/new_taxdump");
            let has_taxdump =
                taxdump_dir.join("nodes.dmp").exists() && taxdump_dir.join("names.dmp").exists();
            let tmp_db = aves_root.join("aves_refprot.tmp.dmnd");
            if tmp_db.exists() {
                let _ = fs::remove_file(&tmp_db);
            }
            let mut cmd = std::process::Command::new(diamond_bin);
            cmd.arg("makedb")
                .arg("--in")
                .arg(&aves_fasta_gz)
                .arg("--db")
                .arg(&tmp_db);
            if aves_ncbi.exists() {
                cmd.arg("--taxonmap").arg(&aves_ncbi);
            }
            if has_taxdump {
                cmd.arg("--taxonnodes")
                    .arg(taxdump_dir.join("nodes.dmp"))
                    .arg("--taxonnames")
                    .arg(taxdump_dir.join("names.dmp"));
            }
            let status = cmd.status().map_err(|e| e.to_string())?;
            if !status.success() {
                let _ = fs::remove_file(&tmp_db);
                return Err(format!(
                    "diamond makedb (refprot) failed with status {}",
                    status
                ));
            }
            let info = run_diamond_dbinfo(diamond_bin, &tmp_db)
                .map_err(|e| format!("refprot: dbinfo failed after makedb: {}", e))?;
            let hash = dbinfo_extract_hash(&info);
            if aves_db.exists() {
                let _ = fs::remove_file(&aves_db);
            }
            fs::rename(&tmp_db, &aves_db)
                .map_err(|e| format!("refprot: rename tmp db failed: {}", e))?;
            if let Some(hash) = hash {
                let hash_path = aves_root.join("aves_refprot.hash");
                let _ = fs::write(hash_path, format!("{}\n", hash));
            }
        }
        Ok(())
    })?;
    Ok(())
}
