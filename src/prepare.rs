use std::fs;
use std::fs::{File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};

use crate::*;

pub(crate) fn run_prepare(p: PrepareArgs) -> Result<(), Box<dyn std::error::Error>> {
    let diamond_bin = p.diamond_bin.unwrap_or_else(|| "diamond".to_string());
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
    let prep_json = matches!(p.log_format, LogFormat::Json);
    enforce_taxonomy_metadata(&diamond_bin, db_out, db_done, "diamond makedb");
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

        let mut cmd = std::process::Command::new(&diamond_bin);
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
                &diamond_bin,
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
            &diamond_bin,
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
                &diamond_bin,
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

    let refprot_dir = Path::new("share/uniprot/reference_proteomes");
    let aves_done = Path::new("share/refprot/aves/refprot_aves.done");
    let aves_root = Path::new("share/refprot/aves");
    let aves_db = aves_root.join("aves_refprot.dmnd");
    let aves_fasta_gz = aves_root.join("aves_refprot.fasta.gz");
    let aves_proteome_map = aves_root.join("aves_refprot.proteome_map.tsv");
    enforce_taxonomy_metadata(&diamond_bin, &aves_db, aves_done, "refprot_aves");
    ensure_refprot_integrity(&diamond_bin, &aves_db, &aves_fasta_gz, aves_done);
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

        let missing_taxonomy = aves_db.exists() && !diamond_db_has_taxonomy(&diamond_bin, &aves_db);
        let db_failed_validation =
            aves_db.exists() && !diamond_db_integrity_ok(&diamond_bin, &aves_db);
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
            let mut cmd = std::process::Command::new(&diamond_bin);
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
            let info = run_diamond_dbinfo(&diamond_bin, &tmp_db)
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

fn build_refprot_proteome_map(
    root: &Path,
    map_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = File::create(map_path)?;
    for entry in std::fs::read_dir(root)? {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        if name == "aves_refprot.fasta.gz" || !name.ends_with(".fasta.gz") {
            continue;
        }
        let proteome_id = name.trim_end_matches(".fasta.gz").to_string();
        let mut reader = parse_fastx_file(&path)
            .map_err(|e| format!("refprot map parse {}: {}", path.display(), e))?;
        while let Some(rec) = reader.next() {
            let rec = rec.map_err(|e| format!("refprot map parse {}: {}", path.display(), e))?;
            let acc = String::from_utf8_lossy(rec.id()).to_string();
            if seen.insert(acc.clone()) {
                use std::io::Write as _;
                writeln!(out, "{}\t{}", acc, proteome_id)?;
            }
        }
    }
    Ok(())
}

fn rebuild_combined_gzip(files: &[PathBuf], out_path: &Path) -> Result<(), String> {
    use std::io::copy;
    let out_file = File::create(out_path).map_err(|e| e.to_string())?;
    let mut encoder = GzEncoder::new(out_file, Compression::default());
    for path in files {
        let file = File::open(path).map_err(|e| e.to_string())?;
        let mut decoder = MultiGzDecoder::new(file);
        copy(&mut decoder, &mut encoder)
            .map_err(|e| format!("refprot: concat {} failed: {}", path.display(), e))?;
    }
    encoder.finish().map_err(|e| e.to_string())?;
    Ok(())
}

fn validate_gzip_file(path: &Path) -> Result<(), String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut sink = io::sink();
    io::copy(&mut decoder, &mut sink).map_err(|e| e.to_string())?;
    Ok(())
}

fn diamond_db_integrity_ok(diamond_bin: &str, db_path: &Path) -> bool {
    match run_diamond_dbinfo(diamond_bin, db_path) {
        Ok(_) => true,
        Err(e) => {
            log::warn!("diamond dbinfo failed for {}: {}", db_path.display(), e);
            false
        }
    }
}

fn run_diamond_dbinfo(diamond_bin: &str, db_path: &Path) -> Result<String, String> {
    if !db_path.exists() {
        return Err("db missing".into());
    }
    let output = std::process::Command::new(diamond_bin)
        .arg("dbinfo")
        .arg("--db")
        .arg(db_path)
        .output()
        .map_err(|e| e.to_string())?;
    if !output.status.success() {
        return Err(format!("status {}", output.status));
    }
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

fn dbinfo_extract_hash(text: &str) -> Option<String> {
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("Database hash") {
            return trimmed.split_whitespace().last().map(|s| s.to_string());
        }
    }
    None
}

fn ensure_refprot_integrity(
    diamond_bin: &str,
    db_path: &Path,
    fasta_path: &Path,
    done_marker: &Path,
) {
    if !done_marker.exists() {
        return;
    }
    let mut needs_rerun = false;
    if fasta_path.exists() {
        if let Err(e) = validate_gzip_file(fasta_path) {
            log::warn!(
                "refprot: detected corrupt combined FASTA {} ({}); scheduling rebuild",
                fasta_path.display(),
                e
            );
            let _ = fs::remove_file(fasta_path);
            needs_rerun = true;
        }
    }
    if db_path.exists() && !diamond_db_integrity_ok(diamond_bin, db_path) {
        log::warn!(
            "refprot: DIAMOND database {} failed validation; scheduling rebuild",
            db_path.display()
        );
        let _ = fs::remove_file(db_path);
        needs_rerun = true;
    }
    if needs_rerun {
        log::warn!(
            "refprot: marking {} stale so step will be rerun",
            done_marker.display()
        );
        let _ = fs::remove_file(done_marker);
    }
}

fn http_get_string(url: &str) -> Result<String, String> {
    use curl::easy::Easy;
    let mut data = Vec::new();
    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    let mut transfer = easy.transfer();
    transfer
        .write_function(|new| {
            data.extend_from_slice(new);
            Ok(new.len())
        })
        .map_err(|e| e.to_string())?;
    transfer.perform().map_err(|e| e.to_string())?;
    drop(transfer);
    Ok(String::from_utf8_lossy(&data).to_string())
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct HttpCacheMeta {
    etag: Option<String>,
    last_modified: Option<String>,
}

impl HttpCacheMeta {
    fn is_empty(&self) -> bool {
        self.etag.is_none() && self.last_modified.is_none()
    }
}

enum DownloadStatus {
    NotModified,
    Downloaded(u64),
}

fn load_http_cache_meta(path: &Path) -> Option<HttpCacheMeta> {
    fs::read_to_string(path)
        .ok()
        .and_then(|txt| serde_json::from_str(&txt).ok())
}

fn save_http_cache_meta(path: &Path, meta: &HttpCacheMeta) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let text = serde_json::to_string(meta).map_err(|e| e.to_string())?;
    fs::write(path, text).map_err(|e| e.to_string())
}

fn remote_not_modified(url: &str, meta: &HttpCacheMeta) -> Result<bool, String> {
    use curl::easy::{Easy, List};
    if meta.is_empty() {
        return Ok(false);
    }
    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    easy.nobody(true).map_err(|e| e.to_string())?;
    let mut headers = List::new();
    let mut has = false;
    if let Some(etag) = &meta.etag {
        headers
            .append(&format!("If-None-Match: {}", etag))
            .map_err(|e| e.to_string())?;
        has = true;
    }
    if let Some(lm) = &meta.last_modified {
        headers
            .append(&format!("If-Modified-Since: {}", lm))
            .map_err(|e| e.to_string())?;
        has = true;
    }
    if !has {
        return Ok(false);
    }
    easy.http_headers(headers).map_err(|e| e.to_string())?;
    easy.perform().map_err(|e| e.to_string())?;
    let code = easy.response_code().map_err(|e| e.to_string())?;
    Ok(code == 304)
}

fn http_download(url: &str, dest: &str) -> Result<DownloadStatus, String> {
    use curl::easy::{Easy, List, WriteError};
    use std::io::Write as _;
    let dest_path = Path::new(dest);
    let tmp_path = PathBuf::from(format!("{}.part", dest));
    let meta_path = PathBuf::from(format!("{}.httpmeta", dest));
    let mut meta = load_http_cache_meta(&meta_path).unwrap_or_default();

    if dest_path.exists() && !tmp_path.exists() && remote_not_modified(url, &meta)? {
        return Ok(DownloadStatus::NotModified);
    }

    if tmp_path.exists() && !meta.is_empty() && !remote_not_modified(url, &meta)? {
        let _ = fs::remove_file(&tmp_path);
        meta = HttpCacheMeta::default();
    }

    let resume_from = if tmp_path.exists() {
        fs::metadata(&tmp_path).map(|m| m.len()).unwrap_or(0)
    } else {
        0
    };
    if resume_from == 0 {
        if let Some(parent) = tmp_path.parent() {
            fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        if tmp_path.exists() {
            let _ = fs::remove_file(&tmp_path);
        }
    }
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&tmp_path)
        .map_err(|e| e.to_string())?;

    let mut easy = Easy::new();
    easy.url(url).map_err(|e| e.to_string())?;
    if resume_from > 0 {
        easy.resume_from(resume_from).map_err(|e| e.to_string())?;
    }
    let mut headers = List::new();
    let mut has_headers = false;
    if resume_from == 0 {
        if let Some(etag) = &meta.etag {
            headers
                .append(&format!("If-None-Match: {}", etag))
                .map_err(|e| e.to_string())?;
            has_headers = true;
        }
        if let Some(lm) = &meta.last_modified {
            headers
                .append(&format!("If-Modified-Since: {}", lm))
                .map_err(|e| e.to_string())?;
            has_headers = true;
        }
    } else if let Some(etag) = &meta.etag {
        headers
            .append(&format!("If-Range: {}", etag))
            .map_err(|e| e.to_string())?;
        has_headers = true;
    } else if let Some(lm) = &meta.last_modified {
        headers
            .append(&format!("If-Range: {}", lm))
            .map_err(|e| e.to_string())?;
        has_headers = true;
    }
    if has_headers {
        easy.http_headers(headers).map_err(|e| e.to_string())?;
    }
    let header_meta = Arc::new(Mutex::new(HttpCacheMeta::default()));
    let header_clone = Arc::clone(&header_meta);
    easy.header_function(move |header| {
        if let Ok(text) = std::str::from_utf8(header) {
            let lower = text.to_ascii_lowercase();
            if lower.starts_with("etag:") {
                if let Ok(mut guard) = header_clone.lock() {
                    guard.etag = Some(
                        text.split_once(':')
                            .map(|(_, v)| v.trim().trim_matches('"').to_string())
                            .unwrap_or_default(),
                    );
                }
            } else if lower.starts_with("last-modified:") {
                if let Ok(mut guard) = header_clone.lock() {
                    guard.last_modified = Some(
                        text.split_once(':')
                            .map(|(_, v)| v.trim().to_string())
                            .unwrap_or_default(),
                    );
                }
            }
        }
        true
    })
    .map_err(|e| e.to_string())?;

    let mut written: u64 = 0;
    {
        let mut transfer = easy.transfer();
        transfer
            .write_function(|new| {
                file.write_all(new).map_err(|_| WriteError::Pause)?;
                written += new.len() as u64;
                Ok(new.len())
            })
            .map_err(|e| e.to_string())?;
        if let Err(e) = transfer.perform() {
            let _ = fs::remove_file(&tmp_path);
            return Err(e.to_string());
        }
    }
    let code = easy.response_code().map_err(|e| e.to_string())?;
    if code == 304 {
        let _ = fs::remove_file(&tmp_path);
        return Ok(DownloadStatus::NotModified);
    }
    if code == 416 && dest_path.exists() {
        let _ = fs::remove_file(&tmp_path);
        return Ok(DownloadStatus::NotModified);
    }
    if !(200..300).contains(&code) {
        let _ = fs::remove_file(&tmp_path);
        return Err(format!("download failed with status {}", code));
    }
    file.flush().map_err(|e| e.to_string())?;
    drop(file);
    if dest_path.exists() {
        fs::remove_file(dest_path).map_err(|e| e.to_string())?;
    }
    fs::rename(&tmp_path, dest_path).map_err(|e| e.to_string())?;
    drop(easy);
    let final_meta = Arc::try_unwrap(header_meta)
        .ok()
        .and_then(|m| m.into_inner().ok())
        .unwrap_or_default();
    if !final_meta.is_empty() {
        meta = final_meta;
    }
    save_http_cache_meta(&meta_path, &meta)?;
    let final_size = fs::metadata(dest_path).map(|m| m.len()).unwrap_or(0);
    if final_size == 0 {
        return Err("download produced empty file".into());
    }
    Ok(DownloadStatus::Downloaded(written))
}

fn diamond_db_has_taxonomy(diamond_bin: &str, db: &Path) -> bool {
    if !db.exists() {
        return false;
    }
    if let Ok(text) = run_diamond_dbinfo(diamond_bin, db) {
        if dbinfo_text_has_taxonomy(&text) {
            return true;
        }
    }
    let acc_taxid_candidates = [
        PathBuf::from(format!("{}.acc_taxid.tsv", db.display())),
        db.with_extension("acc_taxid.tsv"),
    ];
    for acc_taxid in acc_taxid_candidates {
        if let Ok(meta) = fs::metadata(&acc_taxid) {
            if meta.len() > 0 {
                return true;
            }
        }
    }
    false
}

pub(crate) fn dbinfo_text_has_taxonomy(text: &str) -> bool {
    let lower = text.to_lowercase();
    lower.contains("taxon") && !lower.contains("no taxonomy")
}

fn enforce_taxonomy_metadata(diamond_bin: &str, db_path: &Path, done_marker: &Path, label: &str) {
    if db_path.exists() && done_marker.exists() && !diamond_db_has_taxonomy(diamond_bin, db_path) {
        log::warn!(
            "{} missing taxonomy metadata; forcing rebuild",
            db_path.display()
        );
        let _ = fs::remove_file(done_marker);
        log::warn!("{} step will rerun to add taxonomy", label);
    }
}
