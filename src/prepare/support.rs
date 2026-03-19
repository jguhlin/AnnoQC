use std::fs;
use std::fs::{File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};

pub(crate) fn build_refprot_proteome_map(
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

pub(crate) fn rebuild_combined_gzip(files: &[PathBuf], out_path: &Path) -> Result<(), String> {
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

pub(crate) fn validate_gzip_file(path: &Path) -> Result<(), String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut sink = io::sink();
    io::copy(&mut decoder, &mut sink).map_err(|e| e.to_string())?;
    Ok(())
}

pub(crate) fn diamond_db_integrity_ok(diamond_bin: &str, db_path: &Path) -> bool {
    match run_diamond_dbinfo(diamond_bin, db_path) {
        Ok(_) => true,
        Err(e) => {
            log::warn!("diamond dbinfo failed for {}: {}", db_path.display(), e);
            false
        }
    }
}

pub(crate) fn run_diamond_dbinfo(diamond_bin: &str, db_path: &Path) -> Result<String, String> {
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

pub(crate) fn dbinfo_extract_hash(text: &str) -> Option<String> {
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("Database hash") {
            return trimmed.split_whitespace().last().map(|s| s.to_string());
        }
    }
    None
}

pub(crate) fn ensure_refprot_integrity(
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

pub(crate) fn http_get_string(url: &str) -> Result<String, String> {
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

pub(crate) enum DownloadStatus {
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

pub(crate) fn http_download(url: &str, dest: &str) -> Result<DownloadStatus, String> {
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
    let header_meta = Arc::new(Mutex::new(HttpCacheMeta::default()));
    let header_meta_clone = Arc::clone(&header_meta);
    {
        let mut headers = List::new();
        if let Some(etag) = &meta.etag {
            headers
                .append(&format!("If-None-Match: {}", etag))
                .map_err(|e| e.to_string())?;
        }
        if let Some(lm) = &meta.last_modified {
            headers
                .append(&format!("If-Modified-Since: {}", lm))
                .map_err(|e| e.to_string())?;
        }
        if !meta.is_empty() {
            easy.http_headers(headers).map_err(|e| e.to_string())?;
        }
    }
    let mut transfer = easy.transfer();
    transfer
        .header_function(move |header| {
            let line = String::from_utf8_lossy(header).trim().to_string();
            let lower = line.to_ascii_lowercase();
            if let Some(value) = line.split_once(':').map(|(_, v)| v.trim().to_string()) {
                if lower.starts_with("etag:") {
                    if let Ok(mut meta) = header_meta_clone.lock() {
                        meta.etag = Some(value);
                    }
                } else if lower.starts_with("last-modified:") {
                    if let Ok(mut meta) = header_meta_clone.lock() {
                        meta.last_modified = Some(value);
                    }
                }
            }
            true
        })
        .map_err(|e| e.to_string())?;
    transfer
        .write_function(|data| match file.write_all(data) {
            Ok(()) => Ok(data.len()),
            Err(_) => Err(WriteError::Pause),
        })
        .map_err(|e| e.to_string())?;
    transfer.perform().map_err(|e| e.to_string())?;
    drop(transfer);

    let response_code = easy.response_code().map_err(|e| e.to_string())?;
    if response_code == 304 {
        let _ = fs::remove_file(&tmp_path);
        return Ok(DownloadStatus::NotModified);
    }
    if response_code >= 400 {
        return Err(format!("HTTP {}", response_code));
    }

    fs::rename(&tmp_path, dest_path).map_err(|e| e.to_string())?;
    if let Ok(meta_guard) = header_meta.lock() {
        let fresh = meta_guard.clone();
        if !fresh.is_empty() {
            save_http_cache_meta(&meta_path, &fresh)?;
        }
    }
    let bytes = fs::metadata(dest_path).map(|m| m.len()).unwrap_or(0);
    Ok(DownloadStatus::Downloaded(bytes))
}

pub(crate) fn diamond_db_has_taxonomy(diamond_bin: &str, db_path: &Path) -> bool {
    match run_diamond_dbinfo(diamond_bin, db_path) {
        Ok(text) => {
            if dbinfo_text_has_taxonomy(&text) {
                true
            } else {
                log::warn!(
                    "diamond db {} missing taxonomy metadata; scheduling rebuild",
                    db_path.display()
                );
                false
            }
        }
        Err(err) => {
            log::warn!(
                "diamond dbinfo unavailable for {}: {}; assuming missing taxonomy",
                db_path.display(),
                err
            );
            false
        }
    }
}

pub(crate) fn dbinfo_text_has_taxonomy(text: &str) -> bool {
    text.lines().any(|line| {
        let lower = line.to_ascii_lowercase();
        if lower.contains("no taxonomy") {
            return false;
        }
        lower.contains("taxon count") || lower.contains("taxonomy")
    })
}

pub(crate) fn enforce_taxonomy_metadata(
    diamond_bin: &str,
    db_path: &Path,
    done_marker: &Path,
    label: &str,
) {
    if !done_marker.exists() || !db_path.exists() {
        return;
    }
    if diamond_db_has_taxonomy(diamond_bin, db_path) {
        return;
    }
    log::warn!(
        "{}: removing stale marker {} so taxonomy-enabled rebuild can run",
        label,
        done_marker.display()
    );
    let _ = fs::remove_file(done_marker);
}
