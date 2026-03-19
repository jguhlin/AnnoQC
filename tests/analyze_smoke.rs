use std::{
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    process::Command,
};

use serde_json::Value;
use tempfile::TempDir;

/// Build a minimal DIAMOND stub binary in the temp dir for smoke tests.
fn compile_diamond_stub(dir: &Path) -> PathBuf {
    let stub_src = dir.join("diamond_stub.rs");
    const STUB_SOURCE: &str = r#"
        use std::io::Read;

        fn main() {
            let mut args = std::env::args();
            let _bin = args.next();
            match args.next().as_deref() {
                Some("--version") => {
                    println!("diamond stub 0.0");
                }
                Some("blastp") => {
                    let mut buffer = Vec::new();
                    let _ = std::io::stdin().read_to_end(&mut buffer);
                    for _ in args {}
                    // no hits emitted; smoke test verifies pipeline handles empty output
                }
                _ => {
                    // consume remaining args to avoid unused variable warnings
                    for _ in args {}
                }
            }
        }
    "#;
    fs::write(&stub_src, STUB_SOURCE).expect("write stub source");
    let stub_bin = dir.join(if cfg!(windows) {
        "diamond_stub.exe"
    } else {
        "diamond_stub"
    });
    let status = Command::new("rustc")
        .arg(&stub_src)
        .arg("-O")
        .arg("-o")
        .arg(&stub_bin)
        .status()
        .expect("compile stub diamond binary");
    assert!(status.success(), "rustc failed to build diamond stub");
    stub_bin
}

fn invalid_data(message: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message)
}

fn load_fasta_ids(
    path: &Path,
) -> Result<std::collections::HashSet<String>, Box<dyn std::error::Error>> {
    let mut ids = std::collections::HashSet::new();
    let content = fs::read_to_string(path)?;
    for line in content.lines() {
        let line = line.trim();
        if !line.starts_with('>') {
            continue;
        }
        let header = line
            .get(1..)
            .ok_or_else(|| invalid_data("fasta header missing id"))?;
        let id = header
            .split_whitespace()
            .next()
            .ok_or_else(|| invalid_data("fasta header missing id"))?;
        if id.is_empty() {
            return Err(invalid_data("fasta header missing id").into());
        }
        ids.insert(id.to_string());
    }
    if ids.is_empty() {
        return Err(invalid_data("no fasta headers found").into());
    }
    Ok(ids)
}

#[test]
fn analyze_smoke_produces_outputs() -> Result<(), Box<dyn std::error::Error>> {
    let tmp = TempDir::new()?;
    let stub = compile_diamond_stub(tmp.path());
    let aligner = std::env::var("SMOKE_ALIGNER").unwrap_or_else(|_| "mafft".to_string());

    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let input = fixtures.join("input.faa");
    let reference = fixtures.join("reference.faa");
    let expected_gene_ids = load_fasta_ids(&input)?;
    let db_src = fixtures.join("mock.dmnd");
    let db = tmp.path().join("mock.dmnd");
    if db_src.exists() {
        let _ = fs::copy(&db_src, &db);
    }
    if !db.exists() || db.metadata().map(|m| m.len() == 0).unwrap_or(true) {
        fs::write(&db, b"stub-db")?;
    }
    let out_dir = tmp.path().join("results");

    let status = Command::new(env!("CARGO_BIN_EXE_AnnoQC"))
        .args(["analyze", "--fasta"])
        .arg(&input)
        .args(["--db"])
        .arg(&db)
        .args(["--diamond-bin"])
        .arg(&stub)
        .args(["--reference-fasta"])
        .arg(&reference)
        .args(["--out"])
        .arg(&out_dir)
        .args(["--threads", "1", "--aligner"])
        .arg(&aligner)
        .status()?;
    assert!(status.success(), "analyze command should succeed");

    let jsonl_path = out_dir.join("qc_report.jsonl");
    let csv_path = out_dir.join("qc_summary.csv");
    assert!(jsonl_path.exists(), "qc_report.jsonl missing");
    assert!(csv_path.exists(), "qc_summary.csv missing");

    let json_file = fs::File::open(&jsonl_path)?;
    let reader = BufReader::new(json_file);
    let mut seen_gene_ids = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let value: Value = serde_json::from_str(&line)?;
        if value["type"].as_str() == Some("metadata") {
            continue;
        }
        let gene_id = value["gene_id"]
            .as_str()
            .ok_or_else(|| invalid_data("missing gene_id"))?;
        assert!(
            expected_gene_ids.contains(gene_id),
            "unexpected gene_id in output: {gene_id}"
        );
        seen_gene_ids.insert(gene_id.to_string());
        let prov = value["panel_provenance"]
            .as_object()
            .ok_or_else(|| invalid_data("missing panel_provenance"))?;
        let swissprot = prov["swissprot"]
            .as_u64()
            .ok_or_else(|| invalid_data("missing swissprot provenance"))?;
        let refprot = prov["refprot"]
            .as_u64()
            .ok_or_else(|| invalid_data("missing refprot provenance"))?;
        let cluster = prov["cluster"]
            .as_u64()
            .ok_or_else(|| invalid_data("missing cluster provenance"))?;
        assert_eq!(swissprot, 0);
        assert_eq!(refprot, 0);
        assert_eq!(cluster, 0);
        let taxonomy_status = value["taxonomy"]["status"]
            .as_str()
            .ok_or_else(|| invalid_data("missing taxonomy status"))?;
        assert_eq!(taxonomy_status, "disabled");
        let warnings = value["warnings"]
            .as_array()
            .ok_or_else(|| invalid_data("warnings not array"))?;
        let has_diamond_warning = warnings.iter().any(|w| {
            w.as_str()
                .map(|msg| msg.contains("No DIAMOND hits"))
                .unwrap_or(false)
        });
        assert!(has_diamond_warning, "missing expected DIAMOND warning");
        assert!(value["score_components"]["taxonomy"].is_null());
    }
    assert_eq!(
        seen_gene_ids, expected_gene_ids,
        "scorecards did not match input genes"
    );

    let csv_content = fs::read_to_string(&csv_path)?;
    let mut lines = csv_content.lines().filter(|line| !line.trim().is_empty());
    let header = loop {
        let Some(line) = lines.next() else {
            panic!("qc_summary.csv missing header");
        };
        if !line.trim_start().starts_with('#') {
            break line;
        }
    };
    let headers: Vec<&str> = header.split(',').collect();
    let panel_sw_idx = headers
        .iter()
        .position(|h| *h == "panel_swissprot")
        .expect("panel_swissprot column");
    let panel_ref_idx = headers
        .iter()
        .position(|h| *h == "panel_refprot")
        .expect("panel_refprot column");
    let panel_cluster_idx = headers
        .iter()
        .position(|h| *h == "panel_cluster")
        .expect("panel_cluster column");
    let taxonomy_score_idx = headers
        .iter()
        .position(|h| *h == "taxonomy_score")
        .expect("taxonomy_score column");
    let taxonomy_status_idx = headers
        .iter()
        .position(|h| *h == "taxonomy_status")
        .expect("taxonomy_status column");
    for line in lines {
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        assert!(cols.len() >= headers.len(), "unexpected column count");
        assert_eq!(cols[panel_sw_idx], "0");
        assert_eq!(cols[panel_ref_idx], "0");
        assert_eq!(cols[panel_cluster_idx], "0");
        assert!(cols[taxonomy_score_idx].is_empty() || cols[taxonomy_score_idx] == "0.0000");
        assert_eq!(cols[taxonomy_status_idx], "disabled");
    }

    Ok(())
}

#[test]
fn analyze_smoke_spoa_backend() -> Result<(), Box<dyn std::error::Error>> {
    let tmp = TempDir::new()?;
    let stub = compile_diamond_stub(tmp.path());

    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let input = fixtures.join("input.faa");
    let reference = fixtures.join("reference.faa");
    let expected_gene_ids = load_fasta_ids(&input)?;
    let db_src = fixtures.join("mock.dmnd");
    let db = tmp.path().join("mock.dmnd");
    if db_src.exists() {
        let _ = fs::copy(&db_src, &db);
    }
    if !db.exists() || db.metadata().map(|m| m.len() == 0).unwrap_or(true) {
        fs::write(&db, b"stub-db")?;
    }
    let out_dir = tmp.path().join("results");

    // Explicitly test SPOA backend
    let status = Command::new(env!("CARGO_BIN_EXE_AnnoQC"))
        .args(["analyze", "--fasta"])
        .arg(&input)
        .args(["--db"])
        .arg(&db)
        .args(["--diamond-bin"])
        .arg(&stub)
        .args(["--reference-fasta"])
        .arg(&reference)
        .args(["--out"])
        .arg(&out_dir)
        .args(["--threads", "1", "--aligner", "spoa"])
        .status()?;
    assert!(status.success(), "analyze command with SPOA should succeed");

    let jsonl_path = out_dir.join("qc_report.jsonl");
    let csv_path = out_dir.join("qc_summary.csv");
    assert!(jsonl_path.exists(), "qc_report.jsonl missing");
    assert!(csv_path.exists(), "qc_summary.csv missing");

    let json_file = fs::File::open(&jsonl_path)?;
    let reader = BufReader::new(json_file);
    let mut seen_gene_ids = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let value: Value = serde_json::from_str(&line)?;
        if value["type"].as_str() == Some("metadata") {
            continue;
        }
        let gene_id = value["gene_id"]
            .as_str()
            .ok_or_else(|| invalid_data("missing gene_id"))?;
        assert!(
            expected_gene_ids.contains(gene_id),
            "unexpected gene_id in output: {gene_id}"
        );
        seen_gene_ids.insert(gene_id.to_string());
    }
    assert_eq!(
        seen_gene_ids, expected_gene_ids,
        "scorecards did not match input genes"
    );

    Ok(())
}
