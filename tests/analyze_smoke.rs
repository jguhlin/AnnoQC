use std::{
    fs,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    process::Command,
};

use serde_json::Value;
use tempfile::TempDir;

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

#[test]
fn analyze_smoke_produces_outputs() -> Result<(), Box<dyn std::error::Error>> {
    let tmp = TempDir::new()?;
    let stub = compile_diamond_stub(tmp.path());

    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let input = fixtures.join("input.faa");
    let reference = fixtures.join("reference.faa");
    let db = fixtures.join("mock.dmnd");
    let out_dir = tmp.path().join("results");

    let status = Command::new(env!("CARGO_BIN_EXE_AnnoQC"))
        .args([
            "analyze",
            "--fasta",
            input.to_str().unwrap(),
            "--db",
            db.to_str().unwrap(),
            "--diamond-bin",
            stub.to_str().unwrap(),
            "--reference-fasta",
            reference.to_str().unwrap(),
            "--out",
            out_dir.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .status()?;
    assert!(status.success(), "analyze command should succeed");

    let jsonl_path = out_dir.join("qc_report.jsonl");
    let csv_path = out_dir.join("qc_summary.csv");
    assert!(jsonl_path.exists(), "qc_report.jsonl missing");
    assert!(csv_path.exists(), "qc_summary.csv missing");

    let json_file = fs::File::open(&jsonl_path)?;
    let reader = BufReader::new(json_file);
    let mut records: Vec<Value> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let value: Value = serde_json::from_str(&line)?;
        records.push(value);
    }
    assert_eq!(records.len(), 2, "expected two scorecards");
    for (idx, record) in records.iter().enumerate() {
        let gene_id = record["gene_id"].as_str().unwrap();
        let expected = if idx == 0 { "geneA" } else { "geneB" };
        assert_eq!(gene_id, expected);
        let taxonomy_status = record["taxonomy"]["status"].as_str().unwrap();
        assert_eq!(taxonomy_status, "disabled");
        let warnings = record["warnings"].as_array().expect("warnings array");
        assert!(warnings
            .iter()
            .any(|w| w.as_str().unwrap().contains("No DIAMOND hits")));
        assert!(record["score_components"]["taxonomy"].is_null());
    }

    let csv_content = fs::read_to_string(&csv_path)?;
    let mut lines = csv_content.lines();
    let header = lines.next().unwrap();
    let headers: Vec<&str> = header.split(',').collect();
    let taxonomy_score_idx = headers
        .iter()
        .position(|h| *h == "taxonomy_score")
        .expect("taxonomy_score column");
    let taxonomy_status_idx = headers
        .iter()
        .position(|h| *h == "taxonomy_status")
        .expect("taxonomy_status column");
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        assert!(cols.len() >= headers.len(), "unexpected column count");
        assert!(cols[taxonomy_score_idx].is_empty() || cols[taxonomy_score_idx] == "0.0000");
        assert_eq!(cols[taxonomy_status_idx], "disabled");
    }

    Ok(())
}
