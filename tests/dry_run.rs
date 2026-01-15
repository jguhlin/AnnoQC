use std::path::Path;
use std::process::Command;

#[test]
fn analyze_dry_run_prints_rubric() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let input = fixtures.join("input.faa");
    let db = fixtures.join("mock.dmnd");

    let output = Command::new(env!("CARGO_BIN_EXE_AnnoQC"))
        .args([
            "analyze",
            "--fasta",
            input.to_str().unwrap(),
            "--db",
            db.to_str().unwrap(),
            "--dry-run",
        ])
        .output()
        .expect("run analyze --dry-run");

    assert!(output.status.success(), "dry-run should succeed");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("scoring_rubric"));
    assert!(stdout.contains("weights_raw:"));
}
