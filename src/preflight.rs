use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone, Default)]
/// Resolved tool versions discovered during preflight checks.
pub struct ToolVersions {
    /// DIAMOND version string when available.
    pub diamond: Option<String>,
    /// MAFFT version string when available.
    pub mafft: Option<String>,
    /// HMMER hmmscan version string when available.
    pub hmmscan: Option<String>,
}

/// Resolve a binary name or path to an executable path.
fn resolve_executable(bin: &str) -> Result<PathBuf, String> {
    let bin = bin.trim();
    if bin.is_empty() {
        return Err("binary path is empty".to_string());
    }
    let bin_path = Path::new(bin);
    if bin_path.components().count() > 1 {
        return validate_executable_path(bin_path);
    }
    let path_var = std::env::var_os("PATH").ok_or_else(|| "PATH is not set".to_string())?;
    for dir in std::env::split_paths(&path_var) {
        let candidate = dir.join(bin);
        if let Ok(path) = validate_executable_path(&candidate) {
            return Ok(path);
        }
    }
    Err(format!("binary '{}' not found in PATH", bin))
}

fn validate_executable_path(path: &Path) -> Result<PathBuf, String> {
    let meta = path
        .metadata()
        .map_err(|e| format!("binary '{}' is not accessible: {}", path.display(), e))?;
    if !meta.is_file() {
        return Err(format!("binary '{}' is not a file", path.display()));
    }
    if !is_executable(&meta) {
        return Err(format!("binary '{}' is not executable", path.display()));
    }
    Ok(path.to_path_buf())
}

fn is_executable(meta: &std::fs::Metadata) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        meta.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        let _ = meta;
        true
    }
}

fn output_preview(bytes: &[u8], max_len: usize) -> String {
    let mut text = String::from_utf8_lossy(bytes).trim().to_string();
    if text.len() > max_len {
        text.truncate(max_len);
    }
    text
}

/// Run a tool with a version/help flag and return its output (stdout/stderr).
pub fn check_tool_version(bin: &str, flag: &str) -> Result<String, String> {
    let bin_path = resolve_executable(bin)?;
    let output = Command::new(&bin_path)
        .arg(flag)
        .output()
        .map_err(|e| format!("failed to run {} {}: {}", bin_path.display(), flag, e))?;
    if !output.status.success() {
        return Err(format!(
            "{} {} exited with status {} stderr='{}'",
            bin_path.display(),
            flag,
            output.status,
            output_preview(&output.stderr, 400)
        ));
    }
    let stdout = output_preview(&output.stdout, 400);
    if !stdout.is_empty() {
        return Ok(stdout);
    }
    let stderr = output_preview(&output.stderr, 400);
    if !stderr.is_empty() {
        return Ok(stderr);
    }
    Err(format!(
        "{} {} produced no output",
        bin_path.display(),
        flag
    ))
}

/// Run preflight checks and return discovered tool versions.
pub fn preflight(
    diamond_bin: &str,
    mafft_bin: Option<&str>,
    hmmscan_bin: Option<&str>,
) -> ToolVersions {
    let mut t = ToolVersions::default();
    match check_tool_version(diamond_bin, "--version") {
        Ok(v) => t.diamond = Some(v),
        Err(e) => log::warn!("preflight: diamond version check failed: {}", e),
    }
    if let Some(m) = mafft_bin {
        match check_tool_version(m, "--version") {
            Ok(v) => t.mafft = Some(v),
            Err(e) => log::warn!("preflight: mafft version check failed: {}", e),
        }
    }
    if let Some(h) = hmmscan_bin {
        match check_tool_version(h, "-h") {
            Ok(v) => t.hmmscan = Some(v),
            Err(e) => log::warn!("preflight: hmmscan version check failed: {}", e),
        }
    }
    t
}

/// Validate that the DIAMOND database exists and dbinfo can read it.
pub fn check_diamond_db(diamond_bin: &str, db: &str) -> Result<(), String> {
    let path = Path::new(db);
    if !path.exists() {
        return Err(format!("DIAMOND database not found: {}", path.display()));
    }
    if let Ok(meta) = path.metadata() {
        if meta.len() == 0 {
            return Err(format!("DIAMOND database is empty: {}", path.display()));
        }
    }
    let diamond_path = resolve_executable(diamond_bin)?;
    let output = Command::new(&diamond_path)
        .arg("dbinfo")
        .arg("--db")
        .arg(db)
        .output()
        .map_err(|e| format!("failed to run {} dbinfo: {}", diamond_path.display(), e))?;
    if !output.status.success() {
        return Err(format!(
            "DIAMOND dbinfo failed: status={} stderr='{}'",
            output.status,
            output_preview(&output.stderr, 400)
        ));
    }
    Ok(())
}
