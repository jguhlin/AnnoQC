use std::process::Command;

#[derive(Debug, Clone, Default)]
pub struct ToolVersions {
    pub diamond: Option<String>,
    pub mafft: Option<String>,
    pub hmmscan: Option<String>,
}

pub fn check_tool_version(bin: &str, flag: &str) -> Result<String, String> {
    let output = Command::new(bin)
        .arg(flag)
        .output()
        .map_err(|e| format!("failed to run {} {}: {}", bin, flag, e))?;
    if !output.status.success() {
        return Err(format!(
            "{} {} exited with status {}",
            bin, flag, output.status
        ));
    }
    Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
}

pub fn preflight(
    diamond_bin: &str,
    mafft_bin: Option<&str>,
    hmmscan_bin: Option<&str>,
) -> ToolVersions {
    let mut t = ToolVersions::default();
    if let Ok(v) = check_tool_version(diamond_bin, "--version") {
        t.diamond = Some(v);
    }
    if let Some(m) = mafft_bin {
        if let Ok(v) = check_tool_version(m, "--version") {
            t.mafft = Some(v);
        }
    }
    if let Some(h) = hmmscan_bin {
        if let Ok(v) = check_tool_version(h, "-h") {
            t.hmmscan = Some(v);
        }
    }
    t
}
