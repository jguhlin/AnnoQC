use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(Debug)]
struct MermaidBlock {
    source_path: String,
    fence_line: usize,         // 1-based
    content_start_line: usize, // 1-based
    content_end_line: usize,   // 1-based
    text: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Backend {
    Auto,
    Mmdc,
    Npx,
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    let mut root = "book/src".to_string();
    let mut backend = Backend::Auto;
    let mut mmdc_cmd = "mmdc".to_string();
    let mut npx_pkg = "@mermaid-js/mermaid-cli@9.2.2".to_string();
    let mut verbose = false;
    let mut fail_fast = false;
    let mut keep_temp = false;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--root" => {
                i += 1;
                root = args.get(i).ok_or("missing value for --root")?.clone();
            }
            "--backend" => {
                i += 1;
                let value = args.get(i).ok_or("missing value for --backend")?;
                backend = match value.as_str() {
                    "auto" => Backend::Auto,
                    "mmdc" => Backend::Mmdc,
                    "npx" => Backend::Npx,
                    other => {
                        return Err(format!(
                            "invalid --backend value: {other} (expected: auto|mmdc|npx)"
                        ));
                    }
                };
            }
            "--mmdc" => {
                i += 1;
                mmdc_cmd = args.get(i).ok_or("missing value for --mmdc")?.clone();
            }
            "--npx-pkg" => {
                i += 1;
                npx_pkg = args.get(i).ok_or("missing value for --npx-pkg")?.clone();
            }
            "--verbose" => verbose = true,
            "--fail-fast" => fail_fast = true,
            "--keep-temp" => keep_temp = true,
            "-h" | "--help" => {
                print_help();
                return Ok(());
            }
            other if other.starts_with('-') => {
                return Err(format!("unknown flag: {other}"));
            }
            value => {
                // Back-compat positional root
                root = value.to_string();
            }
        }
        i += 1;
    }

    let root_path = Path::new(&root);
    if !root_path.exists() {
        return Err(format!("root path not found: {root}"));
    }

    let mut md_files = Vec::new();
    collect_md_files(root_path, &mut md_files).map_err(|e| e.to_string())?;
    md_files.sort();

    let mut blocks = Vec::new();
    for path in &md_files {
        let text = fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        blocks.extend(extract_mermaid_blocks(&text, &path.to_string_lossy()));
    }

    if blocks.is_empty() {
        if verbose {
            println!("no mermaid blocks found under {}", root_path.display());
        }
        return Ok(());
    }

    validate_blocks(
        &blocks, backend, &mmdc_cmd, &npx_pkg, keep_temp, verbose, fail_fast,
    )
}

fn print_help() {
    println!(
        "\
check_mermaid: validate Mermaid fenced blocks outside mdBook (via Mermaid CLI)

Usage:
  cargo run --bin check_mermaid -- [--root book/src]

Options:
  --root <dir>         Directory to scan (default: book/src)
  --backend <b>        Backend: auto|mmdc|npx (default: auto)
  --mmdc <cmd>         Mermaid CLI executable (default: mmdc)
  --npx-pkg <pkg>      Package for npx backend (default: @mermaid-js/mermaid-cli@9.2.2)
  --verbose            Print success summary
  --fail-fast          Stop after first reported error
  --keep-temp          Keep temporary .mmd/.svg files for debugging
  -h, --help           Show this help
"
    );
}

fn collect_md_files(root: &Path, out: &mut Vec<PathBuf>) -> std::io::Result<()> {
    if root.is_file() {
        if root.extension().and_then(|s| s.to_str()) == Some("md") {
            out.push(root.to_path_buf());
        }
        return Ok(());
    }
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            collect_md_files(&path, out)?;
        } else if path.extension().and_then(|s| s.to_str()) == Some("md") {
            out.push(path);
        }
    }
    Ok(())
}

fn extract_mermaid_blocks(text: &str, source_path: &str) -> Vec<MermaidBlock> {
    let mut out = Vec::new();
    let mut lines = text.lines();
    let mut line_no = 0usize;

    while let Some(line) = lines.next() {
        line_no += 1;
        let trimmed = line.trim_start();
        if trimmed == "```mermaid" {
            let fence_line = line_no;
            let content_start_line = line_no + 1;
            let mut body = String::new();
            let mut end_line = line_no;
            while let Some(body_line) = lines.next() {
                line_no += 1;
                let t = body_line.trim_start();
                if t == "```" {
                    end_line = line_no - 1;
                    break;
                }
                body.push_str(body_line);
                body.push('\n');
            }
            // Trim trailing newline to keep parse messages more stable.
            if body.ends_with('\n') {
                body.pop();
            }
            out.push(MermaidBlock {
                source_path: source_path.to_string(),
                fence_line,
                content_start_line,
                content_end_line: end_line.max(content_start_line),
                text: body,
            });
        }
    }
    out
}

fn validate_blocks(
    blocks: &[MermaidBlock],
    backend: Backend,
    mmdc_cmd: &str,
    npx_pkg: &str,
    keep_temp: bool,
    verbose: bool,
    fail_fast: bool,
) -> Result<(), String> {
    let selected = select_backend(backend, mmdc_cmd, verbose)?;

    let pid = std::process::id();
    let tmp_root = std::env::temp_dir().join(format!("annoqc_check_mermaid_{pid}"));
    fs::create_dir_all(&tmp_root).map_err(|e| format!("create temp dir: {e}"))?;

    let mut error_count = 0usize;
    for (idx, b) in blocks.iter().enumerate() {
        let in_path = tmp_root.join(format!("block_{idx}.mmd"));
        let out_path = tmp_root.join(format!("block_{idx}.svg"));
        fs::write(&in_path, &b.text)
            .map_err(|e| format!("write temp file {}: {e}", in_path.display()))?;

        let out = match selected {
            SelectedBackend::Mmdc => Command::new(mmdc_cmd)
                .arg("-i")
                .arg(&in_path)
                .arg("-o")
                .arg(&out_path)
                .arg("-q")
                // Keep this light: no background rendering tweaks, just parse+render.
                .output()
                .map_err(|e| format!("failed to run {mmdc_cmd}: {e}"))?,
            SelectedBackend::Npx => Command::new("npx")
                .arg("-y")
                .arg(npx_pkg)
                .arg("-i")
                .arg(&in_path)
                .arg("-o")
                .arg(&out_path)
                .arg("-q")
                .output()
                .map_err(|e| format!("failed to run npx: {e}"))?,
        };

        if !out.status.success() {
            error_count += 1;
            eprintln!(
                "{}:{} (mermaid block, fence at line {}, lines {}..={}) -> mermaid-cli failed",
                b.source_path,
                b.content_start_line,
                b.fence_line,
                b.content_start_line,
                b.content_end_line
            );
            let stderr = String::from_utf8_lossy(&out.stderr);
            let stdout = String::from_utf8_lossy(&out.stdout);
            if !stderr.trim().is_empty() {
                eprintln!("{stderr}");
            } else if !stdout.trim().is_empty() {
                eprintln!("{stdout}");
            }
            if fail_fast {
                if !keep_temp {
                    let _ = fs::remove_dir_all(&tmp_root);
                } else {
                    eprintln!("temp preserved at {}", tmp_root.display());
                }
                return Err("mermaid syntax errors found".to_string());
            }
        }
    }

    if !keep_temp {
        let _ = fs::remove_dir_all(&tmp_root);
    } else if verbose {
        println!("temp preserved at {}", tmp_root.display());
    }

    if error_count == 0 {
        if verbose {
            println!("ok: {} mermaid blocks", blocks.len());
        }
        Ok(())
    } else {
        Err(format!("mermaid syntax errors found: {error_count}"))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelectedBackend {
    Mmdc,
    Npx,
}

fn select_backend(
    backend: Backend,
    mmdc_cmd: &str,
    verbose: bool,
) -> Result<SelectedBackend, String> {
    match backend {
        Backend::Mmdc => {
            ensure_mmdc_exists(mmdc_cmd)?;
            Ok(SelectedBackend::Mmdc)
        }
        Backend::Npx => {
            ensure_npx_exists()?;
            Ok(SelectedBackend::Npx)
        }
        Backend::Auto => match Command::new(mmdc_cmd)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
        {
            Ok(_) => {
                if verbose {
                    println!("backend: mmdc (from PATH)");
                }
                Ok(SelectedBackend::Mmdc)
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                ensure_npx_exists()?;
                if verbose {
                    println!("backend: npx (mmdc not found in PATH)");
                }
                Ok(SelectedBackend::Npx)
            }
            Err(e) => Err(format!("failed to probe `{mmdc_cmd} --version`: {e}")),
        },
    }
}

fn ensure_mmdc_exists(mmdc_cmd: &str) -> Result<(), String> {
    match Command::new(mmdc_cmd)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
    {
        Ok(_) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(format!(
            "`{mmdc_cmd}` not found in PATH. Either:\n  - Install globally: npm i -g @mermaid-js/mermaid-cli\n  - Or use the npx backend: cargo run --bin check_mermaid -- --backend npx"
        )),
        Err(e) => Err(format!("failed to run `{mmdc_cmd} --version`: {e}")),
    }
}

fn ensure_npx_exists() -> Result<(), String> {
    match Command::new("npx")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
    {
        Ok(_) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(
            "`npx` not found in PATH (need Node.js/npm). Install Node.js or provide --backend mmdc."
                .to_string(),
        ),
        Err(e) => Err(format!("failed to run `npx --version`: {e}")),
    }
}
