use std::path::Path;
use std::time::Instant;

pub fn run_step<F>(done_path: &Path, resume: bool, name: &str, f: F) -> Result<(), String>
where
    F: FnOnce() -> Result<(), String>,
{
    if resume && done_path.exists() {
        log::info!(
            "resume: skipping step '{}' ({} exists)",
            name,
            done_path.display()
        );
        return Ok(());
    }
    log::info!("step '{}' started", name);
    let t0 = Instant::now();
    f()?;
    let dt = t0.elapsed();
    if let Some(parent) = done_path.parent() {
        if let Err(e) = std::fs::create_dir_all(parent) {
            log::warn!(
                "could not create done dir for {}: {}",
                done_path.display(),
                e
            );
        }
    }
    if let Err(e) = std::fs::write(done_path, b"") {
        log::warn!("could not write done file {}: {}", done_path.display(), e);
    }
    log::info!("step '{}' finished in {:.2?}", name, dt);
    Ok(())
}
