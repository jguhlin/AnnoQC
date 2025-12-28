use std::path::Path;
use std::time::Instant;

use serde::ser::SerializeMap;

fn log_json_event<F>(buf: &mut Vec<u8>, build: F)
where
    F: FnOnce(&mut Vec<u8>) -> Result<(), serde_json::Error>,
{
    buf.clear();
    if let Err(e) = build(buf) {
        log::warn!("failed to serialize log event: {}", e);
        return;
    }
    log::info!("{}", String::from_utf8_lossy(buf));
}

pub fn run_step<F>(
    done_path: &Path,
    resume: bool,
    name: &str,
    log_json: bool,
    f: F,
) -> Result<(), String>
where
    F: FnOnce() -> Result<(), String>,
{
    let mut json_buf = if log_json {
        Vec::with_capacity(128)
    } else {
        Vec::new()
    };

    if resume && done_path.exists() {
        if log_json {
            let done_display = done_path.display().to_string();
            log_json_event(&mut json_buf, |buf| {
                let mut ser = serde_json::Serializer::new(buf);
                let mut map = ser.serialize_map(Some(3))?;
                map.serialize_entry("event", "resume_skip")?;
                map.serialize_entry("step", name)?;
                map.serialize_entry("done", &done_display)?;
                map.end()
            });
        } else {
            log::info!(
                "resume: skipping step '{}' ({} exists)",
                name,
                done_path.display()
            );
        }
        return Ok(());
    }
    if log_json {
        log_json_event(&mut json_buf, |buf| {
            let mut ser = serde_json::Serializer::new(buf);
            let mut map = ser.serialize_map(Some(2))?;
            map.serialize_entry("event", "step_start")?;
            map.serialize_entry("step", name)?;
            map.end()
        });
    } else {
        log::info!("step '{}' started", name);
    }
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
    if log_json {
        let duration_sec = dt.as_secs_f64();
        log_json_event(&mut json_buf, |buf| {
            let mut ser = serde_json::Serializer::new(buf);
            let mut map = ser.serialize_map(Some(3))?;
            map.serialize_entry("event", "step_finish")?;
            map.serialize_entry("step", name)?;
            map.serialize_entry("duration_sec", &duration_sec)?;
            map.end()
        });
    } else {
        log::info!("step '{}' finished in {:.2?}", name, dt);
    }
    Ok(())
}
