use std::path::Path;

use serde::Serialize;

use crate::config_types::FileConfig;
use crate::{filehash_xx64, AnalyzeArgs, CalibrationMode};

#[derive(Debug, Clone, Copy)]
pub(crate) struct CalibrationSettings {
    pub(crate) mode: CalibrationMode,
    pub(crate) min_samples: usize,
    pub(crate) min_unique: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct EffectiveConfig {
    pub(crate) fasta: String,
    pub(crate) db: String,
    pub(crate) out: String,
    pub(crate) threads: usize,
    pub(crate) diamond_bin: String,
    pub(crate) reference_fasta: Option<String>,
    pub(crate) refprot_db: Option<String>,
}

#[derive(Serialize, Clone)]
pub struct Checksums {
    pub fasta_xx64: Option<String>,
    pub db_xx64: Option<String>,
}

pub(crate) fn resolve_effective_config(
    file: &FileConfig,
    args: &AnalyzeArgs,
) -> Result<EffectiveConfig, Box<dyn std::error::Error>> {
    let fasta = args
        .fasta
        .clone()
        .or_else(|| file.fasta.clone())
        .ok_or("--fasta or config fasta required")?;
    let db = args
        .db
        .clone()
        .or_else(|| file.db.clone())
        .ok_or("--db or config db required")?;
    let out = args
        .out
        .clone()
        .or_else(|| file.out.clone())
        .unwrap_or_else(|| "results".to_string());
    let threads = args.threads;
    let diamond_bin = args
        .diamond_bin
        .clone()
        .or_else(|| file.diamond_bin.clone())
        .unwrap_or_else(|| "diamond".to_string());
    let reference_fasta = args
        .reference_fasta
        .clone()
        .or_else(|| file.reference_fasta.clone());
    let refprot_db = if let Some(db) = args.refprot_db.clone() {
        Some(db)
    } else {
        let default = Path::new("share/refprot/aves/aves_refprot.dmnd");
        if default.exists() {
            Some(default.to_string_lossy().to_string())
        } else {
            None
        }
    };

    Ok(EffectiveConfig {
        fasta,
        db,
        out,
        threads,
        diamond_bin,
        reference_fasta,
        refprot_db,
    })
}

pub(crate) fn collect_checksums(
    cfg: &EffectiveConfig,
) -> Result<Checksums, Box<dyn std::error::Error>> {
    let fasta = Path::new(&cfg.fasta);
    let db = Path::new(&cfg.db);
    Ok(Checksums {
        fasta_xx64: if fasta.exists() {
            Some(filehash_xx64(fasta)?)
        } else {
            None
        },
        db_xx64: if db.exists() {
            Some(filehash_xx64(db)?)
        } else {
            None
        },
    })
}

pub(crate) fn resolve_calibration_settings(
    args: &AnalyzeArgs,
    file_cfg: &FileConfig,
) -> CalibrationSettings {
    let mode = args
        .calibration_mode
        .or_else(|| file_cfg.calibration.as_ref().and_then(|c| c.mode))
        .unwrap_or(CalibrationMode::Off);
    let min_samples = file_cfg
        .calibration
        .as_ref()
        .and_then(|c| c.min_samples)
        .unwrap_or(50)
        .max(1);
    let min_unique = file_cfg
        .calibration
        .as_ref()
        .and_then(|c| c.min_unique)
        .unwrap_or(5)
        .max(1);
    CalibrationSettings {
        mode,
        min_samples,
        min_unique,
    }
}
