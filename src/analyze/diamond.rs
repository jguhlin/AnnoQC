use crate::config_types::FileConfig;
use crate::*;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

pub(in crate::analyze) struct DiamondStageResult {
    pub(in crate::analyze) diamond_tsv: PathBuf,
    pub(in crate::analyze) diamond_secs: f64,
    pub(in crate::analyze) refprot_grouped: HashMap<String, Vec<diamond::DiamondHitRow>>,
}

pub(in crate::analyze) fn run_diamond_stage(
    cfg: &EffectiveConfig,
    file_cfg: &FileConfig,
    args: &AnalyzeArgs,
    log_json: bool,
) -> Result<DiamondStageResult, Box<dyn std::error::Error>> {
    let dia_cfg = DiamondConfig {
        bin: cfg.diamond_bin.clone(),
        db: cfg.db.clone(),
        query_fasta: cfg.fasta.clone(),
        threads: cfg.threads,
        out_dir: cfg.out.clone(),
        out_name: "diamond.blastp.tsv".to_string(),
        retries: 2,
        max_hsps: args
            .diamond_max_hsps
            .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.max_hsps))
            .unwrap_or(5),
    };
    log::info!(
        "diamond: started mode={:?} out_dir={} out_name={}",
        args.diamond_mode,
        cfg.out,
        dia_cfg.out_name
    );
    let t_diamond = step_start("diamond", log_json);
    let diamond_tsv = match args.diamond_mode {
        DiamondMode::Single => blastp_once(&dia_cfg)?,
        DiamondMode::Batch => diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?,
        DiamondMode::Auto => {
            let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
            let threshold = args
                .diamond_auto_threshold
                .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold))
                .unwrap_or(200_000usize);
            if n > threshold {
                diamond::blastp_chunked(&dia_cfg, args.batch_size.max(1), log_json)?
            } else {
                blastp_once(&dia_cfg)?
            }
        }
    };
    let diamond_secs = step_finish("diamond", t_diamond, log_json);
    log::info!(
        "diamond: finished out={} seconds={:.3}",
        diamond_tsv.display(),
        diamond_secs
    );

    let refprot_map_path = Path::new("share/refprot/aves/refprot_proteome_map.tsv");
    let refprot_proteome_map = load_refprot_proteome_map(refprot_map_path);
    let mut refprot_grouped: HashMap<String, Vec<diamond::DiamondHitRow>> = HashMap::new();
    if let Some(ref_db) = cfg.refprot_db.as_ref() {
        if Path::new(ref_db).exists() {
            let refprot_ok = if let Err(e) = preflight::check_diamond_db(&cfg.diamond_bin, ref_db) {
                log::warn!("refprot db preflight failed; skipping refprot: {}", e);
                false
            } else {
                true
            };
            if refprot_ok {
                let ref_db_size = fs::metadata(ref_db).map(|m| m.len()).unwrap_or(0);
                log::info!(
                    "refprot: db={} size_bytes={} fasta={} mode={:?} out_dir={}",
                    ref_db,
                    ref_db_size,
                    cfg.fasta,
                    args.diamond_mode,
                    cfg.out
                );
                let t_ref = step_start("diamond_refprot", log_json);
                let mut ref_cfg = dia_cfg.clone();
                ref_cfg.out_name = "diamond.refprot.tsv".into();
                ref_cfg.db = ref_db.clone();
                let mut used_chunked = false;
                let ref_tsv_result = match args.diamond_mode {
                    DiamondMode::Batch => {
                        used_chunked = true;
                        diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json)
                    }
                    DiamondMode::Auto => {
                        let n = diamond::estimate_query_count(&cfg.fasta).unwrap_or(0);
                        let threshold = args
                            .diamond_auto_threshold
                            .or_else(|| file_cfg.diamond.as_ref().and_then(|d| d.auto_threshold))
                            .unwrap_or(200_000usize);
                        if n > threshold {
                            used_chunked = true;
                            diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json)
                        } else {
                            blastp_once(&ref_cfg)
                        }
                    }
                    DiamondMode::Single => blastp_once(&ref_cfg),
                };
                if let Ok(ref_tsv) = ref_tsv_result {
                    let mut ref_bytes = fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                    if ref_bytes == 0 && !used_chunked {
                        log::warn!(
                            "refprot blastp produced empty output; retrying in chunked mode"
                        );
                        let _ = diamond::blastp_chunked(&ref_cfg, args.batch_size.max(1), log_json);
                        ref_bytes = fs::metadata(&ref_tsv).map(|m| m.len()).unwrap_or(0);
                    }
                    if ref_bytes == 0 {
                        log::warn!("refprot blastp output is empty; continuing without refprot");
                    } else {
                        refprot_grouped = diamond::parse_tsv_grouped(&ref_tsv, None, Some(100))
                            .unwrap_or_default();
                        annotate_refprot_hits(&mut refprot_grouped, &refprot_proteome_map);
                    }
                } else if let Err(e) = ref_tsv_result {
                    log::warn!("refprot blastp failed; continuing without refprot: {}", e);
                }
                let _ = step_finish("diamond_refprot", t_ref, log_json);
            }
        } else {
            log::warn!("refprot db '{}' not found; skipping", ref_db);
        }
    }

    Ok(DiamondStageResult {
        diamond_tsv,
        diamond_secs,
        refprot_grouped,
    })
}
