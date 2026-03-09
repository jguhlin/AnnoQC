use std::collections::HashMap;
use std::sync::Arc;

use crate::genomic;
use crate::orchestration::{step_finish, step_start};
use crate::rnaseq::{parse_expression_file, parse_quant_table, RnaseqMetrics};

pub fn load_rnaseq_data(
    rnaseq_path: Option<String>,
) -> (bool, Arc<HashMap<String, RnaseqMetrics>>) {
    if let Some(path) = rnaseq_path {
        log::info!("loading RNA-seq data from: {}", path);
        let map = if path.ends_with(".tsv") || path.ends_with(".txt") {
            parse_quant_table(&path).or_else(|_| parse_expression_file(&path))
        } else {
            parse_expression_file(&path)
        };
        match map {
            Ok(m) => {
                log::info!("loaded RNA-seq data for {} genes", m.len());
                (true, Arc::new(m))
            }
            Err(e) => {
                log::warn!("failed to load RNA-seq data: {}", e);
                (false, Arc::new(HashMap::new()))
            }
        }
    } else {
        (false, Arc::new(HashMap::new()))
    }
}

pub fn load_genomic_context(
    gff: Option<&str>,
    genome: Option<&str>,
    log_json: bool,
) -> Option<Arc<HashMap<String, genomic::GenomicMetrics>>> {
    if let (Some(gff), Some(genome)) = (gff, genome) {
        let t_gff = step_start("genomic_context", log_json);
        let result = match genomic::analyze_gff_context(gff, genome) {
            Ok(map) => {
                log::info!("genomic context: loaded {} genes", map.len());
                Some(Arc::new(map))
            }
            Err(e) => {
                log::warn!("genomic context analysis failed: {}", e);
                None
            }
        };
        let _ = step_finish("genomic_context", t_gff, log_json);
        result
    } else {
        None
    }
}
