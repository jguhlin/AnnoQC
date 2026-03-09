use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::scoring_thresholds::classification_base;
use crate::*;
use std::collections::HashMap;

pub(crate) fn export_high_sequences(
    out_dir: &str,
    override_path: Option<&str>,
    metrics: &[GeneMetrics],
    intrinsic: &HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>,
    scores_map: &HashMap<String, (f64, String)>,
) -> Result<Option<(PathBuf, usize)>, Box<dyn std::error::Error>> {
    let target_path = override_path
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(out_dir).join("high_scoring.faa"));
    let mut writer: Option<BufWriter<File>> = None;
    let mut count = 0usize;
    for m in metrics {
        let Some((_, classif)) = scores_map.get(&m.gene_id) else {
            continue;
        };
        if classification_base(classif) != "High" {
            continue;
        }
        let Some((_im, seq)) = intrinsic.get(&m.gene_id) else {
            continue;
        };
        if writer.is_none() {
            if let Some(parent) = target_path.parent() {
                if !parent.as_os_str().is_empty() {
                    fs::create_dir_all(parent)?;
                }
            }
            writer = Some(BufWriter::new(File::create(&target_path)?));
        }
        if let Some(buf) = writer.as_mut() {
            write_fasta_record(buf, &m.gene_id, seq)?;
            count += 1;
        }
    }
    if let Some(mut buf) = writer {
        buf.flush()?;
        if count == 0 {
            drop(buf);
            fs::remove_file(&target_path).ok();
            Ok(None)
        } else {
            Ok(Some((target_path, count)))
        }
    } else {
        Ok(None)
    }
}

fn write_fasta_record<W: Write>(writer: &mut W, gene_id: &str, seq: &[u8]) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(gene_id.as_bytes())?;
    writer.write_all(b"\n")?;
    for chunk in seq.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
