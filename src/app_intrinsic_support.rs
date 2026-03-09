use std::collections::{HashMap, HashSet};

use needletail::parse_fastx_file;

use crate::metrics::compute_intrinsic;
use crate::{metrics, GeneMetrics};

#[allow(clippy::type_complexity)]
pub(crate) fn compute_intrinsic_for_ids(
    fasta_path: &str,
    metrics: &[GeneMetrics],
) -> Result<HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)>, Box<dyn std::error::Error>> {
    let wanted: HashSet<String> = metrics.iter().map(|m| m.gene_id.clone()).collect();
    let mut map: HashMap<String, (metrics::IntrinsicMetrics, Vec<u8>)> = HashMap::new();
    let mut reader = parse_fastx_file(fasta_path)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let gid = id.split_whitespace().next().unwrap_or("").to_string();
        if !wanted.contains(&gid) {
            continue;
        }
        let seq = rec.seq().to_vec();
        let im = compute_intrinsic(&seq);
        map.insert(gid, (im, seq));
    }
    Ok(map)
}
