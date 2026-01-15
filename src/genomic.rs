use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::gff;
use noodles::gff::record::Strand;

#[derive(Debug, Clone, Default)]
pub struct GenomicMetrics {
    #[allow(dead_code)]
    pub gene_id: String,
    pub introns_total: usize,
    pub splice_canonical: usize,    // GT..AG
    pub splice_major_noncan: usize, // GC..AG
    pub splice_minor: usize,        // AT..AC
    pub splice_weird: usize,        // Everything else
    pub intron_len_min: usize,
    pub intron_len_max: usize,
    pub intron_len_avg: f64,
}

#[derive(Debug)]
struct Exon {
    start: usize, // 1-based inclusive
    end: usize,
    strand: Strand,
}

#[derive(Default)]
struct GeneAccum {
    introns_total: usize,
    splice_canonical: usize,
    splice_major_noncan: usize,
    splice_minor: usize,
    splice_weird: usize,
    intron_len_min: usize,
    intron_len_max: usize,
    intron_len_total: usize,
}

fn split_attr_list(value: &gff::record::attributes::field::Value) -> impl Iterator<Item = &str> {
    value
        .iter()
        .flat_map(|s| s.split(','))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
}

fn as_2mer(b: &[u8]) -> Option<[u8; 2]> {
    let first = *b.first()?;
    let second = *b.get(1)?;
    Some([first, second])
}

fn revcomp_2mer(b: [u8; 2]) -> [u8; 2] {
    fn comp(x: u8) -> u8 {
        match x {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        }
    }
    [comp(b[1]), comp(b[0])]
}

fn query_2mer<R>(
    reader: &mut fasta::io::IndexedReader<R>,
    chrom: &str,
    start: usize,
) -> Option<[u8; 2]>
where
    R: std::io::BufRead + std::io::Seek,
{
    if start == 0 {
        return None;
    }
    let end = start.checked_add(1)?;
    let start_pos = Position::try_from(start).ok()?;
    let end_pos = Position::try_from(end).ok()?;
    let region = Region::new(chrom, start_pos..=end_pos);
    let record = reader.query(&region).ok()?;
    let seq = record.sequence().as_ref();
    if seq.len() == 2 {
        as_2mer(seq)
    } else {
        None
    }
}

/// Analyze GFF gene context against a reference genome and return per-gene metrics.
///
/// `gff_path` must point to a GFF3 file with exon/CDS features. `genome_path`
/// should reference the corresponding FASTA; if a `.fai` index is present it is
/// used, otherwise sequences are loaded into memory.
/// Returns intron and splice-site statistics aggregated by gene, propagating
/// any parsing or I/O errors from the underlying readers.
pub fn analyze_gff_context(
    gff_path: &str,
    genome_path: &str,
) -> Result<HashMap<String, GenomicMetrics>, Box<dyn std::error::Error>> {
    // 1. Index Genome if a FASTA index is available; otherwise load into memory.
    let mut genome_map = HashMap::new();
    let mut indexed_reader = None;
    let fai_path = format!("{}.fai", genome_path);
    if Path::new(&fai_path).exists() {
        if let Ok(reader) =
            fasta::io::indexed_reader::Builder::default().build_from_path(genome_path)
        {
            indexed_reader = Some(reader);
        }
    }
    if indexed_reader.is_none() {
        let mut reader = fasta::io::Reader::new(BufReader::new(File::open(genome_path)?));
        for result in reader.records() {
            let record = result?;
            let name = String::from_utf8_lossy(record.name()).to_string();
            // Keep just the first word of the header as ID
            let id = name.split_whitespace().next().unwrap_or(&name).to_string();
            genome_map.insert(id, record.sequence().as_ref().to_vec());
        }
    }

    // 2. Parse GFF
    let mut reader = gff::io::Reader::new(BufReader::new(File::open(gff_path)?));

    // gene_id -> list of Exons
    // We assume mRNA/transcript structure. We group exons by Parent.
    let mut transcript_exons: HashMap<String, Vec<(String, Exon)>> = HashMap::new(); // TranscriptID -> (Chrom, Exon)
    let mut transcript_to_gene: HashMap<String, String> = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let ty = record.ty(); // e.g. "exon" or "CDS"
        let attributes = record.attributes();
        if ty == "mRNA" || ty == "transcript" {
            if let Some(id_raw) = attributes.get("ID") {
                let parent = attributes
                    .get("Parent")
                    .and_then(|value| split_attr_list(value).next());
                if let Some(parent) = parent {
                    for id in split_attr_list(id_raw) {
                        transcript_to_gene.insert(id.to_string(), parent.to_string());
                    }
                }
            }
        }
        if ty == "exon" || ty == "CDS" {
            if let Some(parent) = attributes.get("Parent") {
                let chrom =
                    String::from_utf8_lossy(record.reference_sequence_name().as_ref()).to_string();
                let start: usize = record.start().into();
                let end: usize = record.end().into();
                let strand = record.strand();
                for parent_id in split_attr_list(parent) {
                    transcript_exons
                        .entry(parent_id.to_string())
                        .or_default()
                        .push((chrom.clone(), Exon { start, end, strand }));
                }
            }
        }
    }

    // 3. Analyze Introns
    let mut transcript_metrics: HashMap<String, (GenomicMetrics, usize)> = HashMap::new();

    for (tx_id, mut exons_with_chrom) in transcript_exons {
        if exons_with_chrom.len() < 2 {
            // Single exon gene, no introns
            transcript_metrics.insert(
                tx_id.clone(),
                (
                    GenomicMetrics {
                        gene_id: tx_id,
                        ..Default::default()
                    },
                    0,
                ),
            );
            continue;
        }

        // Sort by start coord
        exons_with_chrom.sort_by_key(|(_, e)| e.start);

        // Check they are on same chrom
        let chrom = &exons_with_chrom[0].0;
        if !exons_with_chrom.iter().all(|(c, _)| c == chrom) {
            continue; // Skip trans-splicing weirdness
        }

        let mut m = GenomicMetrics {
            gene_id: tx_id.clone(),
            intron_len_min: usize::MAX,
            ..Default::default()
        };

        let seq = genome_map.get(chrom);
        let mut total_len = 0;

        for i in 0..exons_with_chrom.len() - 1 {
            let e1 = &exons_with_chrom[i].1;
            let e2 = &exons_with_chrom[i + 1].1;

            // Intron is between e1.end and e2.start
            // GFF is 1-based inclusive.
            // Intron start = e1.end + 1
            // Intron end = e2.start - 1

            if e2.start <= e1.end {
                continue; // Overlapping exons?
            }

            let i_start = e1.end + 1;
            let i_end = e2.start - 1;
            let len = i_end - i_start + 1;

            m.introns_total += 1;
            total_len += len;
            m.intron_len_max = m.intron_len_max.max(len);
            m.intron_len_min = m.intron_len_min.min(len);

            if let Some(reader) = indexed_reader.as_mut() {
                let donor = query_2mer(reader, chrom, i_start);
                let acceptor = if i_end >= 2 {
                    query_2mer(reader, chrom, i_end - 1)
                } else {
                    None
                };
                if let (Some(donor_raw), Some(acceptor_raw)) = (donor, acceptor) {
                    let strand = if e1.strand == e2.strand {
                        e1.strand
                    } else {
                        Strand::Unknown
                    };
                    let splice = match strand {
                        Strand::Forward => Some((donor_raw, acceptor_raw)),
                        Strand::Reverse => {
                            Some((revcomp_2mer(acceptor_raw), revcomp_2mer(donor_raw)))
                        }
                        Strand::None | Strand::Unknown => None,
                    };
                    if let Some((donor, acceptor)) = splice {
                        if donor == *b"GT" && acceptor == *b"AG" {
                            m.splice_canonical += 1;
                        } else if donor == *b"GC" && acceptor == *b"AG" {
                            m.splice_major_noncan += 1;
                        } else if donor == *b"AT" && acceptor == *b"AC" {
                            m.splice_minor += 1;
                        } else {
                            m.splice_weird += 1;
                        }
                    } else {
                        m.splice_weird += 1;
                    }
                } else {
                    m.splice_weird += 1;
                }
            } else if let Some(s) = seq {
                // Get splice sites (first 2, last 2)
                // 0-based index: start-1, end-1
                let idx_start = i_start - 1;
                let idx_end = i_end - 1;

                let donor_end = idx_start.checked_add(2);
                let acceptor_end = idx_end.checked_add(1);
                if idx_end > 0 {
                    let donor_raw = donor_end
                        .and_then(|end| s.get(idx_start..end))
                        .and_then(as_2mer);
                    let acceptor_raw = acceptor_end
                        .and_then(|end| s.get(idx_end - 1..end))
                        .and_then(as_2mer);
                    let strand = if e1.strand == e2.strand {
                        e1.strand
                    } else {
                        Strand::Unknown
                    };
                    let splice = match (donor_raw, acceptor_raw, strand) {
                        (Some(donor), Some(acceptor), Strand::Forward) => Some((donor, acceptor)),
                        (Some(donor), Some(acceptor), Strand::Reverse) => {
                            Some((revcomp_2mer(acceptor), revcomp_2mer(donor)))
                        }
                        _ => None,
                    };
                    if let Some((donor, acceptor)) = splice {
                        if donor == *b"GT" && acceptor == *b"AG" {
                            m.splice_canonical += 1;
                        } else if donor == *b"GC" && acceptor == *b"AG" {
                            m.splice_major_noncan += 1;
                        } else if donor == *b"AT" && acceptor == *b"AC" {
                            m.splice_minor += 1;
                        } else {
                            m.splice_weird += 1;
                        }
                    } else {
                        m.splice_weird += 1;
                    }
                } else {
                    m.splice_weird += 1; // Out of bounds
                }
            }
        }

        if m.introns_total > 0 {
            m.intron_len_avg = total_len as f64 / m.introns_total as f64;
        } else {
            m.intron_len_min = 0;
        }

        transcript_metrics.insert(tx_id, (m, total_len));
    }

    // 4. Aggregate per gene (or transcript if no mapping)
    let mut gene_accum: HashMap<String, GeneAccum> = HashMap::new();
    for (tx_id, (m, total_len)) in transcript_metrics {
        let gene_id = transcript_to_gene.get(&tx_id).cloned().unwrap_or(tx_id);
        let entry = gene_accum.entry(gene_id).or_insert_with(|| GeneAccum {
            intron_len_min: usize::MAX,
            ..Default::default()
        });
        entry.introns_total += m.introns_total;
        entry.splice_canonical += m.splice_canonical;
        entry.splice_major_noncan += m.splice_major_noncan;
        entry.splice_minor += m.splice_minor;
        entry.splice_weird += m.splice_weird;
        entry.intron_len_total += total_len;
        if m.introns_total > 0 {
            entry.intron_len_min = entry.intron_len_min.min(m.intron_len_min);
            entry.intron_len_max = entry.intron_len_max.max(m.intron_len_max);
        }
    }

    let mut metrics_map = HashMap::new();
    for (gene_id, acc) in gene_accum {
        let intron_len_min = if acc.introns_total > 0 {
            acc.intron_len_min
        } else {
            0
        };
        let intron_len_avg = if acc.introns_total > 0 {
            acc.intron_len_total as f64 / acc.introns_total as f64
        } else {
            0.0
        };
        metrics_map.insert(
            gene_id.clone(),
            GenomicMetrics {
                gene_id,
                introns_total: acc.introns_total,
                splice_canonical: acc.splice_canonical,
                splice_major_noncan: acc.splice_major_noncan,
                splice_minor: acc.splice_minor,
                splice_weird: acc.splice_weird,
                intron_len_min,
                intron_len_max: acc.intron_len_max,
                intron_len_avg,
            },
        );
    }

    Ok(metrics_map)
}

#[cfg(test)]
mod tests {
    use super::analyze_gff_context;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn genomic_context_counts_forward_and_reverse_splice_sites(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let tmp = TempDir::new()?;
        let genome_path = tmp.path().join("genome.fa");
        let gff_path = tmp.path().join("genes.gff3");

        let genome = b">chr1\nAAAAGTAGAAAA\n>chr2\nAAAACTACAAAA\n";
        fs::write(&genome_path, genome)?;

        let gff = b"##gff-version 3
chr1\t.\tgene\t1\t12\t.\t+\t.\tID=gene1
chr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1
chr1\t.\texon\t1\t4\t.\t+\t.\tParent=tx1
chr1\t.\texon\t9\t12\t.\t+\t.\tParent=tx1
chr2\t.\tgene\t1\t12\t.\t-\t.\tID=gene2
chr2\t.\tmRNA\t1\t12\t.\t-\t.\tID=tx2;Parent=gene2
chr2\t.\texon\t1\t4\t.\t-\t.\tParent=tx2
chr2\t.\texon\t9\t12\t.\t-\t.\tParent=tx2
";
        fs::write(&gff_path, gff)?;

        let metrics =
            analyze_gff_context(gff_path.to_str().unwrap(), genome_path.to_str().unwrap())?;

        let g1 = metrics.get("gene1").expect("gene1 metrics");
        assert_eq!(g1.introns_total, 1);
        assert_eq!(g1.splice_canonical, 1);

        let g2 = metrics.get("gene2").expect("gene2 metrics");
        assert_eq!(g2.introns_total, 1);
        assert_eq!(g2.splice_canonical, 1);

        Ok(())
    }
}
