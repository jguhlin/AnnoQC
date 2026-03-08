use needletail::parse_fastx_file;

/// Multi-frame ORF translator for nucleotide input.
///
/// Scans all 6 reading frames (3 forward, 3 reverse complement) and returns
/// the longest open reading frame (ATG → STOP). Each sequence is processed
/// independently and the longest ORF across all frames is returned.
/// Writes a temporary protein FASTA and returns its path; caller must clean it up.
pub fn translate_nt_fasta_to_protein(input_nt_fasta: &str) -> Result<String, String> {
    let tmp = std::path::Path::new(input_nt_fasta).with_extension("translated.fa");
    let mut out = std::fs::File::create(&tmp).map_err(|e| e.to_string())?;
    let mut reader = parse_fastx_file(input_nt_fasta).map_err(|e| e.to_string())?;
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| e.to_string())?;
        let id = String::from_utf8_lossy(rec.id());
        let seq = rec.seq();
        let aa = translate_longest_orf(&seq);
        use std::io::Write as _;
        writeln!(out, ">{}", id).map_err(|e| e.to_string())?;
        writeln!(out, "{}", aa).map_err(|e| e.to_string())?;
    }
    Ok(tmp.to_string_lossy().to_string())
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

fn translate_longest_orf(seq: &[u8]) -> String {
    let mut best = String::new();

    // Check all 3 forward frames
    for frame in 0..3 {
        let orf = find_orf_in_frame(seq, frame);
        if orf.len() > best.len() {
            best = orf;
        }
    }

    // Check all 3 reverse frames
    let rc = reverse_complement(seq);
    for frame in 0..3 {
        let orf = find_orf_in_frame(&rc, frame);
        if orf.len() > best.len() {
            best = orf;
        }
    }

    best
}

fn find_orf_in_frame(seq: &[u8], frame_offset: usize) -> String {
    if seq.len() < frame_offset + 3 {
        return String::new();
    }

    fn normalize_base(base: u8) -> u8 {
        let upper = base.to_ascii_uppercase();
        if upper == b'U' {
            b'T'
        } else {
            upper
        }
    }

    fn translate_codon(a: u8, b: u8, c: u8) -> char {
        match [a, b, c] {
            [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => 'F',
            [b'T', b'T', b'A']
            | [b'T', b'T', b'G']
            | [b'C', b'T', b'T']
            | [b'C', b'T', b'C']
            | [b'C', b'T', b'A']
            | [b'C', b'T', b'G'] => 'L',
            [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => 'I',
            [b'A', b'T', b'G'] => 'M',
            [b'G', b'T', b'T'] | [b'G', b'T', b'C'] | [b'G', b'T', b'A'] | [b'G', b'T', b'G'] => {
                'V'
            }
            [b'T', b'C', b'T']
            | [b'T', b'C', b'C']
            | [b'T', b'C', b'A']
            | [b'T', b'C', b'G']
            | [b'A', b'G', b'T']
            | [b'A', b'G', b'C'] => 'S',
            [b'C', b'C', b'T'] | [b'C', b'C', b'C'] | [b'C', b'C', b'A'] | [b'C', b'C', b'G'] => {
                'P'
            }
            [b'A', b'C', b'T'] | [b'A', b'C', b'C'] | [b'A', b'C', b'A'] | [b'A', b'C', b'G'] => {
                'T'
            }
            [b'G', b'C', b'T'] | [b'G', b'C', b'C'] | [b'G', b'C', b'A'] | [b'G', b'C', b'G'] => {
                'A'
            }
            [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => 'Y',
            [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A'] => '*',
            [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => 'H',
            [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => 'Q',
            [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => 'N',
            [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => 'K',
            [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => 'D',
            [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => 'E',
            [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => 'C',
            [b'T', b'G', b'G'] => 'W',
            [b'C', b'G', b'T']
            | [b'C', b'G', b'C']
            | [b'C', b'G', b'A']
            | [b'C', b'G', b'G']
            | [b'A', b'G', b'A']
            | [b'A', b'G', b'G'] => 'R',
            [b'G', b'G', b'T'] | [b'G', b'G', b'C'] | [b'G', b'G', b'A'] | [b'G', b'G', b'G'] => {
                'G'
            }
            _ => 'X',
        }
    }

    let mut best_orf = String::new();
    let mut current_orf = String::new();
    let mut in_orf = false;

    let mut i = frame_offset;
    while i + 3 <= seq.len() {
        let a = normalize_base(seq[i]);
        let b = normalize_base(seq[i + 1]);
        let c = normalize_base(seq[i + 2]);
        let codon = [a, b, c];

        if codon == [b'A', b'T', b'G'] {
            // Start codon ATG
            if !in_orf {
                // New ORF starting
                in_orf = true;
                current_orf.clear();
                current_orf.push('M');
            } else {
                // ATG within ORF, just add M
                current_orf.push('M');
            }
        } else if in_orf {
            let aa = translate_codon(a, b, c);
            if aa == '*' {
                // Stop codon - end this ORF
                if current_orf.len() > best_orf.len() {
                    best_orf = std::mem::take(&mut current_orf);
                }
                in_orf = false;
            } else {
                current_orf.push(aa);
            }
        }
        i += 3;
    }

    // Check if we're still in an ORF at the end
    if in_orf && current_orf.len() > best_orf.len() {
        best_orf = current_orf;
    }

    best_orf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"ATGAAATAA"), b"TTATTTCAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GGGG"), b"CCCC");
    }

    #[test]
    fn test_forward_frame_0() {
        // ATGAAATAA: ATG at position 0, TAA stop
        // Should translate to MK
        let result = find_orf_in_frame(b"ATGAAATAA", 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_forward_frame_1() {
        // AATGAAATAA: ATG at position 1, TAA stop
        // Should translate to MK
        let result = find_orf_in_frame(b"AATGAAATAA", 1);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_forward_frame_2() {
        // GGATGAAATAA: ATG at position 2 (0=G,1=G,2=A,3=T,4=G)
        // Frame 2 starts at position 2: ATGAAATAA -> MK
        let result = find_orf_in_frame(b"GGATGAAATAA", 2);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_reverse_frame() {
        // TTATTTCAT is reverse complement of ATGAAATAA
        // Should find ORF and translate to MK
        let seq = b"TTATTTCAT";
        let rc = reverse_complement(seq);
        let result = find_orf_in_frame(&rc, 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_no_start_codon() {
        // No ATG start codon = empty result
        let result = find_orf_in_frame(b"AAATAAGGG", 0);
        assert_eq!(result, "");
    }

    #[test]
    fn test_short_sequence() {
        // < 3 bases = empty
        assert_eq!(find_orf_in_frame(b"AT", 0), "");
        assert_eq!(find_orf_in_frame(b"A", 0), "");
        assert_eq!(find_orf_in_frame(b"", 0), "");
    }

    #[test]
    fn test_chooses_longest_orf() {
        // ATGAAATAAGGGATGCCCCCCTAA
        // ORF1: MK (from ATGAAATAA)
        // ORF2: MPP (from ATGCCCCCCTAA) - 3 aa total
        // MPP is longer (3 aa vs 2 aa)
        let result = translate_longest_orf(b"ATGAAATAAGGGATGCCCCCCTAA");
        assert_eq!(result, "MPP");
    }

    #[test]
    fn test_partial_codon_at_end() {
        // ATGAAA: ATG = M, AAA = K, partial codon handled
        let result = find_orf_in_frame(b"ATGAAA", 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_multiple_orfs_same_frame() {
        // ATGAAATAAATGTTTTAA
        // ORF1: MK (ATGAAATAA) - 2 aa
        // ORF2: MF (ATGTTTTAA) - 2 aa (same length, first wins)
        let result = find_orf_in_frame(b"ATGAAATAAATGTTTTAA", 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_orf_without_stop() {
        // ATGAAA - no stop codon, should return MK
        let result = find_orf_in_frame(b"ATGAAA", 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_all_frames_checked() {
        // Sequence with ORF in frame 2
        // GGATGAAATAA = GG (padding) + ATG (start at position 2) + AAATAA
        let seq = b"GGATGAAATAA";
        let result = translate_longest_orf(seq);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_non_standard_bases() {
        // ATGNNNTAA - N should normalize and still find ORF
        let result = find_orf_in_frame(b"ATGNNNTAA", 0);
        // ATG = M, NNN = X (one codon), TAA = stop
        assert_eq!(result, "MX");
    }

    #[test]
    fn test_reverse_complement_finds_orf() {
        // TTATTTCAT -> reverse complement is ATGAAATAA (has ORF)
        let seq = b"TTATTTCAT";
        let result = translate_longest_orf(seq);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_lowercase_bases() {
        // atgaaataa should work same as ATGAAATAA
        let result = find_orf_in_frame(b"atgaaataa", 0);
        assert_eq!(result, "MK");
    }

    #[test]
    fn test_mixed_case() {
        // AtGaaATaA should work
        let result = find_orf_in_frame(b"AtGaaATaA", 0);
        assert_eq!(result, "MK");
    }
}
