use needletail::parse_fastx_file;

/// Very simple longest-ORF translator (one frame) placeholder to enable nt input.
/// Writes a temporary protein FASTA and returns its path.
pub fn translate_nt_fasta_to_protein(input_nt_fasta: &str) -> Result<String, String> {
    let tmp = std::path::Path::new(input_nt_fasta)
        .with_extension("translated.fa");
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

fn translate_longest_orf(seq: &[u8]) -> String {
    // Very naive: translate frame 0 only; stop at '*' and pick longest stretch
    let table = |codon: &[u8]| -> char {
        match codon {
            b"TTT"|b"TTC"=> 'F', b"TTA"|b"TTG"|b"CTT"|b"CTC"|b"CTA"|b"CTG"=> 'L',
            b"ATT"|b"ATC"|b"ATA"=> 'I', b"ATG"=> 'M', b"GTT"|b"GTC"|b"GTA"|b"GTG"=> 'V',
            b"TCT"|b"TCC"|b"TCA"|b"TCG"|b"AGT"|b"AGC"=> 'S', b"CCT"|b"CCC"|b"CCA"|b"CCG"=> 'P',
            b"ACT"|b"ACC"|b"ACA"|b"ACG"=> 'T', b"GCT"|b"GCC"|b"GCA"|b"GCG"=> 'A',
            b"TAT"|b"TAC"=> 'Y', b"TAA"|b"TAG"|b"TGA"=> '*', b"CAT"|b"CAC"=> 'H',
            b"CAA"|b"CAG"=> 'Q', b"AAT"|b"AAC"=> 'N', b"AAA"|b"AAG"=> 'K',
            b"GAT"|b"GAC"=> 'D', b"GAA"|b"GAG"=> 'E', b"TGT"|b"TGC"=> 'C',
            b"TGG"=> 'W', b"CGT"|b"CGC"|b"CGA"|b"CGG"|b"AGA"|b"AGG"=> 'R',
            b"GGT"|b"GGC"|b"GGA"|b"GGG"=> 'G', _ => 'X'
        }
    };
    let mut aa: Vec<char> = Vec::new();
    let mut i = 0usize;
    while i + 3 <= seq.len() {
        let codon = &seq[i..i+3];
        let mut cod = codon.to_ascii_uppercase();
        for b in &mut cod { if *b == b'U' { *b = b'T'; } }
        aa.push(table(&cod));
        i += 3;
    }
    let s: String = aa.into_iter().collect();
    // choose longest segment between stops
    s.split('*').max_by_key(|seg| seg.len()).unwrap_or("").to_string()
}
