use needletail::parse_fastx_file;

/// Very simple longest-ORF translator (one frame) placeholder to enable nt input.
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

fn translate_longest_orf(seq: &[u8]) -> String {
    // Very naive: translate frame 0 only; stop at '*' and pick longest stretch.
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

    let mut best = String::new();
    let mut current = String::new();
    current.reserve(seq.len() / 3);
    let mut i = 0usize;
    while i + 3 <= seq.len() {
        let a = normalize_base(seq[i]);
        let b = normalize_base(seq[i + 1]);
        let c = normalize_base(seq[i + 2]);
        let aa = translate_codon(a, b, c);
        if aa == '*' {
            if current.len() > best.len() {
                best.clear();
                best.push_str(&current);
            }
            current.clear();
        } else {
            current.push(aa);
        }
        i += 3;
    }
    if current.len() > best.len() {
        best = current;
    }
    best
}
