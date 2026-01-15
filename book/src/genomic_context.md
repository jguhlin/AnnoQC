# Genomic Context (GFF3 + Genome)

AnnoQC can compute splice-site and intron statistics from a GFF3 annotation paired with a genome FASTA.

## Enable the pillar

```
annoqc analyze \
  --fasta input.faa \
  --db uniprot_sprot.dmnd \
  --gff annotations.gff3 \
  --genome genome.fa
```

The genomic metrics appear in JSONL under the `genomic` block and contribute to `genomic_score` when the pillar is enabled.

## What is computed

For each gene (aggregated across isoforms):

- `introns_total`
- `splice_canonical` (GT..AG)
- `splice_major_noncan` (GC..AG)
- `splice_minor` (AT..AC)
- `splice_weird` (everything else / unknown)
- `intron_len_min`, `intron_len_max`, `intron_len_avg`

Splice classification respects strand. On `-` strand, donor/acceptor motifs are reverse-complemented before classification.

## How records are interpreted

- Exons are collected from `exon` or `CDS` features.
- Exons are grouped by `Parent` (transcript ID).
- Transcripts map to genes via `mRNA` or `transcript` records that contain `ID` and `Parent`.
- Metrics are aggregated at the gene level by summing introns across transcripts and summarizing lengths.

## Caveats and gotchas

- If exon `Parent` or transcript `ID/Parent` attributes are missing, those records are skipped.
- If exons for a transcript span multiple contigs, the transcript is skipped.
- If strand is inconsistent across exons, the splice motif is counted as `splice_weird`.
- If the genome FASTA lacks a `.fai` index, the entire genome is loaded into memory.
- Contig naming must match between the GFF3 and genome FASTA.

## Tips

- Build a FASTA index to speed random access:
  - `samtools faidx genome.fa`
- Normalize GFF3 attributes (`ID`, `Parent`) before running AnnoQC.
