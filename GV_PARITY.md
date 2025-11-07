# GeneValidator Parity Features

This file tracks features from the original GeneValidator that are desirable for AnnoQC to implement to achieve feature parity.

## 1. Advanced Structural Variation Detection

- **Goal:** Explicitly detect gene fusions, splits, and internal duplications.
- **Current:** The `coverage_delta` metric provides a hint for fusions/splits.
- **Needed:** A more direct analysis, likely by examining the alignment pattern of multiple DIAMOND High-Scoring Pairs (HSPs) against a single query sequence. This would allow distinguishing a fusion (HSPs map to distinct regions of the query) from an internal duplication (HSPs map to overlapping regions of the query).

## 2. Nucleotide Sequence Input and ORF Validation

- **Goal:** Allow users to provide cDNA/nucleotide FASTA files as input, not just protein sequences.
- **Current:** AnnoQC is protein-focused. The `orf_start_score` is a simple check on protein sequences.
- **Needed:**
    -   Logic to handle nucleotide input.
    -   A robust Open Reading Frame (ORF) finding step.
    -   Validation of the ORF against its homologs (e.g., checking for frameshifts, premature stop codons, or retained introns that are not present in the reference proteins).

## 3. Interactive HTML Report

- **Goal:** Generate a user-friendly, graphical HTML report for each run.
- **Current:** AnnoQC produces machine-readable JSONL and CSV files.
- **Needed:** An HTML report that summarizes the run and provides a detailed page for each gene, visualizing its scores, alignments, and domain structure. This is critical for manual curation and for making the results accessible to a wider audience.

## 4. Detailed Alignment-Based Warnings

- **Goal:** Provide more specific warnings about structural issues based on the multiple sequence alignment.
- **Current:** The `start_class` (`LikelyNTruncated`, etc.) is a good start.
- **Needed:** Extend the analysis of the MAFFT alignment to explicitly identify and flag:
    -   Potentially missing exons (indicated by large gaps in the query sequence relative to the consensus of its homologs).
    -   Potentially retained introns (indicated by large insertions in the query that are absent in homologs).
