# Project Ideas

## Why
Speed up and normalize annotation QC

## GeneValidator
https://link.springer.com/protocol/10.1007/978-1-4939-9173-0_16/figures/1

https://academic.oup.com/bioinformatics/article/32/10/1559/1742817

https://github.com/wurmlab/genevalidator

### GV Pre-process
* Input Format Verification
* BLAST against reference
* Extract FASTA Seq for all BLAST hits

### GV Analyses
For each input seq

* Length (cluster and rank validations)
* Coverage
* Conserved Regions - MAFFT Alignment
* Different Genes
* Open Reading Frames (ab initio and simliarity based validations, nt seqs only)

## Proposed AnnoQC

## Status (implemented vs ideas)

This page started as a brainstorming list. As the project matured, many of these ideas are now implemented and documented elsewhere in this book. The table below tracks what is **implemented**, what is **diagnostic-only/partial**, and what is **not yet implemented**.

Legend:
- Implemented: exists in the pipeline and is emitted (and often scored).
- Partial / diagnostics: exists, but may be diagnostics-only or not fully integrated into scoring.
- Not yet: not currently implemented in AnnoQC.

| Area / idea | Status | Where to look (docs + code) |
| --- | --- | --- |
| Input format verification / preflight | Implemented | Docs: `Analyze` • Code: `src/preflight.rs` |
| DIAMOND homology search (vs BLAST) | Implemented | Docs: `Analyze`, `Scoring` • Code: `src/diamond.rs` |
| Homology scoring from hit stats | Implemented | Docs: `Scoring` • Code: `src/scoring.rs::compute_homology_score` |
| DIAMOND HSP aggregation / “consensus coverage” | Implemented | Docs: `Scoring` • Code: `src/consensus.rs::{aggregate_hits_by_subject, select_panel_with_result}` |
| Subject coverage pillar (`top_scov`) | Implemented | Docs: `Scoring` • Code: `src/scoring.rs::compute_subject_cov_score` |
| Protein length consistency vs homolog panel | Implemented | Docs: `Scoring` • Code: `src/length.rs::compute_length_consistency` |
| Expected length range from homology panel (min/max) | Implemented | Folded into the length score (penalty outside expected min/max). Docs: `Scoring` (Length section) • Code: `src/length.rs::compute_length_consistency` + JSON/CSV emit in `src/main.rs::render_gene_record` |
| Conserved regions / alignment evidence | Implemented | Pillar `conserved_regions` (gated on ≥10 aligned seqs). Docs: `Scoring` • Code: `src/scoring.rs::compute_conserved_regions_score` + gating in `src/main.rs::build_scores_map` |
| Consistent start/end sequences (termini concordance) | Implemented | Docs: `Scoring` • Code: termini gating/scoring in `src/main.rs::build_scores_map` |
| Divergence vs homolog panel | Implemented | Docs: `Divergence & Phylo-Lite`, `Scoring` • Code: `src/scoring.rs::compute_divergence_score` |
| Different genes (fusion/split) | Implemented | Penalizes the per-gene score via `structvar_multiplier` when `FusionPossible` / `SplitPossible` is detected. Docs: `Scoring` • Code: `src/structvar.rs::analyze` + applied in `src/main.rs::build_scores_map` |
| Hard caps for catastrophic calls | Implemented (optional) | Optional `[scoring.caps]` can enforce ceilings for fusions/splits/duplications. Docs: `Scoring` • Code: `src/main.rs::ScoringCapsConfigOverride` + applied in `src/main.rs::build_scores_map` |
| Intrinsic sequence QC (ambiguity/low-complexity/homopolymer) | Implemented | Docs: `Scoring` • Code: `src/metrics.rs::compute_intrinsic` + `src/scoring.rs::compute_intrinsic_score` |
| ORF validation (nt) | Not yet | Code stub: `src/orf.rs` (placeholder translator); current ORF-ish heuristics are AA-only in `src/metrics.rs::compute_intrinsic` |
| “All hits align within a single ORF” (nt) | Not yet | Needs nucleotide + ORF-aware alignment stage |
| Genomic context (splice-site quality from genome+GFF3) | Implemented | Docs: `Genomic Context`, `Scoring` • Code: `src/genomic.rs::analyze_gff_context` + `src/scoring.rs::compute_genomic_score_with_cfg` |
| Phylogenetic consistency | Partial | Today’s proxy is taxonomy consensus + divergence (no full tree inference). Docs: `Taxonomy & Domains`, `Divergence & Phylo-Lite` |
| Pfam / domains architecture agreement | Implemented | Docs: `Taxonomy & Domains`, `Scoring` • Code: `src/hmmer.rs::domains_architecture_diagnostics` |
| Domain strength (top e-value) as part of final score | Implemented | Pillar `domains_strength` (weight key `domains_strength`). Score mapping: `src/scoring.rs::compute_domains_strength_score` and weighted in `src/main.rs::build_scores_map`. Emitted as `domains_score` in CSV and `score_components.domains_strength` in JSON (`src/main.rs::render_gene_record`). |
| Intron/exon length distributions | Not yet | Would require additional GFF/genome statistics beyond splice motifs |
| “Nearby genome for comparison” | Not yet | Needs genome neighborhood model/feature extraction |
| Number of transcripts | Not yet | Needs transcript grouping + scoring policy |
| Consistent exon/intron boundaries | Partial | Current proxy: splice motif quality only; deeper boundary consistency not implemented |
| Exon fuses/splits | Not yet | Likely ties into structural-variation and/or transcript evidence |
| % nucleotides/AA found via kmers | Not yet | Needs k-mer reference models and/or read evidence |
| Clusters to expand beyond top hits | Implemented | Prepare step produces `clusters.recluster`; analyze can backfill scarce panels from it. Docs: `Scoring` (Cluster backfill), `Prepare` • Code: `src/clusters.rs::load_cluster_map` + `src/consensus.rs::finalize_with_backfill` + provenance emit in `src/main.rs` |

Newer chapters that cover implemented ideas in detail:
- `Genomic Context`, `Taxonomy & Domains`, `Divergence & Phylo-Lite`, `Scoring`, `Observability`, `Performance & Tuning`.

### Pre-Processing
* Input Format Verification
* DIAMOND blast against reference or pre-filter using syncmers?
* Convert to SFASTA internally

### Analyses
* Length (cluster and rank validations) - distribution as well
* Coverage of consensus
* Conserved Regions - MAFFT Alignment - 10 most sig (following GV), find missing or extra regions
* Different Genes - Gene fusions and gene splits
* Open Reading Frames (ab initio and simliarity based validations, nt seqs only) - More than one major ORF
* ORFs - All BLAST hits to align within a single ORF
* Complexity (Similar to sdust) 
* Can we do number of exons / introns / distribution? log likelihood
* Inconsistent Insertions / Deletions (part of Conserved Regions)
* Consistent Start + End sequences
* GC Content? - Species Specific
* Phylogenetic Consistency?
* Intron Length Distributions - Should be species specific
* Exon Length Distributions - Should be species specific
* Nearby genome for comparison - What is the distribution there of exons/introns/length/etc?
* Number of transcripts?
* Consistent exon/intron boundaries
* Exon fuses and splits
* % Nucleotides found (via kmers)
* % AA's found (via kmers)
* Maybe use clusters? Diamond DeepClust? Looks like a good idea! Then can exceed the 10 most significant hits

### Ideas for Speed
* ECS dataflow
* Use SFASTA internally
* Use syncmers for pre-filtering?
* Use DIAMOND instead of BLAST
* Rust (ofc)

### Other Ideas
* Syncmers + BiWFA instead of diamond?
