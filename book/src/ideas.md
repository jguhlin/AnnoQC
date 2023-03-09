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

### Ideas for Speed
* ECS dataflow
* Use SFASTA internally
* Use syncmers for pre-filtering?
* Use DIAMOND instead of BLAST
* Rust (ofc)

### Other Ideas
* Syncmers + BiWFA instead of diamond?