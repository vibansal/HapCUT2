This directory contains code for FragmentCut, a method for the detection of long DNA fragments from dilution pool sequencing experiments.

#### Input for FragmentCut:

1. sorted BAM file for each pool
2. reference fasta file (indexed)
3. genome mappability file (fasta format)
4. VCF file with heterozygous variants


#### Output from FragmentCut:

1. bedfile with list of fragment intervals 
2. HapCUT formatted fragment file for haplotype assembly


#### NOTES 

1. FragmentCut works for 10X Genomics bam files as well (a bed file) 
