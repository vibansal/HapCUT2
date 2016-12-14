This directory contains code for FragmentCut, a method for the detection of long DNA fragments from dilution pool sequencing experiments.

Input for FragmentCut:

sorted BAM file for each pool
reference fasta file (indexed)
genome mappability file (fasta format)
VCF file with heterozygous variants
Output from FragmentCut:

bedfile with list of fragment intervals
HapCUT formatted fragment file for haplotype assembly



## FragmentCut works for 10X Genomics bam files as well (additional input is a bed file) 
