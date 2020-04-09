HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

## About:
HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed for speed and accuracy across diverse sequencing technologies, including but not limited to:
- NGS short reads (Illumina HiSeq)
- single-molecule long reads (PacBio and Oxford Nanopore)
- Linked-Reads (e.g. 10X Genomics, stLFR or TELL-seq)
- proximity-ligation (Hi-C) reads
- high-coverage sequencing (>40x coverage-per-SNP) using above technologies
- combinations of the above technologies (e.g. scaffold long reads with Hi-C reads)

See below for specific examples of command line options and best practices for some of these technologies.

NOTE: At this time HapCUT2 is for diploid organisms only and can assemble haplotypes for one individual at a time. VCF input should contain variants and genotypes for a single diploid individual. 

## Citation:
If you use HapCUT2 in your research, please cite:

[Edge, P., Bafna, V. & Bansal, V. HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies. Genome Res. gr.213462.116 (2016). doi:10.1101/gr.213462.116](http://genome.cshlp.org/content/early/2016/12/09/gr.213462.116.abstract)

## dependencies:
Requires htslib > 1.2.1. It is assumed that htslib is installed, but otherwise the path can be specified in the Makefile (HTSLIB=path). 
To install htslib directly from source: git clone https://github.com/samtools/htslib.git

## to build:

 ```make ```

## Input:
HapCUT2 requires the following input:
- BAM file for an individual containing reads aligned to a reference genome
- VCF file containing short variant calls (SNVs and indels) and **diploid** genotypes for the same individual with respect to the reference genome

**Note: the program does not accept gzipped VCF files**

### To Run:

Assembling haplotypes requires two steps:

(1) use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2. 
```
./build/extractHAIRS [options] --bam reads.sorted.bam --VCF variants.vcf --out fragment_file
```

(2) use HAPCUT2 to assemble fragment file into haplotype blocks.
```
./build/HAPCUT2 --fragments fragment_file --VCF variants.vcf --output haplotype_output_file
```

If you have data from different technologies or in different bam files for the same individual, run step (1) separately on each input bam file and combine the output fragment files into a single file that can be used as input in step (2). 

Run the programs without arguments to see all options.

### Note about HAPCUT2 options
For the vast majority of use cases (including most short reads, long reads, clone sequences), only the required HAPCUT2 options above are necessary.
Based on user preference, SNV pruning (filtering of low-quality phased SNVs) may be adjusted with ```--threshold <float>``` (closer to 1 prunes more, closer to 0.5 prunes less) or turned off with ```--no_prune 1```.

### Output Format:

There are two output files with the phased variants: 

1. A phased block output file. The format of this file is described [here](outputformat.md)

2. A phased VCF file "output_haplotype_file.phased.vcf". The format of this file follows the standard VCF format. This is a recent addition.


## Pacific Biosciences and Oxford Nanopore Reads

Use the --pacbio 1 and --ont 1 options in extractHAIRS for greatly improved accuracy when using Pacific Biosciences and Oxford Nanopore reads, respectively. 
Here is an example using Pacific Biosciences data (replace --pacbio with --ont for oxford nanopore):
```
./build/extractHAIRS --pacbio 1 --bam reads.sorted.bam --VCF variants.VCF --out fragment_file
./build/HAPCUT2 --ea 1 --fragments fragment_file --VCF variantcalls.vcf --output haplotype_output_file
python3 utilities/prune_haplotype.py -i haplotype_output_file -o haplotype_output_file.pruned --min_mismatch_qual 30 --min_switch_qual 30
# the quality-filtered haplotype is in haplotype_output_file.pruned
```
The --indels option may be used if desired -- the realignment strategy used with these options allows better detection of indel variants in fragments than the previous approach.

## Linked Reads (10X Genomics, stLFR, etc)

Phasing using Linked reads require an extra step to link the short reads together into barcoded molecules. Details are provided [here](linkedreads.md)

## Hi-C (Proximity Ligation) Sequencing Reads

For improved haplotype accuracy with Hi-C reads, use the --HiC 1 option for both extractHAIRS and HapCUT2 steps.

## Calculating Haplotype Statistics
The calculate_haplotype_statistics script in the utilities directory calculates haplotype error rates with respect to a reference haplotype, as well as completeness statistics such as N50 and AN50.

## Example pipelines for various types of sequencing data

The directory **recipes** contains example pipelines to assemble haplotypes from various types of sequencing data.

## [Updates to the code](UPDATES.md)

