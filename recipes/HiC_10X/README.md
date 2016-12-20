Sample pipeline to assemble combined Hi-C + 10X Genomics data with HapCUT2
======

##About:
This is an example pipeline for assembling haplotypes using combined Hi-C
and 10X Genomics data for the same individual. It is written in Snakemake,
a python-based workflow management tool (available using pip):

```
pip3 install snakemake
```
It is assumed that the Hi-C reads are in paired-end fastq format and the
10X reads are in aligned bam format (It is recommended to align
10X reads using e.g. 10X's own LongRanger software suite.)
Before running the pipeline, edit the top of the Snakefile to specify the paths to the
required data and tools.
If the Hi-C reads are already in BAM format they can be placed in
'data/hic_processed.bam' to skip the earlier steps. The pipeline will generate 
haplotype block files and place them in 'output/chr*.hap'

##Caveats:
This pipeline was extracted and simplified from the pipeline used for the HapCUT2
manuscript. It has not been tested extensively since then, so please report any problems
encountered.

##Requirements:
- Snakemake 3.5.5
- a reference genome in fasta format
- a bam file containing aligned, barcoded 10X reads
- paired-end fastq files (ending in *_1.fastq and *_2.fastq) containing Hi-C reads
- VCF files containing variants for the individual, named chr1.vcf, chr2.vcf...
- HapCUT2
- extractHAIRS
- [https://github.com/vibansal/HapCUT2/tree/master/FragmentCut](FragmentCut)
- BWA 0.7.12-r1044
- Samtools 1-2-244, htslib 1.21
- Picard Tools 2.8
- Bamtools 2.4

##Running
To do a dry run and list the executions steps, call:
```
snakemake -n
```
To run the pipeline, call:
```
snakemake
```
It is recommended to specify the number of available cores when calling snakemake
(--cores), or run in cluster mode if applicable.
See the snakemake documentation for more details.
