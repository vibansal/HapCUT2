Sample pipeline to assemble haplotypes from 10X Genomics data with HapCUT2
======

## About:
This is an example pipeline for assembling haplotypes using
10X Genomics data for a single individual. It is written in Snakemake,
a python-based workflow management tool (available using pip):

```
pip3 install snakemake
```
It is assumed that the 10X reads are in aligned bam format (It is recommended to align
10X reads using e.g. 10X's own LongRanger software suite.)
Before running the pipeline, edit the top of the Snakefile to specify the paths to the
required data and tools.
The pipeline will generate haplotype block files and place them in 'output/chr*.hap'

## Requirements:
- Snakemake 3.5.5
- a reference genome in fasta format
- a bam file containing aligned, barcoded 10X reads (BX tag)
- VCF files containing variants for the individual, named chr1.vcf, chr2.vcf...
- HapCUT2 and extractHAIRS
- pysam 0.11.2
- Samtools 1.3, htslib 1.3

For instructions to install HapCUT2 and extractHAIRS, see the README at the root of this project.

## Running
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
