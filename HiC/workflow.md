Sample pipeline to assemble haplotypes from Hi-C data with HapCUT2
======

## About:
This is an example pipeline for assembling haplotypes for a single
individual using Hi-C data. It is written in Snakemake,
a python-based workflow management tool (available using pip):

```
pip3 install snakemake
```
It is assumed that the Hi-C reads are in paired-end fastq format.
If the Hi-C reads are already in BAM format they can be placed in 'data/hic_processed.bam' to skip the earlier steps.
Before running the pipeline, edit the top of the Snakefile to specify the paths to the
required data and tools.
The pipeline will generate haplotype block files and place them in 'output/chr*.hap'

## Requirements:
- Snakemake 3.5.5
- paired-end fastq files (ending in *_1.fastq and *_2.fastq) containing Hi-C reads
- VCF files containing variants for the individual, named chr1.vcf, chr2.vcf...
- HapCUT2
- extractHAIRS
- BWA 0.7.12-r1044
- Samtools 1.3, htslib 1.3
- Picard Tools 2.8
- Bamtools 2.4

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
