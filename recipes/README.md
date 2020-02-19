HapCUT2 Recipes: Sample pipelines to assemble haplotypes from various sequencing data types
======

## About:
This directory contains example pipelines for assembling haplotypes using
various types of sequencing data. Each pipeline is self-contained
except for dependencies on HapCUT2 and related software, various common
bioinformatics tools, and the input data.
Each pipeline is in a separate directory as described here:

- **HiC:** Hi-C (proximity ligation) reads
- **10X:** 10X genomics gemcode linked-reads
- **HiC_10X:** HiC + 10X combination
- **HiC_Longread:** HiC + Generic long read (e.g. PacBio) combination


### Processing Hi-C reads (updated 02/18/2020)

BWA-mem can map Hi-C read pairs using the command "bwa mem -5SPM ref.fa read1.fq read2.fq" that preserves mate-pair information (thanks to @conchoecia)
for highlighting this. There is no need for further processing except for marking PCR duplicates


## how to combine 10X + long read data for phasing

```
./build/extractHAIRS --10X 1 --bam tenX.sorted.bam --VCF variants.VCF --out unlinked_fragment_file
python3 utilities/LinkFragments.py --bam tenX.sorted.bam --VCF variants.VCF --fragments unlinked_fragment_file --out linked_fragment_file
./build/extractHAIRS --pacbio 1 --nf 1 --bam pacbio.sorted.bam --VCF variants.VCF --out pacbio_fragment_file
cat linked_fragment_file pacbio_fragment_file > combined_fragment_file
./build/HAPCUT2 --nf 1 --fragments combined_fragment_file --VCF variants.vcf --output haplotype_output_file
```
