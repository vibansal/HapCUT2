HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

##About:
HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed for speed and accuracy across diverse sequencing technologies, including but not limited to:
- NGS short reads (Illumina HiSeq)
- clone-based sequencing (Fosmid or BAC clones)
- SMRT reads (PacBio)
- proximity-ligation (Hi-C) reads
- high-coverage sequencing (>40x coverage-per-SNP) using above technologies

##to build:

 ```make ```
 
The makefile will attempt to build samtools 1.2 and htslib 1.2.1 as git submodules. 
If you already have samtools 1.2 and htslib 1.2.1 installed, you can optionally edit the SAMTOOLS and HTSLIB variables in the Makefile to point to the directories where they are installed, prior to building.

##to install:

```
sudo make install-hairs
sudo make install-hapcut2
```

##to uninstall:

```
sudo make uninstall-hairs
sudo make uninstall-hapcut2
```
##Input:
HapCUT2 requires the following input:
- BAM file for an individual containing reads aligned to a reference genome
- VCF file containing SNVs for the individual with respect to the reference

##To Run:

Assembling haplotypes requires two steps:

(1) use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2.
```
./build/extractHAIRS [options] --bam reads.sorted.bam --VCF variants.VCF --out fragment_file 
```

(2) use HAPCUT2 to assemble fragment file into haplotype blocks.
```
./build/HAPCUT2 --fragments fragment_file --vcf variantcalls.vcf --output haplotype_output_file
```

Run the programs without arguments to see all options.

###Note about HAPCUT2 options
For the vast majority of use cases (including most short reads, long reads, clone sequences), only the required HAPCUT2 options above are necessary.
In the case of Hi-C reads, it is recommended to use ```--hic 1``` for both extractHAIRS and HAPCUT2.
Based on user preference, SNV pruning (filtering of low-quality phased SNVs) may be adjusted with ```--threshold <float>``` (closer to 1 prunes more, closer to 0.5 prunes less) or turned off with ```--no_prune 1```.

##Output Format:
Haplotype blocks are printed to the output file, each with a block header with the following format:

BLOCK: offset: \<SNV offset\> len: \<SNV span of block\> phased: \<\# SNVs phased\> SPAN: \<base pair span of block\> fragments \<\# of fragments in block\>

Following the header, there is one line per SNV with the following tab-delimited fields:

1. VCF file index (1-based)
2. phased allele 1
3. phased allele 2
4. chromosome (from VCF)
5. position (from VCF)
6. allele 1 from VCF genotype
7. allele 2 from VCF genotype
8. VCF genotype field
9. discrete pruning status (1 == pruned, 0 == phased)
10. log<sub>10</sub>(confidence that there is not a switch error starting at this SNV)
11. log<sub>10</sub>(confidence that there is not a mismatch [single SNV] error at this SNV)

Field 9 describes the status of the SNV under the discrete SNV pruning method introduced by RefHap (an SNV is pruned if there are equal reads consistent and inconsistent with the phase), with the slight difference that reads are assigned to haplotypes using likelihoods in our implementation. Use the option "--discrete_pruning 1" to automatically prune SNPs ('- -' phasing) based on the value of this field.

The values in field 10 and field 11 are the log<sub>10</sub> of a confidence score that ranges from 0.5 to 1 (less confident to more confident). Field 10 is useful for controlling switch errors, especially on data types such as error-prone SMRT reads. It is fixed at 0 (confidence 1.0) unless switch error scores are computed using flags "--split_blocks 1" (compute switch error confidence and automatically split blocks using the value of --split_threshold) or "--error_analysis_mode 1" (compute switch error confidence but leave blocks intact and all SNVs phased for manual pruning later).

Field 11 is useful for controlling mismatch (single SNV) haplotype errors, similarly to field 9. The default behavior of HapCUT2 is to prune individual SNVs for which 10<sup>(field 11)</sup> < 0.8, as these are highly likely to be errors.

##Converting HapCUT2 output to VCF format
Nils Homer has developed a tool HapCutToVcf that will soon support converting HapCUT2-formatted haplotype blocks into VCF format. It will be included with the fgbio tool suite, available [here](https://github.com/fulcrumgenomics/fgbio).

##Coming Soon
- A Snakemake pipeline to reproduce the results of the HapCUT2 manuscript (in preparation).
- Hi-C pre-processing scripts and a complete Hi-C haplotyping pipeline

