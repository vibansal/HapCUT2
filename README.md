HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

## Important Announcement:
For simplicity, the switch confidence and SNP confidence scores in the last two columns of output are now being represented as phred-scaled probabilities of error, like standard quality scores. Note that the old way printed log10(1-P(error)) instead. In the new format, 0 represents low quality. 100 represents high-quality.

```
OLD       => NEW
-0.000434 => 30.00
```

Also the ```--threshold``` parameter is being interpreted in the same way, rather than as an unscaled floating-point probability:

```
OLD               => NEW
--threshold 0.999 => --threshold 30.0
```

## About:
HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed for speed and accuracy across diverse sequencing technologies, including but not limited to:
- NGS short reads (Illumina HiSeq)
- clone-based sequencing (Fosmid or BAC clones)
- SMRT reads (PacBio)
- 10X Genomics Linked-Reads
- proximity-ligation (Hi-C) reads
- high-coverage sequencing (>40x coverage-per-SNP) using above technologies
- combinations of the above technologies (e.g. scaffold long reads with Hi-C reads)

## Citation:
If you use HapCUT2 in your research, please cite:

[Edge, P., Bafna, V. & Bansal, V. HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies. Genome Res. gr.213462.116 (2016). doi:10.1101/gr.213462.116](http://genome.cshlp.org/content/early/2016/12/09/gr.213462.116.abstract)

## to build:

 ```make ```

The makefile will attempt to build samtools 1.2 and htslib 1.2.1 as git submodules.
If you already have samtools 1.2 and htslib 1.2.1 installed, you can optionally edit the SAMTOOLS and HTSLIB variables in the Makefile to point to the directories where they are installed, prior to building.

## to install:

```
sudo make install-hairs
sudo make install-hapcut2
```

## to uninstall:

```
sudo make uninstall-hairs
sudo make uninstall-hapcut2
```
## Input:
HapCUT2 requires the following input:
- BAM file for an individual containing reads aligned to a reference genome
- VCF file containing SNVs for the individual with respect to the reference

## To Run:

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

### Note about HAPCUT2 options
For the vast majority of use cases (including most short reads, long reads, clone sequences), only the required HAPCUT2 options above are necessary.
In the case of Hi-C reads, it is recommended to use ```--hic 1``` for both extractHAIRS and HAPCUT2.
Based on user preference, SNV pruning (filtering of low-quality phased SNVs) may be adjusted with ```--threshold <float>``` (closer to 1 prunes more, closer to 0.5 prunes less) or turned off with ```--no_prune 1```.

## Output Format:
Haplotype blocks are printed to the output file. For a given block, column 2 represents
the allele on one chromosome copy (0 for reference, 1 for variant), while column 3 represents
the allele on the other copy.

Each block starts with a block header with the following format:

BLOCK: offset: \<SNV offset\> len: \<SNV span of block\> phased: \<\# SNVs phased\> SPAN: \<base pair span of block\> fragments \<\# of fragments in block\>

Following the header, there is one line per SNV with the following tab-delimited fields:

1. VCF file index (1-based index of the line in the input VCF describing variant)
2. allele on haploid chromosome copy A (0 means reference allele, 1 means variant allele)
3. allele on haploid chromosome copy B (0 means reference allele, 1 means variant allele)
4. chromosome
5. position
6. reference allele (allele corresponding to 0 in column 2 or 3)
7. variant allele (allele corresponding to 1 in column 2 or 3)
8. VCF genotype field (unedited, directly from original VCF)
9. discrete pruning status (1 means pruned, 0 means phased)
10. switch quality: phred-scaled estimated probability that there is a switch error starting at this SNV (0 means switch error is likely, 100 means switch is unlikely)
11. mismatch quality: phred-scaled estimated probability that there is a mismatch [single SNV] error at this SNV (0 means SNV is low quality, 100 means SNV is high quality)

Field 9 describes the status of the SNV under the discrete SNV pruning method introduced by RefHap (an SNV is pruned if there are equal reads consistent and inconsistent with the phase), with the slight difference that reads are assigned to haplotypes using likelihoods in our implementation. Use the option "--discrete_pruning 1" to automatically prune SNPs ('- -' phasing) based on the value of this field.

The values in field 10 and field 11 are quality scores that range from 0 to 100 (less confident to more confident). Field 10 is useful for controlling switch errors, especially on data types such as error-prone SMRT reads. It is empty by default (".") unless switch error scores are computed using "--error_analysis_mode 1" (compute switch error confidence but leave blocks intact and all SNVs phased for manual pruning later).

Important note: flag "--split_blocks 1" (compute switch error confidence and automatically split blocks using the value of --split_threshold) is currently broken, for the time being use "--error_analysis_mode 1" and manually split using field 10.

Field 11 is useful for controlling mismatch (single SNV) haplotype errors, similarly to field 9. The default behavior of HapCUT2 is to prune individual SNVs for which this confidence is less than 6.98 (probability of error 0.2), as these are highly likely to be errors.

## Calculating Haplotype Statistics
The calculate_haplotype_statistics script in the utilities directory calculates haplotype error rates with respect to a reference haplotype, as well as completeness statistics such as N50 and AN50.

## Converting HapCUT2 output to VCF format
Nils Homer has developed a tool HapCutToVcf that will soon support converting HapCUT2-formatted haplotype blocks into VCF format. It will be included with the fgbio tool suite, available [here](https://github.com/fulcrumgenomics/fgbio).

## Reproducing the HapCUT2 manuscript

The directory **reproduce_hapcut2_paper** contains the source code and pipeline used to obtain the results of the HapCUT2 manuscript (linked above). It is nearly complete except for some early data access and cleaning steps which are not yet integrated into the pipeline, but these will be added soon.

## Example pipelines for various types of sequencing data

The directory **recipes** contains example pipelines to assemble haplotypes from various types of sequencing data.
