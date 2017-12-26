HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

## About:
HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed for speed and accuracy across diverse sequencing technologies, including but not limited to:
- NGS short reads (Illumina HiSeq)
- clone-based sequencing (Fosmid or BAC clones)
- SMRT reads (PacBio)
- Oxford Nanopore reads
- 10X Genomics Linked-Reads
- proximity-ligation (Hi-C) reads
- high-coverage sequencing (>40x coverage-per-SNP) using above technologies
- combinations of the above technologies (e.g. scaffold long reads with Hi-C reads)

See below for specific examples of command line options and best practices for some of these technologies.

NOTE: At this time HapCUT2 is for diploid organisms only. VCF input should contain diploid variants.

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
- VCF file containing **diploid** SNVs for the individual with respect to the reference

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

## 10X Genomics Linked-Reads
10X Genomics Linked Reads require an extra step to link short reads together into barcoded molecules:

(1) use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2.
```
./build/extractHAIRS --10X 1 --bam reads.sorted.bam --VCF variants.VCF --out unlinked_fragment_file
```
(2) use LinkFragments to link fragments into barcoded molecules:
```
python3 utilities/LinkFragments.py --bam reads.sorted.bam --VCF variants.VCF --fragments unlinked_fragment_file --out linked_fragment_file
```
(3) use HAPCUT2 to assemble fragment file into haplotype blocks.
```
./build/HAPCUT2 --nf 1 --fragments linked_fragment_file --vcf variantcalls.vcf --output haplotype_output_file
```

## Hi-C (Proximity Ligation) Sequencing Reads

For improved haplotype accuracy with Hi-C reads, use the --HiC 1 option for both extractHAIRS and HapCUT2 steps.

## Pacific Biosciences and Oxford Nanopore Sequencing Reads

Use the --pacbio 1 and --ont 1 options in extractHAIRS for greatly improved accuracy when using Pacific Biosciences and Oxford Nanopore reads, respectively. It is also highly recommended to split blocks based on the switch quality score, which can be computed using the --ea 1 option in HapCUT2.

Here is an example using Pacific Biosciences data (replace --pacbio with --ont for oxford nanopore):
```
./build/extractHAIRS --pacbio 1 --bam reads.sorted.bam --VCF variants.VCF --out fragment_file
./build/HAPCUT2 --ea 1 --fragments fragment_file --vcf variantcalls.vcf --output haplotype_output_file
python3 utilities/prune_haplotype.py -i haplotype_output_file -o haplotype_output_file.pruned --min_mismatch_qual 30 --min_switch_qual 30
# the quality-filtered haplotype is in haplotype_output_file.pruned
```
The --indels option may be used if desired -- the realignment strategy used with these options allows better detection of indel variants in fragments than the previous approach.
## Calculating Haplotype Statistics
The calculate_haplotype_statistics script in the utilities directory calculates haplotype error rates with respect to a reference haplotype, as well as completeness statistics such as N50 and AN50.

## Converting HapCUT2 output to VCF format
Nils Homer has developed a tool HapCutToVcf that supports converting HapCUT2-formatted haplotype blocks into VCF format. It will be included with the fgbio tool suite, available [here](https://github.com/fulcrumgenomics/fgbio).

## Reproducing the HapCUT2 manuscript

The directory **reproduce_hapcut2_paper** contains the source code and pipeline used to obtain the results of the HapCUT2 manuscript (linked above). It is nearly complete except for some early data access and cleaning steps which are not yet integrated into the pipeline, but these will be added soon.

## Example pipelines for various types of sequencing data

The directory **recipes** contains example pipelines to assemble haplotypes from various types of sequencing data.

## Updates and Announcements:

#### August 14, 2017
Extracthairs now has optimizations for error prone long read technologies (Pacific Biosciences and Oxford Nanopore). The strategy performs a sensitive realignment in a window around the potential variant. A read-window is aligned to both the reference sequence, and the variant sequence (reference sequence modified to contain the variant). An allele call is assigned based on the best alignment, and the "base quality" of the allele call is determined by a bayesian posterior calculated using both alignment scores.

In the case of a cluster of n variants that are nearby one another (distance to nearest variant < 20 bps), the alignments of read-windows to each variant might not be independent. In this case, all (2^n) possible haplotypes for those variants are enumerated and a read-window spanning the cluster is aligned to each haplotype sequence. The maximum scoring haplotype is selected to determine the allele call for each variant position in the cluster, and the "base quality score" of each position is derived from the sum of posterior probabilities of haplotypes that do not contain that allele in that position.

This also allows for more sensitive phasing of indel variants -- previously these were harder to detect from CIGAR string information alone, but the new approach allows better phasing of arbitrary insertions, deletions, MNPs, etc.

What this means for users: use the --pacbio 1 option with Pacific Biosciences reads, and use the --ont 1 option with Oxford Nanopore reads. See newly created section above for a full example execution.

(note: the default alignment parameters for pacbio and oxford nanopore are not yet finalized, but there will be SIGNIFICANT improvements in accuracy regardless, over the previous approach.)

#### July 24, 2017
The pipeline for phasing 10X Genomics linked reads has been updated. If you are currently using the old 10X pipeline, it is recommended to switch to the new one now. The new pipeline uses a new LinkFragments.py script to link haplotype fragments from short reads into long haplotype fragments based on their barcode. See the instructions above and the updated "10X" and "HiC + 10X" workflows in the recipes folder.

There are two significant differences between the LinkFragments pipeline and the old FragmentCut-based pipeline -- firstly, it circumvents a bug present in the FragmentCut code that resulted in loss of variants from 10X molecules. Secondly, it uses corrected BX barcodes rather than RX barcodes to link reads together which will result in more short reads being assigned to the correct molecule. On a practical level, the new approach meshes better with the production extractHAIRS code. Single-step 10X haplotype fragment generation may be integrated directly into the extractHAIRS tool in the future if there is demand.

#### May 4, 2017
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
