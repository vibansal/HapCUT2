## Updates and Announcements:

#### March 24, 2020

New release of HapCUT2 with following updates:

1. uses fast logsum approximation for likelihood calculation that speeds up HapCUT2 by 7-8 times
2. fixes bugs in extracthairs for PacBio CCS/HiFi reads 
3. several other bug fixes
4. uses htslib for reading bam files
5. reorganization of hapcut2 and vcf parsing code

#### April 4, 2019

New release of HapCUT2 with updates to local realignment for long read allelotyping and ability to process specific genomic regions (--region option in extractHAIRS). 

#### March 13, 2018

HapCUT2 now outputs the phased variants to a VCF file ("haplotype_output_file".phased.vcf) in addition to the haplotype blocks. This VCF file preserves most of the information in the original VCF file but removes phasing information (if any) that is present in the input VCF. Standard VCF tags (PS, PQ, PD: https://samtools.github.io/hts-specs/VCFv4.1.pdf) are used to 
annotate the phased variants. This is a new feature and additional changes to the output VCF are likely. 

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
