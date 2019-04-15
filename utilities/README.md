Haplotyping Utilities
======

## calculate_haplotype_statistics.py:

This file is a self-contained script to compute a number of useful statistics on haplotypes
assembled using HapCUT2 or similar tools. It is written in Python 3. To see options run with no arguments:
```
usage: calculate_haplotype_statistics.py [-h] [-v1 VCF1 [VCF1 ...]]
                                         [-v2 VCF2 [VCF2 ...]]
                                         [-h1 HAPLOTYPE_BLOCKS1 [HAPLOTYPE_BLOCKS1 ...]]
                                         [-h2 HAPLOTYPE_BLOCKS2 [HAPLOTYPE_BLOCKS2 ...]]
                                         [-i]

Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.
Error rates for an assembled haplotype (specified by -v1 and optionally -h1
arguments) are computed with respect to a "reference" haplotype (specified by
-v2 and optionally -h2 arguments). All files must contain information for one
chromosome only! To compute aggregate statistics across multiple chromosomes,
provide files for each chromosome/contig as an ordered list, using the same
chromosome order between flags. Note: Triallelic variants are supported, but
variants with more than 2 alternative alleles are currently NOT supported.
These variants are ignored. Also, variants where the ref and alt alleles
differ between the test haplotype and reference haplotype are skipped.

optional arguments:
  -h, --help            show this help message and exit
  -v1 VCF1 [VCF1 ...], --vcf1 VCF1 [VCF1 ...]
                        A phased, single sample VCF file to compute haplotype
                        statistics on.
  -v2 VCF2 [VCF2 ...], --vcf2 VCF2 [VCF2 ...]
                        A phased, single sample VCF file to use as the "ground
                        truth" haplotype.
  -h1 HAPLOTYPE_BLOCKS1 [HAPLOTYPE_BLOCKS1 ...], --haplotype_blocks1 HAPLOTYPE_BLOCKS1 [HAPLOTYPE_BLOCKS1 ...]
                        Override the haplotype information in "-v1" with the
                        information in this HapCUT2-format haplotype block
                        file. If this option is used, then the VCF specified
                        with -v1 MUST be the same VCF used with HapCUT2
                        (--vcf) to produce the haplotype block file!
  -h2 HAPLOTYPE_BLOCKS2 [HAPLOTYPE_BLOCKS2 ...], --haplotype_blocks2 HAPLOTYPE_BLOCKS2 [HAPLOTYPE_BLOCKS2 ...]
                        Override the haplotype information in "-v2" with the
                        information in this HapCUT2-format haplotype block
                        file. If this option is used, then the VCF specified
                        with -v2 MUST be the same VCF used with HapCUT2
                        (--vcf) to produce the haplotype block file!
  -i, --indels          Use this flag to consider indel variants. Default:
                        Indels ignored.
```

The assembled/test haplotype you wish to assess for accuracy should be input using the -v1 option.
If only -v1 is used, it is assumed to be a phased VCF. The "phase set" information (PS flag)
about haplotype blocks MUST be present.
If the -h1 option is used, then a haplotype block file (i.e. HapCUT2 output format)
can be specified to override the phase information in -v1. If this option is used,
then only the genotype information in -v1 will be considered. If -h1 is used,
then the vcf provided with -v1 should be the same VCF used to run HapCUT2 and produce
the haplotype block file.

The haplotype you wish to compare against (the "ground truth" haplotype, for computing error rates)
is specified in the exact same fashion, except specified with -v2 (and optionally -h2) options.

In the case of -v1 and -v2 inputs, only records for a single contig should be in a given file.
If you wish to assess multiple contigs (or the whole genome), every option (with the exception of -c) may be specified as a list as so:

```
python3 calculate_haplotype_statistics.py -h1 chr1.block chr2.block chr3.block -v1 chr1.vcf chr2.vcf chr3.vcf ........
```

If lists of files are provided, the contig order MUST be the same between arguments.

Here is a brief description of the output:
```
switch rate:          switch errors as a fraction of possible positions for switch errors
mismatch rate:        mismatch errors as a fraction of possible positions for mismatch errors
flat rate:            flat errors as a fraction of possible positions for flat errors
phased count:         count of total SNVs phased in the test haplotype
AN50:                 the AN50 metric of haplotype completeness
N50:                  the N50 metric of haplotype completeness
num snps max blk:   the fraction of SNVs in the largest (most variants phased) block
```

A switch error is defined as a position where the phase is switched from that of the previous heterozygous SNV, when compared to
the reference haplotype. Two switch errors in a row are instead counted as a mismatch error, a single position where the phase differs from the reference haplotype.

We use "flat error" to the minimum hamming distance between the two assembled haplotypes (for a given block)
and the reference haplotype. This is an alternative metric to observing switch/mismatch errors in tandem.
In general, this metric is thought to penalize switch errors too harshly.
It may be of interest for a dataset with extremely low incidence of long switch errors.

#### Update 11/27/2018:
The behavior of the N50 calculation has been changed. The N50 is now calculated with respect to the phased portion of the genome instead of the entire genome.
