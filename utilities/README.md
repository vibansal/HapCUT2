Haplotyping Utilities
======

## calculate_haplotype_statistics.py:

This file is a self-contained script to compute a number of useful statistics on haplotypes
assembled using HapCUT2 or similar tools. It is written in Python 3. To see options run with no arguments:
```
usage: calculate_haplotype_statistics.py [-h]
                                         [-h1 HAPLOTYPE_BLOCKS [HAPLOTYPE_BLOCKS ...]]
                                         [-v1 VCF [VCF ...]]
                                         [-f1 FRAGMENTS [FRAGMENTS ...]]
                                         [-pv [PHASED_VCF [PHASED_VCF ...]]]
                                         [-h2 REFERENCE_HAPLOTYPE_BLOCKS [REFERENCE_HAPLOTYPE_BLOCKS ...]]
                                         [-v2 [REFERENCE_VCF [REFERENCE_VCF ...]]]
                                         [-c [CONTIG_SIZE_FILE]]

Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.
Error rates for an assembled haplotype (specified by -h1,-v1,-f1 arguments)
are computed with respect to a "reference" haplotype (specified by -h2, -v2
arguments or -pv argument). All files must contain information for one
chromosome only (except --contig_size_file)! To compute aggregate statistics
across multiple chromosomes, provide files for each chromosome/contig as an
ordered list, using the same chromosome order between flags.

optional arguments:
  -h, --help            show this help message and exit
  -h1 HAPLOTYPE_BLOCKS [HAPLOTYPE_BLOCKS ...], --haplotype_blocks HAPLOTYPE_BLOCKS [HAPLOTYPE_BLOCKS ...]
                        haplotype block file(s) to compute statistics on
  -v1 VCF [VCF ...], --vcf VCF [VCF ...]
                        VCF file(s) that was used to generate h1 haplotype
                        fragments and phase h1 haplotype (--vcf in
                        extractHAIRS and HapCUT2)
  -f1 FRAGMENTS [FRAGMENTS ...], --fragments FRAGMENTS [FRAGMENTS ...]
                        HapCUT2 format fragment file(s) used to generate input
                        haplotype block file (-h1)
  -pv [PHASED_VCF [PHASED_VCF ...]], --phased_vcf [PHASED_VCF [PHASED_VCF ...]]
                        compute errors with respect to this phased single-
                        individual VCF file(s). NOTE: Files must be separated
                        by contig/chromosome! (Use with no arguments to use
                        same VCF(s) from --vcf.)
  -h2 REFERENCE_HAPLOTYPE_BLOCKS [REFERENCE_HAPLOTYPE_BLOCKS ...], --reference_haplotype_blocks REFERENCE_HAPLOTYPE_BLOCKS [REFERENCE_HAPLOTYPE_BLOCKS ...]
                        compute errors with respect to this haplotype block
                        file(s)
  -v2 [REFERENCE_VCF [REFERENCE_VCF ...]], --reference_vcf [REFERENCE_VCF [REFERENCE_VCF ...]]
                        VCF file(s) that was used to generate h2 haplotype
                        fragments and phase h2 haplotype (--vcf in
                        extractHAIRS and HapCUT2). Use with no arguments to
                        use same VCF(s) from --vcf.
  -c [CONTIG_SIZE_FILE], --contig_size_file [CONTIG_SIZE_FILE]
                        Tab-delimited file with size of contigs
                        (<contig>\t<size>). If not provided, N50 will not be
                        calculated.
```

The assembled haplotype you wish to assess should be input using the -h1, -v1, and -f1 options (the output haplotype block file
from HapCUT2, the VCF  used as input to HapCUT2, and haplotype fragment file used as input to HapCUT2, respectively).

The haplotype you wish to compare against (for computing error rates) can be either a VCF with phase information (-pv),
or another HapCUT2-type assembly (requires both -h2 and -v2 options).

The -c option is required if you wish to compute the N50 statistic. For instance, if Hg19 is your reference, you can use hg19.chrom.sizes (downloaded from UCSC and included here for convenience).

Since fragment files and haplotype block files are separated by chromosome, you may wish to compute aggregate statistics across
multiple contigs (or the whole genome). Every option (with the exception of -c) may be specified as a list as so:

```
python3 calculate_haplotype_statistics.py -h1 chr1.block chr2.block chr3.block -v1 chr1.vcf chr2.vcf chr3.vcf ........
```

If lists of files are provided, the contig order MUST be the same between arguments.

Here is a brief description of the output:
```
switch rate:          switch errors as a fraction of possible positions for switch errors
mismatch rate:        mismatch errors as a fraction of possible positions for mismatch errors
flat rate:            flat errors as a fraction of possible positions for flat errors
missing rate:         fraction of positions that were covered by an informative SNP, but were removed from the final haplotype
phased count:         count of total SNVs phased
AN50:                 the AN50 metric of haplotype completeness
N50:                  the N50 metric of haplotype completeness
max block snp frac:   the fraction of SNVs in the largest (most variants phased) block
```

A switch error is defined as a position where the phase is switched from that of the previous heterozygous SNV, when compared to
the reference haplotype. Two switch errors in a row are instead counted as a mismatch error, a single position where the phase differs from the reference haplotype.

We use "flat error" to the minimum hamming distance between the two assembled haplotypes (for a given block)
and the reference haplotype. This is an alternative metric to observing switch/mismatch errors in tandem.
In general, this metric is thought to penalize switch errors too harshly.
It may be of interest for a dataset with extremely low incidence of long switch errors.
