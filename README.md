HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

##About:
HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed to meet the following criteria:
- Support for diverse sequencing technologies, including but not limited to:
    * NGS short reads (Illumina HiSeq)
    * clone-based sequencing (Fosmid or BAC clones)
    * SMRT reads (PacBio)
    * proximity-ligation (Hi-C) reads
    * high-coverage sequencing (>40x coverage-per-SNP) using above technologies
    
- Although not yet tested, technologies similar in design to the above should work very well (synthetic long reads, dovetail Chicago method, etc.).
- Fast performance across all read technologies, with unprecendented read-coverage scalability.
- Excellent accuracy on each sequencing technology (accuracy better than or equivalent to the previous best methods).
- Low memory footprint.
- Simple to use, no parameter tuning necessary in practice

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

##Running extractHAIRS:
The extractHAIRS tool converts a BAM file to a compact format (fragment file) containing only haplotype-informative information. This is a necessary precursor step to running HapCUT2. To run extractHAIRs and see options:

 ```./build/extractHAIRS ``` or  ```extractHAIRS ``` if installed.

```
Extract haplotype informative reads (HAIRS) from coordinate sorted BAM files 

./extractHAIRS [options] --bam reads.sorted.bam --VCF variants.VCF --out fragment_file 

Options:
--qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33 
--mbq <INT> : minimum base quality to consider a base for haplotype fragment, default 13
--mmq <INT> : minimum read mapping quality to consider a read for phasing, default 20
--hic <0/1> : sets default maxIS to 40MB, prints matrix in new HiC format
--new_format, --nf <0/1> : prints matrix in new format. Requires --new_format option when running HapCUT2.
--VCF <FILENAME> : variant file with genotypes for a single individual in VCF format
--variants : variant file in hapCUT format (use this option or --VCF option but not both), this option will be phased out in future releases
--maxIS <INT> : maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000
--minIS <INT> : minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0
--PEonly <0/1> : do not use single end reads, default is 0 (use all reads)
--indels <0/1> : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option
--noquality <INTEGER> : if the bam file does not have quality string, this value will be used as the uniform quality value, default 0 
--ref <FILENAME> : reference sequence file (in fasta format), optional but required for indels, should be indexed using samtools
--out <FILENAME> : output filename for haplotype fragments, if not provided, fragments will be output to stdout
```

##Running HAPCUT2:
The HAPCUT2 tool assembles a fragment file created with extractHAIRS into haplotype blocks. To run HapCUT2 and see options:

 ```./build/HAPCUT2 ``` or  ```HAPCUT2 ``` if installed.

###Note about HAPCUT2 options
For the vast majority of use cases (including most short reads, long reads, clone sequences), only the required options are necessary: ```--fragments```, ```--vcf```, and ```--output```.
In the case of Hi-C reads, it is recommended to use ```--hic 1``` for both extractHAIRS and HAPCUT2.
Based on user preference, pruning may be adjusted with ```--threshold``` or turned off with ```--no_prune```.

```
HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies

USAGE : ./HAPCUT2 --fragments fragment_file --vcf variantcalls.vcf --output haplotype_output_file

Basic Options:
--fragments, --f <FILENAME>:        file with haplotype-informative reads generated using the extracthairs program
--vcf <FILENAME>:                   variant file in VCF format (use EXACT SAME file that was used for the extracthairs program)
--output, --o <FILENAME> :          file to which phased haplotype segments/blocks will be output
--converge, --c <int>:              cut off iterations (global or maxcut) after this many iterations with no improvement. default: 5
--verbose, --v <0/1>:               verbose mode: print extra information to stdout and stderr. default: 0

Read Technology Options:
--hic <0/1> :                       increases accuracy on Hi-C data; models h-trans errors directly from the data. default: 0
--hic_htrans_file, --hf <FILENAME>  optional tab-delimited input file where second column specifies h-trans error probabilities for insert size bins 0-50Kb, 50Kb-100Kb, etc.
--qv_offset, --qo <33/48/64> :      quality value offset for base quality scores, default: 33 (use same value as for extracthairs)
--long_reads, --lr <0/1> :          reduces memory when phasing long read data with many SNPs per read. default: automatic.

Haplotype Post-Processing Options:
--threshold, --t <float>:           threshold for pruning low-confidence SNPs (closer to 1 prunes more, closer to 0.5 prunes less). default: 0.8
--skip_prune, --sp <0/1>:           skip default likelihood pruning step (prune SNPs after the fact using column 11 of the output). default: 0
--split_blocks, --sb <0/1>:         split blocks using simple likelihood score to reduce switch errors. default: 0
--split_threshold, --st <float>:    threshold for splitting blocks (closer to 1 splits more, closer to 0.5 splits less). default: 0.8
--call_homozygous, --ch <0/1>:      call positions as homozygous if they appear to be false heterozygotes. default: 0
--discrete_pruning, --dp <0/1>:     use discrete heuristic to prune SNPs. default: 0
--error_analysis_mode, --ea <0/1>:  compute switch confidence scores and print to haplotype file but don't split blocks or prune. default: 0

Advanced Options:
--new_format, --nf <0/1>:           use new Hi-C fragment matrix file format (but don't do h-trans error modeling). default: 0
--max_iter, --mi <int> :            maximum number of global iterations. Preferable to tweak --converge option instead. default: 10000
--maxcut_iter, --mc <int> :         maximum number of max-likelihood-cut iterations. Preferable to tweak --converge option instead. default: 10000
--htrans_read_lowbound, --hrl <int> with --hic on, h-trans probability estimation will require this many matepairs per window. default: 500
--htrans_max_window, --hmw <int>    with --hic on, the insert-size window for h-trans probability estimation will not expand larger than this many basepairs. default: 4000000


Hi-C-specific Notes:
  (1) When running extractHAIRS, must use --hic 1 option to create a fragment matrix in the new Hi-C format.
  (2) When running HapCUT2, use --hic 1 if h-trans probabilities are unknown. Use --hic_htrans_file if they are known
  (3) Using --hic_htrans_file is faster than --hic and may yield better results at low read coverage (>30x).
  (4) Set --converge to a larger value if possible/reasonable.
```

##Converting HapCUT2 output to VCF format
Nils Homer has developed a tool HapCutToVcf for converting HapCUT2-formatted haplotype blocks into VCF format. It is included with the fgbio tool suite, available [here](https://github.com/fulcrumgenomics/fgbio).

