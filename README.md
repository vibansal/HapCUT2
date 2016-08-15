HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies
======

#About:
HapCUT2 is a tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
We found that previously described haplotype assembly methods are specialized for specific read technologies or protocols, with slow or inaccurate performance on others. With this in mind, HapCUT2 is designed to meet the following criteria:
- Support for diverse sequencing technologies, including but not limited to:
    * NGS short reads (Illumina HiSeq)
    * clone-based sequencing (Fosmid or BAC clones)
    * SMRT reads (PacBio)
    * proximity-ligation (Hi-C) reads
    * high-coverage sequencing (>40x coverage-per-SNP) using above technologies
    
- Although not yet tested, technologies similar in design to the above should work very well (synthetic long reads, dovetail Chicago method, etc.).
- Fast performance across all read technologies, with unprecendented read-depth scalability.
- Excellent accuracy on each sequencing technology (accuracy better than or equivalent to the previous best methods).
- Low memory footprint.
- Simple to use, no parameter tuning necessary in practice

##to build:

 ```make ```
 
The makefile will attempt to build samtools 1.2 and htslib 1.2.1 as git submodules. 
If you already have samtools 1.2 and htslib 1.2.1 installed, you can optionally edit the SAMTOOLS and HTSLIB variables in the Makefile to point to the directories where they are installed, prior to building.

##to install:

```sudo make install-hairs```

```sudo make install-hapcut2```

##to uninstall:

```sudo make uninstall-hairs```

```sudo make uninstall-hapcut2```

##to run:
The extractHAIRS tool converts a BAM file to a compact format (fragment file) containing only haplotype-informative information.
The HAPCUT2 tool assembles a fragment file created with extractHAIRS into haplotype blocks.
### run with no options 
###if not installed:
 ```./build/extractHAIRS ```

 ```./build/HAPCUT2 ```
###if installed:
 ```extractHAIRS ```

 ```HAPCUT2 ```
##Converting HapCUT2 output to VCF format
Nils Homer has developed a tool HapCutToVcf for converting HapCUT2-formatted haplotype blocks into VCF format. It is included with the fgbio tool suite, available [here](https://github.com/fulcrumgenomics/fgbio).

