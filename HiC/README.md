## Subdirectory for Hi-C (Proximity Ligation) data analysis, issues and pipelines 

For improved haplotype accuracy with Hi-C reads, use the --HiC 1 option for both extractHAIRS and HapCUT2 steps.

HapCUT2 implements an EM-like procedure to estimate the trans-error rates for Hi-C data. The default value of the maximum insert size is 40 megabases. 

## Pipeline for using population reference panels to improve accuracy and resolution of Hi-C haplotyping 

see this webpage: https://github.com/vibansal/IntegratedPhasing

## NA12878 HiC BAM files

HiC data for NA12878 from the Rao et al. (Cell 2014) paper is available for running HapCUT2. Currently, only files for chromosome 20 are available from a Dropbox shared folder (send email to vibansal [AT] ucsd dot edu for the link)
