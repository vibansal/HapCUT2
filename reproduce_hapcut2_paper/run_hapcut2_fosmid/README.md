Run HapCUT2 on Duitama et al fosmid data
======

This directory contains a bash script for running HapCUT2 on the fosmid data from
the 2012 Duitama et al study [1]:

Steps in the script:
1. download 1000 Genomes VCFs (phase 1, hg18) for NA12878 trio [2]
2. filter the VCF for only NA12878 heterozygous sites (and separate by chromosome)
3. check that variant indices and genomic coordinates match between filtered VCF and Duitama phased haplotype data
4. download fosmid fragment files and remove first line (matrix dimensions) which is incompatible with HapCUT2
5. run HapCUT2 on each fragment file to produce HapCUT2 haplotype files
6. use calculate_haplotype_statistics.py script to calculate the haplotype accuracy using the 1000G trio-phased VCF as ground truth

More detailed explanation of steps 1-3:
HapCUT2 requires a VCF file along with the fragment data, such that the variant indices
in the fragment file match those in the VCF. The original filtered VCF from 1000 Genomes
used to generate the fosmid fragments is not provided, so it is necessary to produce one. There is also a
need for ground-truth haplotypes to compare against, and the phased 1000g variants
can also be used for this purpose. So, the bash script
downloads VCF files for NA12878 from 1000 Genomes project phase 1, and filters
them manually (for NA12878 heterozygous sites only) and then checks
that the genomic coordinates and variant indices match against refhap-based phased
haplotype data files provided by Duitama et al.

## Steps to run
1. on a linux machine (tested using ubuntu) clone repository and build HapCUT2 using the makefile and instructions in main github README
2. change to this directory (HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid)
3. run run_hapcut2_fosmid_data.sh using bash:

```
bash run_hapcut2_fosmid_data.sh
```
or
```
chmod +x run_hapcut2_fosmid_data.sh
./run_hapcut2_fosmid_data.sh
```

The HapCUT2 haplotypes will be in ```data/hapcut2_haplotypes``` and the error rates compared to 1000G trio haplotypes will be printed to console.

References:

[1] Duitama, J., McEwen, G.K., Huebsch, T., Palczewski, S., Schulz, S., Verstrepen, K., Suk, E.K. and Hoehe, M.R., 2011. Fosmid-based whole genome haplotyping of a HapMap trio child: evaluation of Single Individual Haplotyping techniques. Nucleic acids research, 40(5), pp.2041-2053.

[2] 1000 Genomes Project Consortium, 2010. A map of human genome variation from population-scale sequencing. Nature, 467(7319), p.1061.
