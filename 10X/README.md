

## pipeline for generating haplotype fragment file from 10x BAM file that can be used to assemble haplotypes using HapCUT2 or other tools 

1. extract intervals that correspond to boundaries of long molecules using barcode information from the 10X sorted BAM file. Each molecule corresponds to an interval and a barcode

        python getMolecules.py -d 20000 NA12878.chr20.bam NA12878.chr20.molecules.bed
    

2. For each molecule, extract reads that overlap variants using FragmentCut (FragmentCut source code is in a separate directory): 

        ./FragmentCut --bam NA12878.chr20.bam --VCF NA12878.chr20.hets.vcf --bed NA12878.chr20.molecules.bed --barcode 1 --out chr20.frags --ref reference.fa > chr20.log

3. run HapCUT2 using fragment file from step 2 

        ./HAPCUT2 --fragments chr20.frags --VCF NA12878.chr20.hets.vcf 
        
        
        
## NOTES 

1. This pipeline should be run on each chromosome separately for efficiency
2. The "-d" parameter in getMolecules.py is set to 20000 bases
