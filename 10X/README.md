

### pipeline for generating haplotype fragment file from 10x BAM file that can be used to assemble haplotypes using HapCUT2 or other tools 

0. Compile FragmentCut code by running 'make all' in the FragmentCut source directory 

1. extract intervals (4-tuple = (chr,start,end,barcode) ) that correspond to boundaries of long molecules using barcode information from the 10X sorted BAM file.

        python getMolecules.py -d 20000 NA12878.chr20.bam NA12878.chr20.molecules.bed
    
2. For each molecule, extract reads that overlap variants using FragmentCut: 

        ./FragmentCut --bam NA12878.chr20.bam --VCF NA12878.chr20.hets.vcf --bed NA12878.chr20.molecules.bed --barcode 1 --out chr20.frags --ref reference.fa > chr20.log

3. run HapCUT2 using fragment file 'chr20.frags' (for instructions, see the HapCUT2 README file) 

	./build/HAPCUT2 --fragments chr20.frags --vcf NA12878.chr20.hets.vcf  --output haplotype_output_file 

               
        
### NOTES 

1. This pipeline should be run on each chromosome separately 
2. The "-d" parameter in getMolecules.py defines the distance between adjacent molecules. Increasing this parameter can result in adjacent molecules (fragments) being merged into a single molecule. 
