

## pipeline for getting haplotype fragment matrix from 10x BAM file

1. extract single chromosome (20) data

        python getMolecules.py -d 20000 NA12878.chr20.bam NA12878.chr20.molecules.bed
    

2. run extractHAIRS with --barcode 1 and --bed option

        ./FragmentCut --bam NA12878.chr20.bam --VCF NA12878_phased_variants.vcf.chr20.hets --bed NA12878.chr20.molecules.bed --barcode 1 --out chr20.frags --ref ~/Public/tools/reference-genomes/hg19.geuvadis.new.fa > chr20.log

3. run HapCUT2 using fragment file from step 2 
