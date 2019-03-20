#!/bin/bash

HAPCUT2=../../build/HAPCUT2
HAP_STATISTICS=../../utilities/calculate_haplotype_statistics.py


mkdir -p data/temp
# download the thousand genomes VCFs
echo "DOWNLOADING 1000 GENOMES VCFS FOR NA12878 TRIO"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/CEU.trio.2010_03.genotypes.vcf.gz \
      -O data/temp/NA12878_trio_genotypes_original.vcf.gz
# unzip the thousand genomes VCFs
gunzip -c data/temp/NA12878_trio_genotypes_original.vcf.gz > data/temp/NA12878_trio_genotypes_original.vcf

echo "FILTERING 1000 GENOMES VCFS FOR JUST NA12878 HETEROZYGOUS SITES"
# process the thousand genomes VCFs to be just NA12878 heterozygous sites
mkdir data/NA12878_hg18_VCFs
python3 create_NA12878_hg18_vcfs.py

echo "CHECKING THAT NA12878 VARIANT INDICES/COORDINATES MATCH DUITAMA PHASED DATA"
# download the Duitama et al 'phased matrix'. It has genomic coordinates for the variants
# we basically want to make sure that our processed NA12878 VCF
# has matching variant indices and genomic coordinates
# to Duitama's own Refhap-based phased data.
wget http://www.molgen.mpg.de/~genetic-variation/SIH/Data/haplotypes.tar.gz \
     -O data/temp/duitama_haplotypes.tar.gz
mkdir data/temp/duitama_haplotypes
tar -xzf data/temp/duitama_haplotypes.tar.gz -C data/temp/duitama_haplotypes
python3 sanity_check_variant_indices.py

echo "DOWNLOADING AND PROCESSING DUITAMA ET AL FRAGMENT FILES"
# download and unzip the Duitama et al phasing matrices (fragment files)
wget http://www.molgen.mpg.de/~genetic-variation/SIH/Data/phasing_matrices.tar.gz \
     -O data/temp/phasing_matrices.tar.gz
mkdir data/NA12878_fosmid_data_original
mkdir data/NA12878_fosmid_data_formatted
tar -xzf data/temp/phasing_matrices.tar.gz -C data/NA12878_fosmid_data_original

# remove first line of fragment files
# Duitama et al phasing matrices have the first line as matrix dimensions
# this does not work with HapCUT2
for i in {1..22} X
    do
    tail -n +2 data/NA12878_fosmid_data_original/chr${i}.matrix.SORTED \
    > data/NA12878_fosmid_data_formatted/chr${i}.matrix.SORTED
done

echo "RUNNING HAPCUT2 ON EACH CHROMOSOME OF DUITAMA ET AL DATA"
mkdir data/hapcut2_haplotypes
# run HapCUT2 on each chromosome
for i in {1..22} X; do
    $HAPCUT2 --fragments data/NA12878_fosmid_data_formatted/chr${i}.matrix.SORTED \
        --vcf data/NA12878_hg18_VCFs/chr${i}.vcf \
        --out data/hapcut2_haplotypes/chr${i}.hap
done

echo "COMPARING ASSEMBLED HAPLOTYPE TO 1000G PHASE DATA FOR ACCURACY"
python3 $HAP_STATISTICS -v1 data/NA12878_hg18_VCFs/chr{1..22}.vcf \
                            data/NA12878_hg18_VCFs/chrX.vcf \
                        -h1 data/hapcut2_haplotypes/chr{1..22}.hap \
                            data/hapcut2_haplotypes/chrX.hap \
                        -v2 data/NA12878_hg18_VCFs/chr{1..22}.vcf \
                            data/NA12878_hg18_VCFs/chrX.vcf
