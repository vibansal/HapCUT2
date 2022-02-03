
### first coded June 10 2019, Vikas Bansal

program that uses the freebayes log file (-d option) to allelotype reads for variants in an input VCF file

 for each window (identified using 'refr-seq'): 
   1. identify start and end and use vcf to generate all possible haplotype sequences for variants overlapping this window
   2. match each allele to a haplotype sequence and convert to individual variant alleles 
 need to account for overlapping reads, output PE reads together

- implemented for linked-read data but can also be useful in general

