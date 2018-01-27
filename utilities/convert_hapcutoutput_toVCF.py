#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math

# last modified dec 6 2017 

# author Vikas Bansal, script to convert output file from HapCUT2 to VCF (phased format using 0|1:block_no)  
# python hapcut_to_vcf.py phased-haplotypes/chr1.filtered chr1.vcf > chr1.vcf.phased

# need to phase homozygotes, can be done after running hapcut 

if len(sys.argv) < 3: print 'python convert_output_vcf.py hapcut.out variants.vcf > variants.vcf.phased'; sys.exit();

# https://samtools.github.io/hts-specs/VCFv4.3.pdf
## PS = A phase set is defined as a set of phased genotypes to which this genotype belongs. Phased genotypes for an individual that are on the same chromosome and have the same PS value are in the same phased set. PS = position of first variant in each haplotype block 
## PQ (integer) = phasing quality, output by hapCUT2 
## PD (integer) = phase-depth, # of haplotype informative reads for variant, not yet added 

################## read hapcut output file for phased genome ##################

hapcut = open(sys.argv[1],'r'); phased_table = {}; blockfirst = -1;
for line in hapcut:
	if '****' in line: # end of previous block, print statistics 
		blockfirst = -1; 
	elif 'BLOCK' in line: # start of new block  
		blockinfo = line.strip();
		var = line.split(); blockid = var[2]; length = var[4]; 
	else:
		var = line.split(); #info = var[0].split('_'); 
		if var[1] != '-' and var[2] != '-': 
			if blockfirst ==-1: blockfirst = int(var[4]);
			phased_table[(var[3],int(var[4]))] = [var[1] + '|' + var[2],var[10],blockfirst]; 
			## (phased genotype,phasing quality,phased-block-id,original-genotype)
hapcut.close();



########################## read VCF file ###############################

VCF_file = open(sys.argv[2],'r');
for line in VCF_file: 
	if line[0]== '#': 
		print line, 
		continue;

	var = line.strip().split('\t'); 
	chrom = var[0]; position = var[1]; rsid = var[2]; ref = var[3]; alleles = var[4].split(','); allele2 = alleles[0]; 

	not_phased = 1;
	try: 
		phasedvar = phased_table[(chrom,int(position))];
		not_phased= 0;
	except KeyError: 	
		#print >>sys.stderr, "unphased variant",chrom,position;
		pass; 

	first_8_cols = "\t".join([var[i] for i in xrange(8)]);  
	C9 = var[8].split(':'); C10 = var[9].split(':');
	newC9 = ':'.join(C9[1:]); newC10 = ':'.join(C10[1:]);

	if not_phased ==1: 
		print first_8_cols + '\tGT:PQ:PS:' + newC9 +  '\t' + C10[0] + ':' +  '.' + ':' + '.' + ':' + newC10;
		continue; 
	else:
		## print phased output
		print first_8_cols + '\tGT:PQ:PS:' + newC9 +  '\t' + phasedvar[0] + ':' + phasedvar[1]  + ':' + `phasedvar[2]` + ':' + newC10;

VCF_file.close();

