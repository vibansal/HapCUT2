#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math,argparse

# last modified 02/06/2018
# author Vikas Bansal, script to convert output file from HapCUT2 to VCF (phased format using 0|1:block_no)  

"""
TODO 
1. if VCF file is already phased, need to remove phased genotypes
2. option to exclude variants below certain phasing quality
"""

# https://samtools.github.io/hts-specs/VCFv4.3.pdf
## PS = A phase set is defined as a set of phased genotypes to which this genotype belongs. Phased genotypes for an individual that are on the same chromosome and have the same PS value are in the same phased set. PS = position of first variant in each haplotype block 
## PQ (integer) = phasing quality, output by hapCUT2 
## PD (integer) = phase-depth, # of haplotype informative reads for variant, not yet added 

################## read hapcut output file for phased genome ##################

def get_phased_blocks(hapcutfile):
	print >>sys.stderr, "reading hapcut blocks file",hapcutfile;
	hapcut = open(hapcutfile,'r'); phased_table = {}; blockfirst = -1;
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
	return phased_table;

########################## read VCF file ###############################

def read_VCF(vcffile,phased_table,outfile):
	OUTFILE = open(outfile,'w');
	VCF_file = open(vcffile,'r');
	for line in VCF_file: 
		if line[0]== '#': 
			print >>OUTFILE, line, 
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
			print >>OUTFILE, first_8_cols + '\tGT:PQ:PS:' + newC9 +  '\t' + C10[0] + ':' +  '.' + ':' + '.' + ':' + newC10;
			continue; 
		else:
			## print phased output
			print >>OUTFILE, first_8_cols + '\tGT:PQ:PS:' + newC9 +  '\t' + phasedvar[0] + ':' + phasedvar[1]  + ':' + `phasedvar[2]` + ':' + newC10;
	VCF_file.close();
	OUTFILE.close();

def parseargs():
    parser = argparse.ArgumentParser(description='')
    #parser.add_argument('-f', '--fragments', nargs='?', type = str, help='file with unlinked hapcut2 fragments (generate using --10X 1 option in extractHAIRS)')
    parser.add_argument('-v', '--vcf', nargs='?', type = str, help='vcf file used as input to hapcut2')
    parser.add_argument('-i', '--hapcut', nargs='?', type = str, help='hapcut2 output file with phased haplotype blocks')
    parser.add_argument('-o', '--outvcf', nargs='?', type = str, help='output VCF file with phased blocks')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# main function
# parse input and run function to call alleles
if __name__ == '__main__':
	args = parseargs()
	phased_table = get_phased_blocks(args.hapcut);
	read_VCF(args.vcf,phased_table,args.outvcf);

