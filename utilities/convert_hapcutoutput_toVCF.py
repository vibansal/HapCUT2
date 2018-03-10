#! /usr/bin/env python
import sys, os, glob, string, subprocess,time,math,argparse
from hap_fragments import read_hairs_file

# last modified 03/09/2018
# author Vikas Bansal

## take the HapCUT2 (or hapCUT) phased output file and the input VCF file 
## output = all variants in the VCF file, phased and unphased (phased format using 0|1:block_no)  

UNPHASE=1; ## don't keep the phasing information from original VCF 
PRESERVE_INFO=0;
HAIR_FILE = 0;
ADD_BARCODES=0; ## output barcodes associated for each variant (only for 10X or LFR data)
KEEP_IDS = ['GQ','AD','DP'];  # GT info fields to keep from original VCF
## PS = A phase set is defined as a set of phased genotypes to which this genotype belongs PS = position of first variant in each haplotype block 

################## read hapcut output file for phased genome ##################

def get_phased_blocks(hapcutfile):
	print >>sys.stderr, "reading hapcut blocks file",hapcutfile;
	phased_table = {};  ## index based on (chrom,pos,ref,alt)
	blockfirst = -1; blocks=0; 
	variant_index = {}; ## index based on integer ID in hairs file
	hapcut = open(hapcutfile,'r'); 
	for line in hapcut:
		if '****' in line: # end of previous block, print statistics 
			blockfirst = -1; 
		elif 'BLOCK' in line: # start of new block  
			blockinfo = line.strip();
			var = line.split(); blockid = var[2]; length = var[4];  blocks +=1;
		else:
			var = line.split(); #info = var[0].split('_'); 
			if var[1] != '-' and var[2] != '-':  ## PRUNED variant
				if blockfirst ==-1: blockfirst = int(var[4]);
				phased_table[(var[3],int(var[4]),var[5],var[6])] = [var[1] + '|' + var[2],var[10],blockfirst,int(var[0])]; 
				variant_index[int(var[0])] = [var[4],var[1],var[2],var[5],var[6],var[7],var[8],blocks-1,int(blockid),[]];
				## (phased genotype,phasing quality,phased-block-id,original-genotype)
	print >>sys.stderr,"\nfinished reading hapcut output file with",blocks,"blocks\n";
	hapcut.close();
	return [phased_table,variant_index]; 

########################## read VCF file ###############################

def read_VCF(vfile):
	VCF_file = open(vfile,'r');
	for line in VCF_file:
		if line[0]== '#': continue;
		var = line.split();
		chrom = var[0]; position = var[1]; rsid = var[2]; allele1 = var[3]; alleles = var[4].split(','); allele2 = alleles[0];
		genotypes = var[9].split(':');

		if genotypes[0] == '0|1': snptable[(var[0],int(var[1]),'c')] = [0,1,genotypes,1,0]; phased +=1;
		elif genotypes[0] == '1|0': snptable[(var[0],int(var[1]),'c')] = [1,0,genotypes,1,0]; phased +=1;
		elif genotypes[0] == '0/1' or genotypes[0] == '1/0': snptable[(var[0],int(var[1]),'c')] = [-1,-1,genotypes,0,0]; unphased +=1; # in this case, genotype is unphased 
		else: pass;
	VCF_file.close();
	print >>sys.stderr, 'comparing phase from VCF file to HAPCUT phase','unphased hets',unphased,'phased hets',phased;
	return snptable;

def print_vcf_headers(outfile):
	#print >>outfile,"##fileformat=VCFv4.0";
	print >>outfile,"##source=HapCUT2 phased haplotype blocks";
	#print >>outfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	print >>outfile,"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ID of Phase Set for Variant\">";
	print >>outfile,"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability at this variant is incorrectly phased relative to the haplotype\">";
	#print >>outfile,"##FORMAT=<ID=JQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability of a phasing switch error in gap prior to this variant\">";
	print >>outfile,"##FORMAT=<ID=BX,Number=.,Type=String,Description=\"Barcodes for this variant\">";
	print >>outfile,"##FORMAT=<ID=PD,Number=1,Type=Integer,Description=\" phased Read Depth\">";
	#print >>outfile,"##commandline=\"....\"";



def output_phased_VCF(vcffile,phased_table,varindex,outfile):
	OUTFILE = open(outfile,'w');
	VCF_file = open(vcffile,'r'); ## hapcut2 VCF
	unphased=0; phased=0;
	print >>sys.stderr, "reading VCF file, genotypes will be unphased:",vcffile;
	for line in VCF_file: 
		if line[0]== '#': 
			v = line.strip().split(); 
			if v[0] == '#CHROM':
				## print information about PS and PQ tags 
				print_vcf_headers(OUTFILE);
				
			print >>OUTFILE, line, 
			continue;

		var = line.strip().split('\t'); 
		chrom = var[0]; position = var[1]; rsid = var[2]; ref = var[3]; alleles = var[4].split(','); allele2 = alleles[0]; 
		genotype = var[9].split(':'); 
		het=0;
		if genotype[0][0] == '0' and genotype[0][2] == '1': het = 1;
		if genotype[0][0] == '1' and genotype[0][2] == '0': het = 1;

		not_phased = 1;
		try: 
			phasedvar = phased_table[(chrom,int(position),var[3],var[4])]; vid = phasedvar[-1]; 
			links= len(varindex[vid][-1]); fragments = varindex[vid][-1]; 
			#print vid,position,links,BXlist;
			ugenotype = phasedvar[0]; 
			not_phased= 0;
			if het ==1: phased +=1;
		except KeyError: 	
			if het ==1: unphased +=1;
			ugenotype = genotype[0][0] + '/' + genotype[0][2];
			#print >>sys.stderr, "unphased variant",chrom,position;
			pass; 

		C9 = var[8].split(':'); C10 = var[9].split(':'); 
		newC9 = ['GT']; newC10 = [ugenotype];
		for i in xrange(1,len(C9)):
			# GT:AD:DP:GQ:PL preserve these values from original VCF 
			if C9[i] in KEEP_IDS: newC9.append(C9[i]); newC10.append(C10[i]); 
			elif C9[i] != 'PS' and C9[i] != 'PQ' and PRESERVE_INFO ==1 : newC9.append(C9[i]); newC10.append(C10[i]); 
		if not_phased ==0: 
			newC9.append('PS'); newC9.append('PQ'); 
			newC10.append(phasedvar[1]); newC10.append(`phasedvar[2]`); 
			if HAIR_FILE ==1: newC9.append('PD');newC10.append(`links`);
			if ADD_BARCODES ==1 and HAIR_FILE ==1: 
				BXlist = ';'.join([fragments[i][1].split(':')[2] for i in xrange(links)]);
				newC9.append('BX'); newC10.append(BXlist);
		
		print >>OUTFILE, "\t".join([var[i] for i in xrange(8)]) + '\t' + ':'.join(newC9) + '\t' + ':'.join(newC10);

	VCF_file.close(); OUTFILE.close();
	print >>sys.stderr, "finished writing new phased VCF",outfile,"phased heterozygous variants",phased,"unphased",unphased;


#############################################################################################################

def parseargs():
    parser = argparse.ArgumentParser(description='## Convert hapcut2 phased block to a VCF file, adds Barcodes associated with each variant if fragment file is provided (TBI)\n')
    parser.add_argument('-f', '--fragments', nargs='?', default="None",type = str, help='file with hapcut2 fragments (generated using extractHAIRS)')
    parser.add_argument('-v', '--VCF', nargs='?', type = str, help='VCF file used as input to hapcut2')
    parser.add_argument('-i', '--hapcut', nargs='?', type = str, help='hapcut2 output file with phased haplotype blocks')
    parser.add_argument('-o', '--pVCF', nargs='?', type = str, help='output VCF file with phased blocks')
    parser.add_argument('-p', '--keep', nargs='?', default = "0",type = str, help='keep additional genotype information in phased VCF, default 0');
    #parser.add_argument('-b', '--addbarcodes', nargs='?', default = 0,type = int, help='print list of barcodes (10X or LFR) for each variant, default 0');

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# main function
# parse input and run function to call alleles
if __name__ == '__main__':
	args = parseargs()
	if int(args.keep) ==1: PRESERVE_INFO = 1;
	#if int(args.addbarcodes) ==1: ADD_BARCODES= 1; 
	[phased_table,varindex] = get_phased_blocks(args.hapcut);
	if args.fragments != "None": hairlist = read_hairs_file(args.fragments,varindex); HAIR_FILE = 1; 

	output_phased_VCF(args.VCF,phased_table,varindex,args.pVCF);

