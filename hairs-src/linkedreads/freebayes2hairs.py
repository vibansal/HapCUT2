from __future__ import print_function
import os,sys,math
import pysam
import copy
import subprocess
import argparse
from fragment import Fragment
from window import Window

COMPEX_INDELS=1
SINGLE_READS=1;
HETS_ONLY=1; # output alleles for heterozygous variants only
LINKED_READS=1;
FREEBAYES='~/CODE/JOINTCODE-coral/HapCUT2/build/freebayes'
runfb=0; # run freebayes

"""
first coded June 10 2019, Vikas Bansal

program that uses the freebayes log file (-d option) to allelotype reads for variants in an input VCF file

 for each window (identified using 'refr-seq'): 
   1. identify start and end and use vcf to generate all possible haplotype sequences for variants overlapping this window
   2. match each allele to a haplotype sequence and convert to individual variant alleles 
 need to account for overlapping reads, output PE reads together

v = [[0,1,2],[1],[0,1]]; 
for i in itertools.product(*v): print i  

"""
def read_vcf(vcf,region):
	varindex = {}; n=0;
	varlist = []
	vcx = pysam.VariantFile(vcf);
        for line in vcx.fetch(region=region):
		#v= line.strip().split('\t');
		#alleles = [v[3]] + v[4].split(',');
		for sample in line.samples.items(): 
			genotype=sample[1]['GT'];
			#print(line.alleles,sample[1]['GT'],file=sys.stdout);
		if genotype[0] == genotype[1]: het = 0;
		elif genotype[0] != genotype[1] and genotype[0] < 2 and genotype[1] < 2: het = 1; 
		else: het =0;
		#genotype = v[9].split(':'); 
		#if genotype[0] == '0/1' or genotype[0] == '1/0' or genotype == '0|1' or genotype == '1|0': het =1; 
		#else: het =0;
		varlist.append([line.pos,line.ref,line.alleles,het]); 
		#varlist.append([int(v[1]),v[3],alleles,het,line.strip()]); 
		n +=1;
		if n%50000==0: print(line,line.alleles)
	print ('VARIANTS',n,file=sys.stderr);
        vcx.close();
	return varlist;


###################################################################################################################

def run_freebayes(bamfile,vcf,region,fasta):
	logfile = 'fb.' + region + '.log';
	fbfile = 'fb.' + region + '.vcf';
	command=FREEBAYES + ' -f ' + fasta + ' --haplotype-basis-alleles ' + vcf + ' -@ ' + vcf + ' ' + bamfile + ' -r ' + region + ' 2> ' + logfile + ' 1> ' + fbfile;
	print('running command',command,file=sys.stderr);
	subprocess.call(command,shell=True); 
	# ~/CODE/JOINTCODE-coral/HapCUT2/build/freebayes -f ~/Public/tools/reference-genomes/hg38.giab.fasta --haplotype-basis-alleles longline.vcf.gz -@ longline.vcf.gz output.bam.rg.bam -r chr20 2> chr20.fb.log 1> chr20.fb.vcf
	return logfile;

def extract_fragments_logfile(varlist,options):
	fragments = {}; ## indexed by read-id 
	File = open(options.logfile);
	len_varlist = len(varlist);
	v=0;
	windows=0;
	lines=0;
	for line in File: 
		read = line.strip().split('\t');
		if read[0] == 'refr_seq': 
			# print previous window
			if windows > 0 and W.vars > 0: 
				W.allelotype_reads(varlist,min_bq=options.mbq,hets_only=options.hets_only,indels=options.indels);
				W.print(varlist);
				W.update_fragments(fragments);

			start = int(read[1])+1; refallele = read[3]; end = start + len(refallele);
			W = Window(start,end,refallele); 
			W.find_overlapping_vars(varlist,len_varlist,v); 
			W.create_alleles(varlist);
			#W.create_haplotypes(varlist);
			v = W.v1; 
			windows +=1;	

		elif read[0] == 'haplo_obs': 
			allele = read[3]; 
			try: W.alleles[allele]+=1; 
			except KeyError: W.alleles[allele] = 1; 
			W.reads +=1; 
			W.readlist.append(read);

		lines+=1; 
		if lines%100000 ==0: 
			#print ('exiting early for test purposes...remove this for production code',file=sys.stderr);
			#break; 
			pass;
	File.close();
	print('processed freebayes logfile to get fragments',lines,file=sys.stderr);
	return fragments;


def print_fragments(fragments,outfilename='fragments.txt',singlereads=1,linked_reads=0):
	outfile = open(outfilename,'w');
	for frag in fragments.keys(): 	
		fragments[frag].prepare_for_hapcut(outfile);
		if fragments[frag].unique > 1 or (singlereads==1 and fragments[frag].unique > 0): fragments[frag].print(outfile,linked_reads);
		print(fragments[frag],frag);
	outfile.close();


def add_barcodes(bamfile,region,fragments):
	print('reading bam file to add barcodes for linked-reads',bamfile,file=sys.stderr)
	pybamfile=pysam.AlignmentFile(bamfile,'rb');
	n=0;
	for read in pybamfile.fetch(region=region):
		try:
			barcode =read.get_tag('BX');
			try: 
				fragments[read.query_name].barcode = barcode;
				#print('found match',read.query_name,barcode,file=sys.stderr);
			except KeyError: pass; 
			#fragments.update({read.query_name: barcode}); 
			#print (read.query_name,barcode,file=sys.stderr);
		except KeyError:pass;
		n +=1;
		if n%1000000==0: print ('processed',n,'reads',file=sys.stderr);
	pybamfile.close();

	
#######################################################################################################################	


def add_program_options(parser): 
	required = parser.add_argument_group('required arguments')
	parser.add_argument('-b', '--bam', action='store', required=True,default="",dest='bam',help='Bam file')
	parser.add_argument('--region', action='store', dest='region',required=True,default='', help='Region in format chr16:1270000-1300000');
	parser.add_argument('--VCF', action='store',required=True, default="",dest='VCFfile',help='VCF file for genotyping')
	parser.add_argument('-f','--fasta',action='store',required=True, type=str,default="", dest='fasta',help='Reference fasta file')
	parser.add_argument('-o','--out', action='store', default="temp",type=str,dest='outfile',help='output file with fragments')
	parser.add_argument('--runfb', action='store', default=0,type=int,help='run freebayes')
	parser.add_argument('--logfile', action='store', default=None,type=str,help='freebayes log file')
	parser.add_argument('--indels', action='store', default=0,type=int,help='include indels [0/1]')
	parser.add_argument('--singlereads', action='store', default=1,type=int,help='output reads covering only a single variant')
	parser.add_argument('--hets_only', action='store', default=1,type=int,help='output alleles for only heterozygous variants, default=1')
	parser.add_argument('--linked', action='store', default=0,type=int,help='linked-reads [0/1]')
	parser.add_argument('--mbq', action='store', default=10,type=int,help='minimum base quality for alleles')
	parser.add_argument('--maxbq', action='store', default=40,type=int,help='maximum base quality for alleles')
	parser.add_argument('--debug', action='store', default=0,type=int,help='debug, verbose output')

#####################################################################################################

parser = argparse.ArgumentParser(); add_program_options(parser); options = parser.parse_args(); 

varlist = read_vcf(options.VCFfile,options.region);

if options.logfile == None and options.runfb ==0: 
		print ("logfile needs to be provided or freebayes needs to be run",file=sys.stderr);
		sys.exit()

if options.runfb==1: options.logfile = run_freebayes(options.bam,options.VCFfile,options.region,options.fasta); 

fragments=extract_fragments_logfile(varlist,options);

if options.linked ==1: 
	add_barcodes(options.bam,options.region,fragments);
	#python3 ~/CODE/JOINTCODE-coral/HapCUT2/utilities/LinkFragments.py --fragments fragments.txt --VCF phasing-40kb/chr20.vcf --bam output.bam.rg.bam --distance  40000 --outfile chr20.linked.fragments.new

print_fragments(fragments,outfilename=options.outfile,singlereads=options.singlereads,linked_reads=options.linked);



