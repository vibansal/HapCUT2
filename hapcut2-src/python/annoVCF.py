import os,sys,math,argparse
import pysam
from variant import Variant

def add_program_options(parser):
        required = parser.add_argument_group('required arguments')
        parser.add_argument('-v', action='store', default="",dest='vcf',help='vcf file filtered')
        parser.add_argument('-c', action='store', dest='chrom',default='chr20', help='name of chromosome');
        parser.add_argument('--bed', action='store', dest='bed',default='/media/NAS/IndividualGenomes/NA12878/GIABcalls/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz', help='GIAB bed file of valid regions');
        parser.add_argument('--giab', action='store', dest='giab_vcf',default='/media/NAS/IndividualGenomes/NA12878/GIABcalls/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz', help='GIAB VCF gold standard');


## annotate variants using VCF as true SNV or near indel (within 5-10 bp) 
def anno_vars(varlist,VARIANTS,CHROM,GIABVCF,BEDFILE):
        vcx = pysam.VariantFile(GIABVCF);
        varindex = {}; i=0;
        for row in vcx.fetch(CHROM):
                for sample in row.samples.items(): varindex[row.pos] = [row.alleles,sample[1]['GT'],'True'];
		i +=1;
        vcx.close();
	print >>sys.stderr,"finished building true variant index, variants=",i;

	print >>sys.stderr,"chrom=",CHROM,"GIAB",GIABVCF;

	## intersect bedfile with VCF 
	tbx = pysam.TabixFile(BEDFILE);
	intervals= []; I=0;
	for row in tbx.fetch(CHROM,parser=pysam.asBed()): intervals.append([row.start,row.end]); I+=1;
	tbx.close();
	i=0;
	for var in varlist:
		if VARIANTS[var].pos < intervals[i][0]: 
			#print VARIANTS[var].pos,'GIAB=0';
			VARIANTS[var].GIAB = 0; continue
		if VARIANTS[var].pos > intervals[i][0] and VARIANTS[var].pos < intervals[i][1]: 
			#print VARIANTS[var].pos,'GIAB=1';
			VARIANTS[var].GIAB = 1; continue; 
		while i< I and VARIANTS[var].pos > intervals[i][1]: i +=1; 
		if i==I: break;
		#print VARIANTS[var].pos,intervals[i],VARIANTS[var].GIAB,i

	TP=0;FP=0;missing=0
	for n in VARIANTS.iterkeys():
		pos = VARIANTS[n].pos;
		if VARIANTS[n].GIAB ==1: 
			try:
				VARIANTS[n].true = varindex[pos];
				TP +=1;
			except KeyError: 
				VARIANTS[n].true = []; 
				FP +=1; 
		else: 
			VARIANTS[n].true = [];
			missing +=1;
	print >>sys.stderr, "TP",TP,"FP",FP,'missing',missing;


def read_vcf(vcf,QVTHRESH=0):
	varlist = []; VARIANTS= {}; n=0; f=0;
	header = []
	File = open(vcf);
	for line in File: 
		v= line.strip().split('\t');
		if line[0] == '#': 
			header.append(line);
			continue;
		genotype = v[9].split(':');PS = v[9].split(':')[2]; QUAL = float(v[5]); ref = v[3]; alt = v[4]; pos = int(v[1]); 
		n +=1;
		varlist.append(n);
		VARIANTS[n] = Variant(v[0],pos,ref,alt,QUAL,genotype[0],v);
		VARIANTS[n].vid = n; 
	File.close()
	print >>sys.stderr, 'VARIANTS',f,'total',n;
	return [varlist,VARIANTS,header]; 

def out_vcf(varlist,VARIANTS,header):
	## output annotated VCF
	for line in header: print line,
	for var in varlist: VARIANTS[var].printvarfull();

###########################################################################################################
	
parser = argparse.ArgumentParser(); add_program_options(parser); r = parser.parse_args();
if r.vcf == "": 
	parser.print_help(sys.stderr)
	sys.exit();	
[varlist,VARIANTS,header] = read_vcf(r.vcf);
anno_vars(varlist,VARIANTS,r.chrom,r.giab_vcf,r.bed);
out_vcf(varlist,VARIANTS,header);
