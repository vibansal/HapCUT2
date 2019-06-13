from __future__ import print_function
import os,sys,math
#import pysam
import itertools

"""
program that uses the freebayes log file (-d option) to allelotype reads for variants in an input VCF file

 for each window (identified using 'refr-seq'): 
   1. identify start and end and use vcf to generate all possible haplotype sequences for variants overlapping this window
   2. match each allele to a haplotype sequence and convert to individual variant alleles 

 need to account for overlapping reads, output PE reads together

 first coded June 10 2019
 duplicates some of the functions in extracthairs but avoids realignment, anchor finding, etc

v = [[0,1,2],[1],[0,1]]; 
for i in itertools.product(*v): print i  

"""
def read_vcf(vcf,QVTHRESH=0):
	varindex = {}; n=0;
	varlist = []
	File = open(vcf);
	for line in File: 
		if line[0] == '#': continue;
		v= line.strip().split('\t');
		n +=1;
		alleles = [v[3]] + v[4].split(','); 
		varindex[(int(v[1]),v[3])] = [alleles,n,line.strip()];
		varlist.append([int(v[1]),v[3],alleles,line.strip()]); 
		if n%50000==0: print(line,alleles)
	File.close()
	print ('VARIANTS',n);
	return varlist;
	# return varindex;



class Window(object):
	def __init__(self,start,end,refallele): 
		self.start = start; 
		self.end = end;
		self.ref = refallele;
		self.alleles = {}; 
		self.readlist = []; 
		self.reads=0;
		self.v1 = 0; self.v2 = 0; # index of variants overlapping this window
		self.vars=0;
		self.haplotypes = {}; ## 2^n list of combinations of alleles where n = self.vars
		self.multiset = []

	def free(self): ## empty the data structures
		self.alleles.clear();
		#self.readlist.clear();

	def print(self,varlist):
		print(self.start,self.end,self.ref,self.alleles,self.v1,self.v2);
		if self.vars ==0: print ('no variants');
		for var in xrange(self.v1,self.v2): print('VARIANT',var,varlist[var][0],varlist[var][0]+len(varlist[var][1]),varlist[var][1],varlist[var][2][1:]);
		print (self.haplotypes);
		print (self.multiset)
		if self.v2-self.v1 > 1: print ('multiple variants\n');
		else: print ();
		
		
	def find_overlapping_vars(self,varlist,len_varlist,s):
		self.v1 = s; 
		while self.v1 < len_varlist and varlist[self.v1][0] < self.start: self.v1 +=1;
		self.v2 = self.v1; 
		while self.v2 < len_varlist and varlist[self.v2][0] < self.end: self.v2 +=1; 
		self.vars = self.v2-self.v1;
	
	def create_alleles(self,varlist): ## alternate approach
		offset = self.start;
		hapset = [(list(self.ref),[0]*self.vars)]; 
		print(hapset);
		for i in xrange(self.v1,self.v2):
			newset = hapset[:]; ## copy of current hapset, apply variant allele to each element of this set
			#s = varlist[i][0]-offset; e = s + len(varlist[
		## 

	def create_haplotypes(self,varlist): ## what to do when a single variant has three alleles... 
		offset = varlist[self.v1][0]; 
		binaryset = []
		e = offset; 
		for v in xrange(self.v1,self.v2): 
			coding = []
			if varlist[v][0] > e: 
				#print('extra',self.ref[e-offset:varlist[v][0]-offset]);
				self.multiset.append([self.ref[e-offset:varlist[v][0]-offset]]);
			if varlist[v][0] < e: print ('OVERLAPPING variants',self.ref);

			s=varlist[v][0]; e = varlist[v][0] + len(varlist[v][1]);
			for i in xrange(len(varlist[v][2])): coding.append(str(i));
			self.multiset.append(varlist[v][2]); 
			binaryset.append(coding)
			#print('VAR',s,e,varlist[v][2]);
		if self.end > e: 
			seq = self.ref[e-offset:self.end-offset];
			self.multiset.append([seq]);
		binary = []; ## binary strings corresponding to each haplotype (00,11) useful for outputting the fragments
		for i in itertools.product(*binaryset): binary.append(''.join(i));
		#for i in itertools.product(*multiset): self.haplotypes.append(''.join(i));
		n=0;
		for i in itertools.product(*self.multiset): 
			hap = ''.join(i);
			self.haplotypes[hap] = binary[n]; 
			n +=1;

	def allelotype_reads(self):
		for read in self.readlist: 
			qual = int(-10*float(read[2])/math.log(10));
			hapseq = read[3];
			try: altype = self.haplotypes[hapseq]; c =1;
			except KeyError: altype = ''; c= 0;
			info = read[4].split(':'); l = len(info); readid = ':'.join(info[1:l-10]);
			print(readid,self.v1,self.v2,altype,'C',c,qual,hapseq);


###################################################################################################################
		
def parse_logfile(logfile,varlist):
	File = open(logfile);
	len_varlist = len(varlist);
	v=0;
	windows=0;
	for line in File: 
		read = line.strip().split('\t');

		if read[0] == 'refr_seq': 

			# print previous window
			if windows > 0: 
				W.allelotype_reads();
				W.print(varlist);

			start = int(read[1])+1; refallele = read[3]; end = start + len(refallele);
			W = Window(start,end,refallele); 
			W.find_overlapping_vars(varlist,len_varlist,v); 
			W.create_alleles(varlist);
			W.create_haplotypes(varlist);
			v = W.v1; 
			windows +=1;	

		elif read[0] == 'haplo_obs': 
			allele = read[3]; 
			try: W.alleles[allele]+=1; 
			except KeyError: W.alleles[allele] = 1; 
			W.reads +=1; 
			W.readlist.append(read);
		
	File.close();
		

if len(sys.argv) < 3: 
	print >>sys.stderr,'python reduced_matrix.py logfile vcf_file';
	sys.exit()
varlist = read_vcf(sys.argv[2]);
parse_logfile(sys.argv[1],varlist);



