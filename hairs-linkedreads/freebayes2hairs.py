from __future__ import print_function
import os,sys,math
import pysam
import itertools
import copy
import subprocess

DEBUG=False;
INDELS=1;
COMPEX_INDELS=1
MIN_BQ=10; MAX_BQ=40;
MIN_MQ=20;
SINGLE_READS=0;
HETS_ONLY=1; # output alleles for heterozygous variants only
LINKED_READS=1;
FREEBAYES='~/CODE/JOINTCODE-coral/HapCUT2/build/freebayes'

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
def read_vcf(vcf,QVTHRESH=0):
	varindex = {}; n=0;
	varlist = []
	File = open(vcf);
	for line in File: 
		if line[0] == '#': continue;
		v= line.strip().split('\t');
		n +=1;
		alleles = [v[3]] + v[4].split(','); 
		genotype = v[9].split(':'); 
		if genotype[0] == '0/1' or genotype[0] == '1/0' or genotype == '0|1' or genotype == '1|0': het =1; 
		else: het =0;
		#print (v[1],v[3],v[4],genotype)
		#varindex[(int(v[1]),v[3])] = [alleles,n,line.strip()];
		varlist.append([int(v[1]),v[3],alleles,het,line.strip()]); 
		if n%50000==0: print(line,alleles)
	File.close()
	print ('VARIANTS',n);
	return varlist;
	# return varindex;

## given a list of haplotype sequences, generate a new list with 'ref' replaced by 'alt' starting at position 'pos'
def add_allele(haplist,pos,ref,alt,b1,b2):
	newset = copy.deepcopy(haplist); #haplist[:]; ## copy of current hapset, apply variant allele to each element of this set
	l1 = len(ref); 
	flag=0;
	for i in xrange(len(newset)): 
		hap = newset[i];
		if pos + l1 > len(hap[0]): 
			print ('DEBUG',len(hap[0]),pos,l1,len(hap[0]),ref,alt)
			flag =1;
			continue; 
		for i in xrange(1,l1): hap[0][i+pos] = ''; # delete ref allele if longer than 1 
		if hap[0][pos] == '': hap[0][pos] = alt[1:]; ## if deleted, don't add the first base
		else: hap[0][pos] = alt; 
		hap[1][b1]= b2; # update bit vector 
	#	print(haplist,newset);
	if flag ==0: return newset
	else: return []
		

### genomic window from freebayes log file corresponding to 1 or more variants
class Window(object):
	def __init__(self,start,end,refallele): 
		self.start = start; 
		self.end = end;
		self.ref = refallele;
		self.alleles = {}; 
		self.readlist = []; self.fraglist = []
		self.reads=0;
		self.v1 = 0; self.v2 = 0; # index of variants overlapping this window
		self.vars=0;
		self.haplotypes = {}; ## 2^n list of combinations of alleles where n = self.vars
		self.ignore=0;

	def free(self): ## empty the data structures
		self.alleles.clear();
		#self.readlist.clear();

	def print(self,varlist):
		print('WINDOW',self.start,self.end,self.ref,self.alleles,self.v1,self.v2);
		if self.vars ==0: print ('no variants');
		for var in xrange(self.v1,self.v2): print('VARIANT',var,varlist[var][0],varlist[var][0]+len(varlist[var][1]),varlist[var][1],varlist[var][2][1:]);
		print (self.haplotypes);
		#print (self.multiset)
		if self.v2-self.v1 > 1: print ('multiple variants',self.vars,'indels',self.indels,'\n');
		else: print ('singleton');
		
		
	def find_overlapping_vars(self,varlist,len_varlist,s):
		self.v1 = s; 
		while self.v1 < len_varlist and varlist[self.v1][0] < self.start: self.v1 +=1;
		self.v2 = self.v1; 
		while self.v2 < len_varlist and varlist[self.v2][0] < self.end: self.v2 +=1; 
		self.vars = self.v2-self.v1;
	
	def create_alleles(self,varlist): ## alternate approach
		offset = self.start;
		hapset = [(list(self.ref),[0]*self.vars)]; 

		self.indels=0;
		for v in xrange(self.v1,self.v2):
			indel=0;
			to_add = [];
			a=1;
			for alt in varlist[v][2][1:]: 
				if len(varlist[v][1]) != len(alt): indel=1;
				#print ('input',hapset);
				newset=add_allele(hapset,varlist[v][0]-offset,varlist[v][1],alt,v-self.v1,a);
				if DEBUG: print('editing',varlist[v][1],'to',alt,varlist[v][0],v-self.v1,a,'OUTPUT',newset)
				a +=1;
				for elem in newset: to_add.append(elem);
			for elem in to_add: hapset.append(elem);
			#s = varlist[i][0]-offset; e = s + len(varlist[
			self.indels +=indel;
			
		for hap in hapset: 
			seq = ''.join(hap[0]);
			binary = ''.join(str(x) for x in hap[1]);
			try: self.haplotypes[seq].append(binary)
			except KeyError: self.haplotypes[seq] = [binary]
			#print(seq,hap[1]),

	def create_haplotypes(self,varlist): ## not used anymore
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

	def allelotype_reads(self,varlist):
		missing=0.0; nreads=0.0;
		self.fraglist = []
		for read in self.readlist: 
			nreads+=1;
			qual = int(-10*float(read[2])/math.log(10));
			hapseq = read[3];
			try: altype = self.haplotypes[hapseq]; c =1;
			except KeyError: altype = []; c= 0; missing +=1;
			info = read[4].split(':'); l = len(info); readid = ':'.join(info[1:l-10]);
			if c ==1: 
				self.fraglist.append([readid]);
				for v in range(self.v1,self.v2): 
					if qual >= MIN_BQ and (varlist[v][3] == 1 or HETS_ONLY ==0): self.fraglist[-1].append([v+1,altype[0][v-self.v1],qual]); 
			if DEBUG or c==0: print(readid,self.v1,self.v2,altype,'C',c,qual,hapseq);
		if missing > 0.2*nreads and missing >= 2: 
			self.ignore=1;
			print ('SUMMARY',missing,nreads);

	def update_fragments(self,fragments):
		if self.ignore ==1: return -1;
		#for frag in self.fraglist: print ('FRAGMENT',frag);
		for frag in self.fraglist: 
			for var in frag[1:]: 
				try: 
					fragments[frag[0]].add(var);
				except KeyError: 
					fragments[frag[0]] = Fragment(frag[0],0);
					fragments[frag[0]].add(var);
		return 1;


class Fragment:
	def __init__(self,fid,n):
		self.fid = fid;
		self.n = n; ## number of variants
		self.unique=0;
		self.varlist = []
		self.blocks=0;  # 3 
		self.qvstring = None; # FFFFFFF
		self.offsets = []; ## 564,570,572
		self.lengths =[]; # 4,2,1
		self.alleles = []; ## 0001,10,1
		self.barcode = None
	
	def add(self,vartuple):
		self.n +=1;
		self.varlist.append(vartuple);

	def prepare_for_hapcut(self,outfile=sys.stdout):
		self.varlist.sort();
		# identify duplicates and merge into unique list, then generate the offsets, alleles
		ulist =[]; c =0;
		for i in xrange(0,self.n):
			if c ==0 or self.varlist[i][0] != ulist[-1][0]: 
				ulist.append([self.varlist[i][0],self.varlist[i][1],self.varlist[i][2]]);
				c +=1;
			else: ## same variant read twice
				if self.varlist[i][1] == ulist[-1][1]: 
					newqv = self.varlist[i][2]+ulist[-1][2]; 
					ulist[-1][2] = min(newqv,MAX_BQ);
				else: 
					if (DEBUG): print('MISMATCH',file=outfile);
					ulist.pop(-1); c -=1; 
		self.unique = c;
		if self.unique > 0: 
			self.blocks =1; self.lengths.append(1);
			for i in range(self.unique-1): 
				if ulist[i][0]+1 != ulist[i+1][0]: self.blocks +=1; self.lengths.append(1);
				else: self.lengths[-1] +=1;
			c=0;
			for i in range(self.blocks):
				self.offsets.append(ulist[c][0]); 
				self.alleles.append(''.join([ulist[a][1] for a in range(c,c+self.lengths[i])]))
				c += self.lengths[i];
			self.qvstring = ''.join([chr(ulist[a][2]+33) for a in range(self.unique)]) 
			print(self.fid,self.varlist,'\n',self.fid,ulist,file=sys.stdout);
			#print('UNIQUE',self.blocks,self.lengths,self.offsets,self.alleles,ulist,self.qvstring,file=outfile);
	
	def print(self,outfile=sys.stdout):
		
		## output format is slightly modified 2 A00741:30:H7JKTDRXX:1:2130:22263:14340 2 TGTCATGTTGCTGCTTCG -1 1479 00 1482 0 FF5
		print (self.blocks,self.fid,end=' ',file=outfile);
		if LINKED_READS ==1: print ('2',self.barcode,-1,end=' ',file=outfile);
		for i in range(self.blocks): print (self.offsets[i],self.alleles[i],end=' ',file=outfile)
		print (self.qvstring,end='\n',file=outfile)



###################################################################################################################

def run_freebayes(bamfile='/media/drive2/HAPCUT-related/TELLseq-project/NovaSeqRun_April2019/output.bam.rg.bam',vcf='/media/drive2/HAPCUT-related/TELLseq-project/NovaSeqRun_April2019/longline.vcf.gz',fasta='~/Public/tools/reference-genomes/hg38.giab.fasta',region='chr19'):
	logfile = 'fb.' + region + '.log';
	fbfile = 'fb.' + region + '.vcf';
	command=FREEBAYES + ' -f ' + fasta + ' --haplotype-basis-alleles ' + vcf + ' -@ ' + vcf + ' ' + bamfile + ' -r ' + region + ' 2> ' + logfile + ' 1> ' + fbfile;
	print('running command',command,file=sys.stderr);
	subprocess.call(command,shell=True); 
	# ~/CODE/JOINTCODE-coral/HapCUT2/build/freebayes -f ~/Public/tools/reference-genomes/hg38.giab.fasta --haplotype-basis-alleles longline.vcf.gz -@ longline.vcf.gz output.bam.rg.bam -r chr20 2> chr20.fb.log 1> chr20.fb.vcf

def extract_fragments_logfile(logfile,varlist):
	fragments = {}; ## indexed by read-id 
	File = open(logfile);
	len_varlist = len(varlist);
	v=0;
	windows=0;
	lines=0;
	for line in File: 
		read = line.strip().split('\t');
		if read[0] == 'refr_seq': 
			# print previous window
			if windows > 0 and W.vars > 0: 
				W.allelotype_reads(varlist);
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


def print_fragments(fragments,outfilename='fragments.txt'):
	outfile = open(outfilename,'w');
	if LINKED_READS ==1: 
		pass;
		## add barcodes to each fragment in 'fragments'
	for frag in fragments.keys(): 	
		fragments[frag].prepare_for_hapcut(outfile);
		if fragments[frag].unique > 1 or (SINGLE_READS ==1 and fragments[frag].unique > 0): fragments[frag].print(outfile);
	outfile.close();


def add_barcodes(fragments,bamfile='/media/drive2/HAPCUT-related/TELLseq-project/NovaSeqRun_April2019/output.bam.rg.bam',region='chr19'):
	print('reading bam file to add barcodes for linked-reads',bamfile,file=sys.stderr)
	pybamfile=pysam.AlignmentFile(bamfile,'rb');
	n=0;
	for read in pybamfile.fetch(region):
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


if len(sys.argv) < 3: 
	print >>sys.stderr,'python reduced_matrix.py logfile vcf_file';
	sys.exit()
run_freebayes(); sys.exit();
varlist = read_vcf(sys.argv[2]);
fragments=extract_fragments_logfile(sys.argv[1],varlist);
add_barcodes(fragments)
print_fragments(fragments);



