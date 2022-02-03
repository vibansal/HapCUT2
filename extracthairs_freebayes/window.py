from __future__ import print_function
import os,sys,math
import itertools
import copy
from fragment import Fragment

## given a list of haplotype sequences, generate a new list with 'ref' replaced by 'alt' starting at position 'pos'
def add_allele(haplist,pos,ref,alt,b1,b2,DEBUG=0):
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
	
	def create_alleles(self,varlist,DEBUG=0): ## alternate approach
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

	def allelotype_reads(self,varlist,min_bq=10,hets_only=1,indels=1,DEBUG=0):
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
					if indels ==0 and (len(varlist[v][1]) != len(varlist[v][2][1])): continue; # ignore indels 
					if qual >= min_bq and (varlist[v][3] == 1 or hets_only ==0): self.fraglist[-1].append([v+1,altype[0][v-self.v1],qual]); 
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




