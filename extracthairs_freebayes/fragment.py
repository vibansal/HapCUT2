from __future__ import print_function
import os,sys,math

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
		self.mq = 0;
		self.isize = 0;
		self.filter =False;
	
	def add(self,vartuple):
		self.n +=1;
		self.varlist.append(vartuple);

	def prepare_for_hapcut(self,outfile=sys.stdout,max_bq=40,DEBUG=0):
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
					ulist[-1][2] = min(newqv,max_bq);
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
	
	def print(self,outfile=sys.stdout,LINKED_READS=0):
		
		## output format is slightly modified 2 A00741:30:H7JKTDRXX:1:2130:22263:14340 2 TGTCATGTTGCTGCTTCG -1 1479 00 1482 0 FF5
		if not self.filter: 
			print (self.blocks,self.fid,end=' ',file=outfile);
			if LINKED_READS ==1: print ('2',self.barcode,-1,end=' ',file=outfile);
			for i in range(self.blocks): print (self.offsets[i],self.alleles[i],end=' ',file=outfile)
			print (self.qvstring,end='\n',file=outfile)



