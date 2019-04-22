
import os,sys,pysam
from scipy  import stats

class Graph: 
	def __init__(self,nvar):
		self.n = nvar; 
		EDGES = {};    ### (v1,v2) = Edge object
		VARIANTS = {};  ### VARIANTS[12] = Variant object 
		varlist = []; # varlist[0] = 10; varlist[1] = 12; integer 

class Edge: ## between two variants 
	def __init__(self,v1,v2):
		self.v1 = v1; self.v2 = v2; 
		self.n = 0;
		self.counts = [0,0]; # 00+11,01+10
		self.fragments = []; ## ('0',13,'1',20) 
		self.G00=0.0; self.G01 = 0.0; ## likelihoods of two possible pairings 00/11  and 01/10 
		self.unphased = 0.0; self.delta =0.0;

## individual variant 
class Variant: 
	def __init__(self,contig,pos,ref,alt,qual,genotype,array):
		self.contig = contig; self.pos = pos; self.ref = ref; self.alt = alt; self.qual = qual; self.genotype= genotype;
		self.alleles = [];
		self.counts = [[0,0],[0,0]];
		self.SB = 1; # strand bias
		self.GLL = [0.0,0.0,0.0]; # genotype likelihoods 0/0 0/1 1/1 
		self.pGLL = [0.0,0.0,0.0,0.0,0.0]; ## 0|0 0|1 1|0 1|1 0/1
		self.LLreads = 0.0; ## likelihood of all reads covering this variant assuming unphased genotype
		self.GIAB = 0; self.true=0; ## valid region, variant present 
		self.fragments = []; ### list of fragments covering the variant, indexes into hairlist 
		self.sequence = '';
		self.vid = -1; 
		self.motif = ''
		self.het = 0;
		self.goodedges = 0; self.badedges=0; self.score = 0.0; ## for pairwise
		self.edges = []; ## list of variants with which there is an edge to this variant 
		self.phased=0; ## variant is phased or not 
		self.neighbors = [];
		self.array = array;  ## split of line 

	def calcSB(self):
	        [OR,pv] = stats.fisher_exact(self.counts); self.SB = pv; 
 
	def strandbias(self):
		self.counts = [[0,0],[0,0]];
		for allele in self.alleles:
			if allele[0] == '0' and allele[2] == '+': self.counts[0][0] +=1;
			elif allele[0] == '0' and allele[2] == '-': self.counts[0][1] +=1;
			elif allele[0] == '1' and allele[2] == '+': self.counts[1][0] +=1;
			elif allele[0] == '1' and allele[2] == '-': self.counts[1][1] +=1;
	        [OR,pv] = stats.fisher_exact(self.counts); self.SB = pv; 
		
	def getseq(self,FASTA):
		s = FASTA.fetch(self.contig,self.pos-1 -10,self.pos-1+11)
		self.sequence = s[0:10] + ':' + s[10] + ':' + s[11:]
		self.motif = s[6:16];

	def printvar(self):
		if len(self.ref) == len(self.alt): vtype = 'SNV';
		else: vtype = 'INDEL';
		print self.vid,self.pos,self.ref,self.alt,vtype,self.qual,self.genotype,self.GLL,self.counts;
	
	def printvarfull(self):
		if len(self.true) > 0: self.array[2] = 'GIAB=' + str(self.GIAB) + ':TP=1'
		else : self.array[2] = 'GIAB=' + str(self.GIAB) + ':TP=0'
		print '\t'.join(self.array);



