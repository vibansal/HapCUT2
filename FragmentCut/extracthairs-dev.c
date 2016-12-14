
/* program to extract haplotype fragments from dilution pool sequencing, author Vikas Bansal */
/* main code is in file fosmidbam_hairs.c */
/* this code can easily be merged with extracthairs.c since it calls a different function and arguments are same */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
//#define _GNU_SOURCE

#include "hashtable.h"
#include "readfasta.h"
#include "bamread.h"
#include "sam.h"
#include "readvariant.h"
#include "hapfragments.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MISSING_QV =0;
int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS =  2000; // maximum insert size
int MIN_IS =  0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs 
int BSIZE = 500;
int IFLAG = 0;
int MAXFRAG = 500000;
int VARIANTS = 0;
int VCFformat = 0;
int PARSEINDELS =0;
int SINGLEREADS =0;
int FOSMIDS = 1;
//int QVoffset = 33; declared in samread.h
FILE* logfile;
int PFLAG = 1;
int PRINT_FRAGMENTS = 1;
FILE* fragment_file;
int TRI_ALLELIC = 0;
int BED_FILE = 0;
int READ_CHROM_IND = 1;  // read the sequence for each chromosome one by one or in a single pass 
char* regions; // global varaible pointing to --regions chr1:534333-54324500 option  
float SAMPLE_READS = 1.0;
int OVERLAP_FRAGS =0;

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "parsebamread.c"
#include "fosmidbam_hairs.c" // code for parsing fosmid pooled sequence data 

//int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist);

void print_options()
{
	fprintf(stderr,"\n PROGRAM TO extract haplotype informative reads (HAIRS) from coordinate sorted BAM files \n\n");
	fprintf(stderr,"./extract_hairs [options] --bam reads.sorted.bam --VCF variants.VCF   > output.fragments \n\n");
	fprintf(stderr,"=============== PROGRAM OPTIONS ======================================== \n\n");
	fprintf(stderr,"--qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33 \n");
	fprintf(stderr,"--mbq <INT> : minimum base quality to consider a base for haplotype fragment, default 13\n");
	fprintf(stderr,"--mmq <INT> : minimum read mapping quality to consider a read for phasing, default 20\n");
	fprintf(stderr,"--VCF <FILENAME> : variant file with genotypes for a single individual in VCF format\n");
	fprintf(stderr,"--variants : variant file in hapCUT format (use this option or --VCF option but not both), this option will be phased out in future releases\n");
	fprintf(stderr,"--maxIS <INT> : maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n");
	fprintf(stderr,"--minIS <INT> : minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n");
	fprintf(stderr,"--PEonly <0/1> : do not use single end reads, default is 0 (use all reads)\n");
	fprintf(stderr,"--indels <0/1> : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n");
	fprintf(stderr,"--noquality <INTEGER> : if the bam file does not have quality string, this value will be used as the uniform quality value, default 0 \n");
	//fprintf(stderr,"--triallelic <0/1> : print information about , default 0 \n");
	fprintf(stderr,"--ref <FILENAME> : reference sequence file (in fasta format), optional but required for indels, should be indexed using samtools\n");
	fprintf(stderr,"--mask <FILENAME> : reference sequence file (in fasta format) with mappability information (see http://lh3lh3.users.sourceforge.net/snpable.shtml), should be indexed using samtools\n");
	fprintf(stderr,"--out <FILENAME> : output filename for haplotype fragments, if not provided, fragments will be output to stdout\n\n");
	fprintf(stderr,"--bsize : BLOCK_SIZE for dynamic programming\n\n");
	fprintf(stderr,"--minsep : maximum distance between consecutive fragments to merge them into single cluster\n");
	fprintf(stderr,"--sample 0.5 : sub-sample reads from BAM-file\n");
}


int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist)
{
	int i = read->tid;
	if (reflist->ns > 0)
	{
		reflist->current = i;
		if (i >= reflist->ns || i < 0 || strcmp(reflist->names[i],read->chrom) !=0)
		{
			reflist->current = -1;
			for (i=0;i<reflist->ns;i++)
			{
				if (strcmp(reflist->names[i],read->chrom) ==0) { reflist->current = i; break; }
			}
		}
	}
	return 1;
}


int init_fastafile(char* fastafile,REFLIST* reflist)
{
        reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
        if (strcmp(fastafile,"None") == 0) return -1;
        if (read_fastaheader(fastafile,reflist) > 0) // all data structures are allocated
        {
		if (READ_CHROM_IND ==0) read_fasta(fastafile,reflist); // this will not work for genome mask file due to upper case 
		else { 
	                fprintf(stderr,"opening fasta file %s \n",fastafile);
        	        reflist->fp = fopen(fastafile,"r");
		}
        }
        return 1;
}

// new function added 08/26/16, print haplotype fragments using list of intervals from bedfile
int print_fragments_bedfile(struct alignedread** readlist,int s,int e,FRAGMENT* flist,VARIANT* varlist,REFLIST* reflist)
{
	find_matepair(readlist,s,e,MAX_IS); // link the paired-end reads together 
	struct FOSMID fosmid; 
	int fi=reflist->first_interval_chrom[readlist[s]->tid]; int chrom = reflist->intervallist[fi].chrom; 
	int firstread=s,lastread=s,start,end,intervals=0,firstread0=0; 
	while (fi < reflist->intervals && reflist->intervallist[fi].chrom == chrom)
	{
		start = reflist->intervallist[fi].start; end = reflist->intervallist[fi].end; 
		while (firstread < e && readlist[firstread]->position < start) firstread++; firstread0 = firstread; 
		lastread = firstread;  while (lastread < e && readlist[lastread]->position < end) lastread++; 

		if (BARCODE ==1) // additional check 
		{
			// find first read that has same barcode as the interval we are looking at, this read's barcode will be used as reference in generate_single_fragment function (no need to pass additional parameter)
			while (firstread < lastread)
			{
				//if (readlist[firstread]->barcode != NULL && fi ==1 )fprintf(stderr,"barcodes |%s| |%s| \n",readlist[firstread]->barcode,reflist->intervallist[fi].annotation);
				if (readlist[firstread]->barcode != NULL && strcmp(readlist[firstread]->barcode,reflist->intervallist[fi].annotation) == 0 ) break;
				else firstread++;
			}
		}
		if (lastread-firstread >= 1) 
		{
			if (reflist->intervallist[fi].annotation != NULL) fprintf(stdout,"%d-%d %d %d-%d %s #ofreads %d\n",start,end,firstread0,firstread,lastread,reflist->intervallist[fi].annotation,lastread-firstread);
			print_reads_window(readlist,firstread,lastread,flist,varlist,1); 
			fosmid.firstread = firstread; fosmid.lastread = lastread; fosmid.ws = end-start; fosmid.reads_window =lastread-firstread; 
			fosmid.start = start; fosmid.end = end; 
			generate_single_fragment(readlist,flist,varlist,&fosmid,"");
		}
		firstread = firstread0; // go back since intervals are sorted based on start position not barcodes
		fi++; intervals++;
	}
	fprintf(stderr,"done processing intervals for chrom %d, intervals %d \n",readlist[s]->tid,intervals);
	return 1;
}

void init_bamfile_index(char* bamfile,samfile_t *fp,char* regions,bam_index_t *idx,int reg[])  // for reading indexed bam files 
{
	int ref=-1,beg=0,end=0;
	bam_parse_region(fp->header,regions,&ref,&beg,&end);
	if (ref < 0)
	{
		fprintf(stderr,"invalid region for bam file %s \n",regions); exit(0);
	}
	reg[0]= ref; reg[1] = beg; reg[2] = end; 

}

// process all reads for a single chromosome since calculation of read-density requires data... extract haplotype informative reads from sorted bam file // 
int parse_bamfile_fosmid(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist,char* maskfile)
{
	// if bedfile is provided, we should simply print fragments 
	// print_reads_window(readlist,fosmidlist[i].firstread,fosmidlist[i].lastread,flist,varlist,1); generate_single_fragment(readlist,flist,varlist,&fosmidlist[i],"");

	REFLIST* genomemask = (REFLIST*)malloc(sizeof(REFLIST));  
	if (BED_FILE ==0) init_fastafile(maskfile,genomemask); 
	else  // maskfile actually points to bed file 
	{
		init_fastafile("None",genomemask); char* bedfile = maskfile; 
		if (read_bedfile(bedfile,reflist) == -1) 
		{
			BED_FILE = 0; fprintf(stderr," bed file not found.. exiting \n\n\n"); return -1;
		}
	}

	fprintf(stderr,"reading sorted bamfile %s for fosmid pool\n",bamfile);
	int reads=0;
	int MAX_READS = 2000000; // 10 million for now
	int READS_INCREMENT = 1000000;
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
	struct alignedread** readlist = calloc(MAX_READS,sizeof(struct alignedread*));
	for (reads=0;reads<MAX_READS;reads++) readlist[reads] = calloc(1,sizeof(struct alignedread)); 
	struct alignedread* read_pt;

	int max_fragments = MAX_READS/5; 
	FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*10000);
	FRAGMENT* flist = (FRAGMENT*)malloc(sizeof(FRAGMENT)*max_fragments); int fragments =0;

	int chrom=0; int r=0,i=0,f=0;
	int prevchrom=-1; int prevtid = -1; int prevposition = -1; // position of previous read in sorted bam file
	int lastread = 0; int terminate = 0; int lastvariantread=-1; 

	int filtered_reads[3] ={0,0,0};
	int flag1 = BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|2048; int bitflag=0; //fprintf(stderr,"flag1 %d %d\n",flag1,1028&flag1);

	samfile_t *fp; 
	bam_index_t *bam_idx; // bam file index
        bam_iter_t bam_iter; // for reading indexed bam files 
	bam1_t *b;
	int rf=0;
	if ((fp = samopen(bamfile, "rb", 0)) == 0) { fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; }

	if (regions == NULL) b = bam_init1(); 
	else  // read sub-set of indexed bam file
	{
		int reg[3] = {-1,0,0}; int ref=-1,beg=0,end=0;
		if ( (bam_idx  = bam_index_load(bamfile)) ==0) { fprintf(stderr,"unable to load bam index for file %s\n",bamfile); exit(0); }
		bam_parse_region(fp->header,regions,&ref,&beg,&end); reg[0]= ref; reg[1] = beg; reg[2] = end; 
	//	init_bamfile_index(bamfile,fp,regions,bam_idx,reg);  // for reading indexed bam files 
		b = bam_init1(); 
		bam_iter = bam_iter_query(bam_idx,reg[0],reg[1],reg[2]);
		fprintf(stderr,"region for bam file is %d:%d-%d \n",reg[0],reg[1],reg[2]);
	}

	while (1)
	{
		if (regions ==NULL) rf =samread(fp,b); else rf  = bam_iter_read(fp->x.bam,bam_iter,b);
		if (rf < 0) break;
		if (drand48() > SAMPLE_READS) continue; 

		fetch_func(b, fp,readlist[r]); 
		bitflag = readlist[r]->flag & flag1; //fprintf(stderr,"read flag %d %d \n",readlist[r]->flag,bitflag);
		if (bitflag)  // unmapped reads, PCR/optical dups are ignored
		{
			filtered_reads[0] +=1;
			free_readmemory(readlist[r]); continue;
		}
		else if (readlist[r]->mquality < MIN_MQ)  // ignore poorly mapped reads 
		{
			filtered_reads[1] +=1;
			free_readmemory(readlist[r]); continue;
		}
		else if (FOSMIDS ==2 && readlist[r]->clipped*2 >= readlist[r]->alignedbases) // excess soft clipped bases
		{
			//fprintf(stderr,"read %s %d\n",readlist[r]->readid,readlist[r]->clipped);
			filtered_reads[2] +=1;
			free_readmemory(readlist[r]); continue;
		}
		// find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
		// if too many reads, break off when distance between adjacent reads is large 
		if (readlist[r]->tid != prevtid) 
		{

			if (prevtid >=0) // process reads from previous chromosome 
			{
				fprintf(stderr,"reads in buffer %d process reads to identify fragments\n",r);
				if (r-lastread >= 1000) 
				{ 
					if (BED_FILE ==1) print_fragments_bedfile(readlist,lastread,r,flist,varlist,reflist);
					else process_reads_chrom(readlist,lastread,r,flist,varlist,reflist,genomemask); 
				}
				// free up reads in list and move up index 
				for (i=lastread;i<r;i++) free_readmemory(readlist[i]); 
				read_pt = readlist[0]; readlist[0] = readlist[r]; readlist[r] = read_pt; r = 0;
				for (i=0;i<fragments;i++) 
				{
					free(flist[i].alist); free(flist[i].id); 
				}
				fprintf(stderr,"free memory for reads from chrom %d cleaning up of fragment list %d\n",prevtid,fragments);
				fragments =0; 
				if (reflist->ns > 0 && prevtid < reflist->ns) free(reflist->sequences[prevtid]); 
				if (genomemask->ns > 0 && prevtid < reflist->ns) free(genomemask->sequences[prevtid]); 
				//return 1;  // comment this out to process all chromosomes
			}
			if (reflist->ns > 0 && readlist[r]->tid < reflist->ns && readlist[r]->tid>=0) 
			{
				if (READ_CHROM_IND ==1) read_chromosome(reflist,readlist[r]->tid,reflist->fp); // read the next chromosome
				reflist->current = readlist[r]->tid; // assign current variable to tid
			}
			// causing segfault when reading 'NT' chromosomes... 08/04/16
			if (genomemask->ns > 0 && readlist[r]->tid >=0)  // segfault here if the two genomes are not identical 
			{
				if (READ_CHROM_IND ==0 || read_chromosome_mask(genomemask,readlist[r]->tid,genomemask->fp) >= 0) genomemask->current = readlist[r]->tid; // assign current variable to tid
				else genomemask->current = -1; 
			}

			chrom = getindex(ht,readlist[r]->chrom); 
			get_chrom_name(readlist[r],ht,reflist); // reflist->current has chromosome index, -1 if no reflist 
			lastread = r;
		}
		else 	chrom = prevchrom;

		fragment.variants =0; //fragment.id = readlist[r]->readid;
		if (chrom >=0)  extract_variants_read(readlist[r],ht,chromvars,varlist,1,&fragment,chrom,reflist);
		if (fragment.variants > 0) 
		{ 
			add_fragment(flist,&fragment,readlist[r],fragments); readlist[r]->findex = fragments; fragments++; 
			if (fragments >= max_fragments-10) 
			{
				// realloc larger sized list 
				fprintf(stderr,"need to increase size of fragment list %d ",fragments);
				flist = (FRAGMENT*)realloc(flist,sizeof(FRAGMENT)*(max_fragments+READS_INCREMENT/5)); 
				max_fragments += READS_INCREMENT/5;
				fprintf(stderr,"new size is %d \n",max_fragments);
			}
		}
		else readlist[r]->findex = -1; // we could free read memory here if needed, don't need sequence, quality, etc non-informative reads 

		reads+=1; if (reads%2000000 ==0) fprintf(stderr,"processed %d reads, useful fragments | filtered reads %d,%d,%d \n",reads,filtered_reads[0],filtered_reads[1],filtered_reads[2]);
		prevchrom = chrom; prevtid = readlist[r]->tid; 
		if (readlist[r]->IS == 0) prevposition = readlist[r]->position; 
		else if (readlist[r]->IS > 0) prevposition = readlist[r]->position + readlist[r]->IS; // outer end of fragment r1....r2

		if (BED_FILE ==1) // discard reads that don't cover a variant
		{
			if (fragment.variants ==0) free_readmemory(readlist[r]); 
			else r++; 
		}
		else r++;
		if (r >= MAX_READS) // increase size of buffer 
		{
			fprintf(stderr,"limit on max reads %d exceeded.. need to increase size of buffer \n",MAX_READS);
			readlist = (struct alignedread**)realloc(readlist,(MAX_READS+READS_INCREMENT)*sizeof(struct alignedread*));
			for (reads=MAX_READS;reads<MAX_READS+READS_INCREMENT;reads++) readlist[reads] = calloc(1,sizeof(struct alignedread)); 
			MAX_READS += READS_INCREMENT;
		}
	}

	if (r-lastread >= 100) 
	{
		if (BED_FILE ==1) print_fragments_bedfile(readlist,lastread,r,flist,varlist,reflist); // bed file of fragment intervals is provided as input 
		else process_reads_chrom(readlist,lastread,r,flist,varlist,reflist,genomemask); 
	}
	for (i=lastread;i<r;i++) free_readmemory(readlist[i]); 
	for (reads=0;reads<MAX_READS;reads++) free(readlist[reads]); free(readlist);
	if (reflist->ns > 0 && prevtid >=0 && reflist->sequences[prevtid] != NULL) { free(reflist->sequences[prevtid]); fclose(reflist->fp); } 
	if (genomemask->ns > 0 && prevtid >=0 && genomemask->sequences[prevtid] != NULL) { free(genomemask->sequences[prevtid]); fclose(genomemask->fp); } 
	free(genomemask);free(reflist);
	bam_destroy1(b);
	return 1;
}


int main (int argc, char** argv)
{
	char bamfile[1024]; strcpy(bamfile,"/media/drive2/Haplotyping/NA12878-SOLID-fosmid/aligned-bams/pool_SPA1.novoalign.sorted.bam"); 
	char variantfile[1024]; strcpy(variantfile,"/media/drive2/Haplotyping/NA12878-SOLID-fosmid/NA12878.hg18.snps.vcf.hets");
	char fastafile[1024]; strcpy(fastafile,"/home/vbansal/Public/tools/reference-genomes/1000genomes-huref/human_b36_male.fa");
	char maskfile[1024]; strcpy(maskfile,"/media/drive2/Haplotyping/genome-mask/hg18-50bp-files/hg18.50bp.mask.fa"); 
	char bedfile[1024]; strcpy(maskfile,"None"); // for intervals 

	int readsorted = 0;
	char* sampleid = (char*)malloc(1024); sampleid[0] = '-'; sampleid[1] = '\0';
	int samplecol=10; // default if there is a single sample in the VCF file
	int i=0,variants=0,hetvariants=0;
	char** bamfilelist = NULL; int bamfiles =0;  regions = NULL;

	logfile = NULL; fragment_file = stdout; // write fragments to this file if it is present
	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)        bamfiles++; 
		else if (strcmp(argv[i],"--variants") ==0)        strcpy(variantfile,argv[i+1]);
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--mask") ==0 || strcmp(argv[i],"--mappability") ==0)        strcpy(maskfile,argv[i+1]);
		else if (strcmp(argv[i],"--VCF") ==0 || strcmp(argv[i],"--vcf") ==0)    {     strcpy(variantfile,argv[i+1]); VCFformat =1; }
		else if (strcmp(argv[i],"--mbq") ==0)       MINQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mmq") ==0)       MIN_MQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--maxIS") ==0)       MAX_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--minIS") ==0)       MIN_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--PEonly") ==0)       PEONLY = 1;  // discard single end mapped reads 
		else if (strcmp(argv[i],"--indels") ==0)       PARSEINDELS = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--pflag") ==0)      IFLAG  = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--qvoffset") ==0)       QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--out") == 0 || strcmp(argv[i],"-o") ==0) fragment_file = fopen(argv[i+1],"w");
		else if (strcmp(argv[i],"--logfile")==0 || strcmp(argv[i],"--log") ==0) logfile = fopen(argv[i+1],"w");  
		else if (strcmp(argv[i],"--singlereads")==0) SINGLEREADS = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--maxfragments")==0) MAXFRAG = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--noquality")==0) MISSING_QV = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--triallelic")==0) TRI_ALLELIC = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--fosmids") == 0 || strcmp(argv[i],"--fosmid") ==0) 
		{
			FOSMIDS = atoi(argv[i+1]); 
		}
		else if (strcmp(argv[i],"--firstpass") == 0) FIRSTPASS = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--comparephase") == 0 || strcmp(argv[i],"--compare") ==0) COMPARE_PHASE = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--bsize") == 0) BLOCK_SIZE = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--bs1") == 0) BS1 = atoi(argv[i+1]);  // number of sub-bins per block 
		else if (strcmp(argv[i],"--minreads") == 0) MIN_READS_PER_FRAGMENT = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--minsep") == 0) MIN_SEPARATION = atoi(argv[i+1]); 
		else if (strcmp(argv[i],"--fragpenalty") == 0) block_penalty_global = atof(argv[i+1]); 
		else if (strcmp(argv[i],"--readdensity") == 0) read_density_global = atof(argv[i+1]); 
		else if (strcmp(argv[i],"--barcode") == 0) { BARCODE = atoi(argv[i+1]); fprintf(stderr,"reads have barcodes \n"); } 
		else if (strcmp(argv[i],"--bed") == 0) { strcpy(bedfile,argv[i+1]); fprintf(stderr,"bed file provided \n"); BED_FILE = 1; } 
		else if (strcmp(argv[i],"--regions") == 0) { regions = argv[i+1]; fprintf(stderr,"regions arg provided \n"); } 
		else if (strcmp(argv[i],"--sample") == 0) { SAMPLE_READS = atof(argv[i+1]);  } 
	}
	if (bamfiles > 0 && strcmp(variantfile,"None") !=0)
	{
		bamfilelist = (char**)malloc(sizeof(char*)*bamfiles); 
		for (i=0;i<bamfiles;i++) bamfilelist[i] = (char*)malloc(1024);
		bamfiles=0;
		for (i=1;i<argc;i+=2)
		{
			if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)     strcpy(bamfilelist[bamfiles++],argv[i+1]);
		}
		fprintf(stderr,"\n extracting haplotype informative reads from bamfiles %s minQV %d minMQ %d maxIS %d \n\n",bamfilelist[0],MINQ,MIN_MQ,MAX_IS);
	}
	else
	{
		print_options(); return -1;
	}
	init_parameters(FOSMIDS);  // initialize parameters based on data type 09/04/16

	HASHTABLE ht; ht.htsize = 7919;  init_hashtable(&ht);
	VARIANT* varlist;
	int chromosomes=0;

	// if we can modify code to extract subset of VCF file corresponding to --regions, it would make processing individual chromosomes much simpler... 
	// easily done in python rather than VCF... 
	if (VCFformat ==1)
	{
		variants = count_variants(variantfile,sampleid,&samplecol); 
		if (variants < 0) return -1; 
		varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
		chromosomes = read_variantfile(variantfile,varlist,&ht,&hetvariants,samplecol); 
	}
	else
	{
		variants = count_variants_oldformat(variantfile);
		if (variants < 0) return -1; 
		varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
		chromosomes = read_variantfile_oldformat(variantfile,varlist,&ht,variants);
	}
	// variants is set to hetvariants only, but this is not correct since 
	VARIANTS = variants;  
	// there are two options, we include all variants in the chromvars datastructure but only use heterozygous variants for outputting HAIRS 
	// variant-id should correspond to line-number in VCF file since that will be used for printing out variants in Hapcut 

	//	fprintf(stderr,"read %d variants from file %s chromosomes %d\n",snps,argv[1],chromosomes);
	CHROMVARS* chromvars  = (CHROMVARS*)malloc(sizeof(CHROMVARS)*chromosomes);
	build_intervalmap(chromvars,chromosomes,varlist,VARIANTS);

	// read reference fasta file for INDELS, 
	REFLIST* reflist = (REFLIST*)malloc(sizeof(REFLIST)); init_fastafile(fastafile,reflist);
 
	if (readsorted ==0 && bamfiles > 0)
	{
		for (i=0;i<bamfiles;i++) 
		{
			//if (FOSMIDS ==0) parse_bamfile_sorted(bamfilelist[i],&ht,chromvars,varlist,reflist);
			if (BED_FILE ==0) parse_bamfile_fosmid(bamfilelist[i],&ht,chromvars,varlist,reflist,maskfile); // fosmid pool bam file 
			else parse_bamfile_fosmid(bamfilelist[i],&ht,chromvars,varlist,reflist,bedfile);// bedfile is provided
		}
	}
	if (logfile != NULL) fclose(logfile);
	if (fragment_file != NULL && fragment_file != stdout) fclose(fragment_file);
	return 0;
}


