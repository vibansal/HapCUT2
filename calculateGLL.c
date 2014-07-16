#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "hashtable.h"
#include "readfasta.h"
#include "bamread.h"
#include "sam.h"
#include "readvariant.h"
#include "hapfragments.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS =  1000; // maximum insert size
int MIN_IS =  0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs 
int BSIZE = 500;
int IFLAG = 0;
int MAXFRAG = 500000;
int VARIANTS = 0;
int VCFformat = 0;
int PARSEINDELS =0;
int SINGLEREADS =1;
int FOSMIDS = 0;
//int QVoffset = 33; declared in samread.h
FILE* logfile;
int PFLAG = 1;
int PRINT_FRAGMENTS = 0;
int POOL_SIZE = 2;// 
char* GROUPNAME;

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "parsebamread.c"

//disabled sam file reading
//#include "samhairs.c" // has two functions that handle sam file parsing 

double AGILENT_REF_BIAS = 0.5;
int update_GLLs(FRAGMENT* fragment,VARIANT* varlist)
{
	int i=0,v=0,j=0; double ep =0,e=0,f=0;
	for (i=0;i<fragment->variants;i++)
	{
		if ((int)fragment->alist[i].qv < (MINQ + QVoffset)) continue;  // filter on base quality already done in extract_variants_

		v = fragment->alist[i].varid; ep = pow(0.1,((double)fragment->alist[i].qv-QVoffset)/10);
		
		for (j=0;j<=POOL_SIZE;j++)
		{
			f = (double)j*(1.0-AGILENT_REF_BIAS); f /= f + (double)(POOL_SIZE-j)*AGILENT_REF_BIAS;
			//f = (double)j/POOL_SIZE;
                        e = (1.0-ep)*(1.0-f) + ep*f; 
			if (fragment->alist[i].allele =='0') varlist[v].GLL[j] += log10(e);
			else if (fragment->alist[i].allele =='1') varlist[v].GLL[j] += log10(1.0-e);
		//fprintf(stdout,"varlist %d %d %f %f %f %f %d\n",i,v,varlist[v].GLL[0],varlist[v].GLL[1],varlist[v].GLL[2],qv,fragment->alist[i].allele);
		}
	}
	return 1;
}

void print_options();
int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist);

void print_options()
{
	fprintf(stderr,"\n PROGRAM TO extract haplotype informative reads (HAIRS) from coordinate sorted BAM files \n\n");
	fprintf(stderr,"./extract_hairs [options] --bam reads.sorted.bam --VCF variants.VCF   > output.fragments \n\n");
	fprintf(stderr,"=============== PROGRAM OPTIONS ======================================== \n\n");
	fprintf(stderr,"--qvoffset : quality value offset, 33/64 depending on how quality values were encoded, default is 33 \n");
	fprintf(stderr,"--mbq  : minimum base quality to consider a base for haplotype fragment, default 13\n");
	fprintf(stderr,"--mmq : minimum read mapping quality to consider a read for phasing, default 20\n");
	fprintf(stderr,"--VCF : variant file with genotypes for a single individual in VCF format\n");
	fprintf(stderr,"--variants : variant file in hapCUT format (use this option or --VCF option but not both), this option will be phased out in future releases\n");
	fprintf(stderr,"--maxIS : maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n");
	fprintf(stderr,"--minIS : minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n");
	fprintf(stderr,"--PEonly 0/1 : do not use single end reads, default is 0 (use all reads)\n");
	fprintf(stderr,"--indels 0/1 : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n");
	fprintf(stderr,"--ref : reference sequence file (in fasta format), optional but required for indels, should be indexed using samtools\n");
	//fprintf(stderr,"--out : output file for haplotype informative fragments (hairs)\n\n");
}



// extract haplotype informative reads from sorted bam file //
// need to discard reads that are marked as duplicates using flag //
int parse_bamfile_sorted(char* bamfile,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,REFLIST* reflist)
{
	fprintf(stderr,"reading sorted bamfile %s \n",bamfile);
	int reads=0;
	struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
	
	int i=0; int sl=0; int chrom=0;
	int v1,v2; int absIS;
	int prevchrom=-1; int prevtid = -1;

	FRAGMENT* flist = (FRAGMENT*)malloc(sizeof(FRAGMENT)*MAXFRAG); int fragments =0; int prevfragments =0;
	FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*1000);

	samfile_t *fp;
	if ((fp = samopen(bamfile, "rb", 0)) == 0) { fprintf(stderr, "Fail to open BAM file %s\n", bamfile); return -1; }
	bam1_t *b = bam_init1();

	while (samread(fp, b) >= 0)
	{
		fetch_func(b, fp,read);
		if ((read->flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) || read->mquality < MIN_MQ) 
		{
			free_readmemory(read); continue;
		}
		// find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
		if (read->tid != prevtid)
		{
			chrom = getindex(ht,read->chrom); // doing this for every read, can replace this by string comparison ..april 4 2012
			i = read->tid;
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
		}
		else chrom = prevchrom;

		absIS = (read->IS < 0) ? -1*read->IS: read->IS; 
		// add check to see if the mate and its read are on same chromosome, bug for contigs, july 16 2012
		if (((read->flag & 4) ==0)) // single read
		{
			fragment.variants =0; v1 =0; v2=0; 
			if (chrom >=0) 
			{
				fragment.id = read->readid;
				v1 = extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
				if (fragment.variants > 0) 
				{
					if (read->IS > 0 && (( read->flag & 8) == 0)) 
					{
						for (i=0;i<fragment.variants;i++) 
						{
							if (varlist[fragment.alist[i].varid].position > read->mateposition) fragment.alist[i].qv = 33;
							else break;
						}
					}
					//if (read->IS > 0 && varlist[fragment.alist[0].varid].position > read->mateposition) fprintf(stdout,"OPE read %s %d %d \n",read->readid,varlist[fragment.alist[0].varid].position,read->mateposition);
					update_GLLs(&fragment,varlist); 
					//fprintf(stdout," useful frag \n");
				}
			}
		}

		reads+=1; if (reads%2000000 ==0) fprintf(stderr,"processed %d reads, useful fragments %d\n",reads,fragments);
		prevchrom = chrom; prevtid = read->tid;
		free_readmemory(read);
	}
	if (fragments > 0) 
	{
		fprintf(stderr,"final cleanup of fragment list: %d current chrom %s %d \n",fragments,read->chrom,read->position);
		clean_fragmentlist(flist,&fragments,varlist,-1,read->position,prevchrom);
	}
	bam_destroy1(b);
}


int main (int argc, char** argv)
{
	char samfile[1024]; char bamfile[1024]; char variantfile[1024]; char fastafile[1024];
	strcpy(samfile,"None"); strcpy(bamfile,"None"); strcpy(variantfile,"None"); strcpy(fastafile,"None");
	GROUPNAME = NULL;
	int readsorted = 0;
	char* sampleid = (char*)malloc(1024); sampleid[0] = '-'; sampleid[1] = '\0';
	int samplecol=10; // default if there is a single sample in the VCF file
	int i=0,j=0,variants=0,hetvariants=0;
	char** bamfilelist = NULL; int bamfiles =0; 

	logfile = NULL;
	for (i=1;i<argc;i+=2)
	{
		if (strcmp(argv[i],"--bam") ==0 || strcmp(argv[i],"--bamfile") ==0)        bamfiles++; 
		else if (strcmp(argv[i],"--variants") ==0)        strcpy(variantfile,argv[i+1]);
		else if (strcmp(argv[i],"--reffile") ==0 || strcmp(argv[i],"--ref") ==0)        strcpy(fastafile,argv[i+1]);
		else if (strcmp(argv[i],"--VCF") ==0 || strcmp(argv[i],"--vcf") ==0)    {     strcpy(variantfile,argv[i+1]); VCFformat =1; }
		else if (strcmp(argv[i],"--sorted") ==0)       readsorted = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mbq") ==0)       MINQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--mmq") ==0 || strcmp(argv[i],"--minmq") ==0 )       MIN_MQ = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--maxIS") ==0)       MAX_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--minIS") ==0)       MIN_IS = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--PEonly") ==0)       PEONLY = 1;  // discard single end mapped reads 
		else if (strcmp(argv[i],"--indels") ==0)       PARSEINDELS = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--pflag") ==0)      IFLAG  = atoi(argv[i+1]);  // allow indels in hairs
		else if (strcmp(argv[i],"--qvoffset") ==0)       QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--logfile")==0 || strcmp(argv[i],"--out") ==0) logfile = fopen(argv[i+1],"w");  
		else if (strcmp(argv[i],"--singlereads")==0) SINGLEREADS = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--maxfragments")==0) MAXFRAG = atoi(argv[i+1]);  
		else if (strcmp(argv[i],"--poolsize")==0 || strcmp(argv[i],"-p")==0) 
		{
			POOL_SIZE = atoi(argv[i+1]);  
			if (POOL_SIZE > 2) AGILENT_REF_BIAS = 0.50;
		}
		else if (strcmp(argv[i],"--groupname") == 0) 
		{
			GROUPNAME = (char*)malloc(1024); strcpy(GROUPNAME,argv[i+1]); 
		}
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
		fprintf(stderr,"\n extracting reads covering SNPs from bamfile %s minQV %d minMQ %d maxIS %d \n\n",bamfilelist[0],MINQ,MIN_MQ,MAX_IS);
	}
	else
	{
		print_options(); return -1;
	}

	HASHTABLE ht; ht.htsize = 7919;  init_hashtable(&ht);
	VARIANT* varlist;
	int chromosomes=0;


	variants = count_variants_oldformat(variantfile);
	if (variants < 0) return -1; else varlist = (VARIANT*)malloc(sizeof(VARIANT)*variants);
	chromosomes = read_variantfile_oldformat(variantfile,varlist,&ht,variants);
	for (i=0;i<variants;i++)
	{
		varlist[i].GLL = calloc(POOL_SIZE,sizeof(double));
		for (j=0;j<POOL_SIZE;j++) varlist[i].GLL[j] = 0.0;
	}

	// variants is set to hetvariants only, but this is not correct since 
	VARIANTS = variants;  
	// there are two options, we include all variants in the chromvars datastructure but only use heterozygous variants for outputting HAIRS 
	// variant-id should correspond to line-number in VCF file since that will be used for printing out variants in Hapcut 

	//	fprintf(stderr,"read %d variants from file %s chromosomes %d\n",snps,argv[1],chromosomes);
	CHROMVARS* chromvars  = (CHROMVARS*)malloc(sizeof(CHROMVARS)*chromosomes);
	build_intervalmap(chromvars,chromosomes,varlist,VARIANTS);

	// read reference fasta file for INDELS//
	REFLIST* reflist = (REFLIST*)malloc(sizeof(REFLIST)); 
	reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
	if (strcmp(fastafile,"None") != 0)
	{
		if (read_fastaheader(fastafile,reflist) > 0) 
		{
			reflist->sequences = (char**)malloc(sizeof(char*)*reflist->ns);
			for (i=0;i<reflist->ns;i++)
			{
				reflist->sequences[i] = (char*)malloc(reflist->lengths[i]+1);
				if (i < 5) fprintf(stderr,"contig %s length %d\n",reflist->names[i],reflist->lengths[i]);
			}
			read_fasta(fastafile,reflist);
		}
	}
	//return 1;
	if (readsorted ==0 && bamfiles > 0)
	{
		for (i=0;i<bamfiles;i++) 
		{
			parse_bamfile_sorted(bamfilelist[i],&ht,chromvars,varlist,reflist);
		}
	}
	if (logfile != NULL) fclose(logfile);

	int xor = pow(2,16)-1; 
	double sum =0; int counts[4];
	for (i=0;i<variants;i++)
	{
		//if (varlist[i].type ==0) continue;
		//if (varlist[i].genotype[0] == varlist[i].genotype[2]) continue;
		counts[0] = varlist[i].A1>>16; counts[1] = varlist[i].A1 & xor; counts[2] = varlist[i].A2>>16; counts[3] = varlist[i].A2 & xor;
		if (counts[0] + counts[1]+counts[2]+counts[3] < 1) continue; 
		sum = varlist[i].GLL[0]; 
		for (j=1;j<=POOL_SIZE;j++) 
		{ 
			if (varlist[i].GLL[j] < sum) sum += log10(1.0+pow(10,varlist[i].GLL[j]-sum)); 
			else sum = varlist[i].GLL[j] + log10(1.0+pow(10,sum-varlist[i].GLL[j]));
		}
		for (j=0;j<=POOL_SIZE;j++) fprintf(stdout,"%0.4f ",varlist[i].GLL[j]-sum);
		fprintf(stdout,"%s %d %s %s %d:%d %d:%d\n",varlist[i].chrom,varlist[i].position,varlist[i].allele1,varlist[i].allele2,counts[0],counts[1],counts[2],counts[3]);
		//fprintf(stdout,"variant %d %s %d %d %s %s %d:%d %d:%d\n",i+1,varlist[i].chrom,varlist[i].position-1,varlist[i].type,varlist[i].allele1,varlist[i].allele2,varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor);
	}
	return 0;
}


