//  program to extract haplotype informative reads from sorted BAM file, input requirements: bamfile and variantfile
//#include "extracthairs.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<limits.h>
//#define _GNU_SOURCE

#include "hashtable.h"
#include "readfasta.h"
#include "bamread.h"
#include "htslib/sam.h"
#include "readvariant.h"
#include "hapfragments.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MISSING_QV = 0;
int MINQ = 13; // minimum base quality
int MIN_MQ = 20; // minimum read mapping quality
int MAX_IS = 1000; // maximum insert size
int MIN_IS = 0; // maximum insert size
int PEONLY = 0; // if this is set to 1, reads for which only one end is mapped are not considered for hairs
int BSIZE = 500;
int IFLAG = 0;
int MAXFRAG = 500000;
int VARIANTS = 0;
int VCFformat = 0;
int PARSEINDELS = 0;
int SINGLEREADS = 0;
int LONG_READS = 0;
int REALIGN_VARIANTS = 0;
int ESTIMATE_PARAMS = 0; // estimate realignment parameters from BAM 
//int QVoffset = 33; declared in samread.h
FILE* logfile;
int PFLAG = 1;
int PRINT_FRAGMENTS = 1;
FILE* fragment_file;
int TRI_ALLELIC = 0;
int VERBOSE = 0;
int PACBIO = 0; 
int USE_SUPP_ALIGNMENTS =0; // use supplementary alignments, flag = 2048
int SUM_ALL_ALIGN =0; // if set to 1, use sum of all alignments scoring forr local realignment 
int HOMOZYGOUS = 0; // also output alleles for homozygous variants, by default such variants are ignored

int* fcigarlist; // global variable

float MATCH = 1;
float MISMATCH = -1;
float INSERTION_OPEN = -1;
float INSERTION_EXTEND = -1;
float DELETION_OPEN = -1;
float DELETION_EXTEND = -1;


// DATA TYPE
// 0 : generic reads
// 1 : HiC
// 2: 10X ??
int DATA_TYPE = 0;

// NEW FORMAT FOR HIC
// same as old format but with more columns.
// column 3 is the data type (0 for normal, 1 for HiC)
// column 4 is the first SNP index of mate 2 (-1 if no mate 2)
// column 5 is the absolute insert size
int NEW_FORMAT = 0;

int PRINT_COMPACT = 1; // 1= print each fragment block by block, 0 = print variant by variant

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "estimate_hmm_params.c"
#include "realign_pairHMM.c" // added 11/29/2018
#include "parsebamread.c"
#include "realignbamread.c"
//#include "fosmidbam_hairs.c" // code for parsing fosmid pooled sequence data

Align_Params* AP; // global alignmnet params

//disabled sam file reading
//#include "samhairs.c" // has two functions that handle sam file parsing

void print_options();
int parse_bamfile_sorted(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist,char* regions);

void print_options() {
	fprintf(stderr, "\nExtract haplotype informative reads (HAIRS) from coordinate sorted BAM files \n\n");
	fprintf(stderr, "./extractHAIRS [options] --bam reads.sorted.bam --VCF variants.VCF --out fragment_file \n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "--qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33 \n");
	fprintf(stderr, "--mbq <INT> : minimum base quality to consider a base for haplotype fragment, default 13\n");
	fprintf(stderr, "--mmq <INT> : minimum read mapping quality to consider a read for phasing, default 20\n");
	fprintf(stderr, "--realign_variants <0/1> : Perform sensitive realignment and scoring of variants.\n");
	fprintf(stderr, "--hic <0/1> : sets default maxIS to 40MB, prints matrix in new HiC format\n");
	fprintf(stderr, "--10X <0/1> : 10X reads. NOTE: Output fragments MUST be processed with LinkReads.py script after extractHAIRS to work with HapCUT2.\n");
	fprintf(stderr, "--pacbio <0/1> : Pacific Biosciences reads. Similar to --realign_variants, but with alignment parameters tuned for PacBio reads.\n");
	fprintf(stderr, "--ONT, --ont <0/1> : Oxford nanopore technology reads. Similar to --realign_variants, but with alignment parameters tuned for Oxford Nanopore Reads.\n");
	fprintf(stderr, "--new_format, --nf <0/1> : prints matrix in new format. Requires --new_format option when running HapCUT2.\n");
	fprintf(stderr, "--VCF <FILENAME> : variant file with genotypes for a single individual in VCF format (unzipped) \n");
	//fprintf(stderr, "--variants : variant file in hapCUT format (use this option or --VCF option but not both), this option will be phased out in future releases\n");
	fprintf(stderr, "--maxIS <INT> : maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n");
	fprintf(stderr, "--minIS <INT> : minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n");
	fprintf(stderr, "--PEonly <0/1> : do not use single end reads, default is 0 (use all reads)\n");
	fprintf(stderr, "--indels <0/1> : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n");
	fprintf(stderr, "--noquality <INTEGER> : if the bam file does not have quality string, this value will be used as the uniform quality value, default 0 \n");
	fprintf(stderr,"--triallelic <0/1> : include variants with genotype 1/2 for parsing, default 0 \n");
	fprintf(stderr, "--ref <FILENAME> : reference sequence file (in fasta format, gzipped is okay), optional but required for indels, should be indexed\n");
	fprintf(stderr, "--out <FILENAME> : output filename for haplotype fragments, if not provided, fragments will be output to stdout\n");
	fprintf(stderr, "--region <chr:start-end> : chromosome and region in BAM file, useful to process individual chromosomes or genomic regions \n");
	fprintf(stderr, "--ep <0/1> : set to 1 to estimate HMM parameters from aligned reads (only with long reads), default = 0\n");
	fprintf(stderr, "--hom <0/1> : set to 1 to include homozygous variants for processing, default = 0 (only heterozygous) \n\n");
	//fprintf(stderr, "--sumall <0/1> : set to 1 to use sum of all local alignments approach (only with long reads), default = 1 \n\n");
}

void check_input_0_or_1(char* x){
	if (!(strcmp(x, "0") == 0 || strcmp(x, "1") == 0)){
		fprintf(stderr, "\nERROR: Invalid input \"%s\" for <0/1> option flag.\n",x);
		exit(1);
	}
}

// extract haplotype informative reads from sorted bam file //
// need to discard reads that are marked as duplicates using flag //

int parse_bamfile_sorted(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist,char* regions) {
	fprintf(stderr, "reading sorted bamfile %s \n", bamfile);
	int reads = 0;
	struct alignedread* read = (struct alignedread*) malloc(sizeof (struct alignedread));

	int rvalue =0;
	// estimate alignment parameters from BAM file ONT/pacbio reads only 12/3/18
	if (REALIGN_VARIANTS && ESTIMATE_PARAMS) rvalue = realignment_params(bamfile,reflist,regions,AP); 
	if (rvalue < -1) return -1;
	if (REALIGN_VARIANTS) fcigarlist = calloc(sizeof(int),400000);
	int i = 0;
	int chrom = 0; //int sl=0;
	// int v1,v2;
	int prevchrom = -1;
	int prevtid = -1;

	FRAGMENT* flist = (FRAGMENT*) malloc(sizeof (FRAGMENT) * MAXFRAG);
	int fragments = 0,VOfragments[2]={0,0}; // fragments overlapping variants
	int prevfragments = 0;
	FRAGMENT fragment;
	fragment.variants = 0;
	fragment.alist = (allele*) malloc(sizeof (allele)*16184);

    samFile *fp;
    bam_hdr_t *hdr;
    bam1_t *b;
    int ret; // for reading indexed bam files 
    hts_idx_t *idx; // bam file index
    hts_itr_t *iter = NULL;
    
    if ((fp = sam_open(bamfile, "r")) == NULL) {
        fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
        return -1;
    }
    if ((hdr = sam_hdr_read(fp)) == NULL) {
        fprintf(stderr, "Fail to read header from file %s\n", bamfile);
        sam_close(fp);
        return -1;
    }
    if (regions != NULL &&  (idx = bam_index_load(bamfile)) ==0) { fprintf(stderr,"unable to load bam index for file %s\n",bamfile); return -1; }
    if (regions == NULL) b = bam_init1();  // sam_read1(fp,hdr,b)
    else 
    {
        b = bam_init1();
        iter = sam_itr_querys(idx,hdr,regions);
        if (iter == NULL) { fprintf(stderr,"invalid region for bam file %s \n",regions); return -1; }
    }

    while (1) {
	if (!iter) ret = sam_read1(fp,hdr,b);  // read full bam file
	else ret = sam_itr_next(fp,iter,b);  // specific region
	//fprintf(stderr,"here %d %d\n",ret,iter);
	if (ret < 0) break;  
        fetch_func(b, hdr, read);
        if ((read->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) || read->mquality < MIN_MQ || (USE_SUPP_ALIGNMENTS == 0 && (read->flag & 2048))) {
            free_readmemory(read);
            continue;
        }
        // find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
        if (read->tid != prevtid) {
	    chrom = getindex(ht,read->chrom);  // this will return -1 if the contig name is not  in the VCF file
	    if (chrom < 0) fprintf(stderr,"chrom \"%s\" not in VCF file, skipping all reads for this chrom.... \n",read->chrom);
	    else 
	    {
		fprintf(stderr,"processing reads mapped to chrom \"%s\" \n",read->chrom);
	    }
		// doing this for every read, can replace this by string comparison ..april 4 2012
            i = read->tid;
            if (reflist->ns > 0) {
                reflist->current = i;
                if (i >= reflist->ns || i < 0 || strcmp(reflist->names[i], read->chrom) != 0) {
                    reflist->current = -1;
                    for (i = 0; i < reflist->ns; i++) {
                        if (strcmp(reflist->names[i], read->chrom) == 0) {
                            reflist->current = i;
                            break;
                        }
                    }
                }
            }
        } else chrom = prevchrom;
        //if (chrom_missing_index ==1) { prevtid = read->tid; free_readmemory(read); continue; }

        //////////
        	
		fragment.absIS = (read->IS < 0) ? -1 * read->IS : read->IS;
		// add check to see if the mate and its read are on same chromosome, bug for contigs, july 16 2012
		if ((read->flag & 8) || fragment.absIS > MAX_IS || fragment.absIS < MIN_IS || read->IS == 0 || !(read->flag & 1) || read->tid != read->mtid) // single read
		{
			fragment.variants = 0; // v1 =0; v2=0;
			if ( (read->flag & 16) ==16) fragment.strand = '-'; else fragment.strand = '+';
			if (chrom >= 0 && PEONLY == 0) {

				fragment.id = read->readid;
				fragment.barcode = read->barcode;
				if (REALIGN_VARIANTS){
					realign_and_extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
				}else{
					extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
				}
				if (fragment.variants >=2) VOfragments[0]++;
				else if (fragment.variants >=1) VOfragments[1]++;
				if (fragment.variants >= 2 || (SINGLEREADS == 1 && fragment.variants >= 1)) print_fragment(&fragment, varlist, fragment_file);
			}
		} else // paired-end read
		{
			fragment.variants = 0;
			fragment.id = read->readid; //v1 =0; v2=0;
			fragment.barcode = read->barcode;

			if (chrom >=0) extract_variants_read(read,ht,chromvars,varlist,1,&fragment,chrom,reflist);

			//fprintf(stderr,"paired read stats %s %d flag %d IS %d\n",read->chrom,read->cigs,read->flag,read->IS);
			if (fragment.variants > 0) {
				//fprintf(stderr,"variants %d read %s %s \n",fragment.variants,read->chrom,read->readid);
				add_fragment(flist, &fragment, read, fragments);
				fragments++;
				if (fragments >= MAXFRAG) {
					fprintf(stderr, "exceeded max #cached fragments: %d,increase MAXFRAGMENTS using --maxfragments option \n", MAXFRAG);
					return -1;
				}
			}
		}
		if ((fragments - prevfragments >= 100000) || fragments >= MAXFRAG - 10000 || (chrom != prevchrom && prevchrom != -1 && fragments > 0)) // chrom of current read is not the same as previous read's chromosome...
		{
			if (PFLAG == 1 && chrom == prevchrom) fprintf(stderr, "cleaning buffer: current chrom %s position %d fragments %d\n", read->chrom, read->position, fragments);
			else if (PFLAG == 1 && chrom != prevchrom) fprintf(stderr, "cleaning buffer for prev chrom fragments %d\n", fragments);
			if (fragments > 0) clean_fragmentlist(flist, &fragments, varlist, chrom, read->position, prevchrom);
			prevfragments = fragments;
			//fprintf(stderr,"remaining %d\n",fragments);
		}

		reads += 1;
		if ( (reads % 2000000 == 0 && chrom >=0) || (PACBIO ==1 && reads % 20000==0 && chrom >= 0) ) { 
			fprintf(stderr, "processed %d reads", reads);
			if (DATA_TYPE != 2) fprintf(stderr, " PE-fragments %d SE-fragments %d %d", fragments,VOfragments[0],VOfragments[1]);
			fprintf(stderr, "\n");
		}
		prevchrom = chrom;
		prevtid = read->tid;
		free_readmemory(read);
	} // end of main while loop
	if (fragments > 0)  // still fragments left in buffer
	{
		fprintf(stderr, "final cleanup of fragment list: %d current chrom %d prev %d \n", fragments, chrom,prevchrom);
		if (prevchrom >=0) clean_fragmentlist(flist, &fragments, varlist, -1, read->position, prevchrom); // added extra filter 03/08/18
	}
	bam_destroy1(b);

	free(flist); free(read); free(fragment.alist);
	if (REALIGN_VARIANTS) free(fcigarlist);
	return 0;
}

int main(int argc, char** argv) {

	char samfile[1024];
	char bamfile[1024];
	char variantfile[1024];
	char fastafile[1024];
	char* regions= NULL;  
	strcpy(samfile, "None");
	strcpy(bamfile, "None");
	strcpy(variantfile, "None");
	strcpy(fastafile, "None");
	int readsorted = 0;
	char* sampleid = (char*) malloc(1024);
	sampleid[0] = '-';
	sampleid[1] = '\0';
	int samplecol = 10; // default if there is a single sample in the VCF file
	int i = 0, variants = 0, hetvariants = 0;
	char** bamfilelist = NULL;
	int bamfiles = 0;

	logfile = NULL;
	fragment_file = stdout; // write fragments to this file if it is present

	if (argc % 2 != 1){
		fprintf(stderr, "\nERROR: Invalid number of arguments specified.\n");
		exit(1);
	}

	AP = init_params();

	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "--bam") == 0 || strcmp(argv[i], "--bamfile") == 0) bamfiles++;
		else if (strcmp(argv[i], "--variants") == 0) strcpy(variantfile, argv[i + 1]);
		else if (strcmp(argv[i], "--reffile") == 0 || strcmp(argv[i], "--ref") == 0) strcpy(fastafile, argv[i + 1]);
		else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
			strcpy(variantfile, argv[i + 1]);
			VCFformat = 1;
			int l=strlen(variantfile); 
			if (variantfile[l-1] == 'z' && variantfile[l-2] == 'g' && variantfile[l-3] == '.') 
			{
				fprintf(stderr,"The input VCF file appears to be gzipped (.gz extension), hapcut only accepts unzipped VCF files as input\n Please provide an unzipped VCF file or make sure that the file doesn't have the .gz extension\n\n");
				exit(0);
			}
		} else if (strcmp(argv[i], "--sorted") == 0){
			check_input_0_or_1(argv[i + 1]);
			readsorted = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i],"--region") ==0 || strcmp(argv[i],"--regions") ==0) { // added 11/30/17
			regions = (char*)malloc(strlen(argv[i+1])+1); strcpy(regions,argv[i+1]); 
		}
		else if (strcmp(argv[i], "--mmq") == 0) MIN_MQ = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--HiC") == 0 || strcmp(argv[i], "--hic") == 0){
			check_input_0_or_1(argv[i + 1]);
			if (atoi(argv[i + 1])){
				MAX_IS = 40000000;
				NEW_FORMAT = 1;
				DATA_TYPE = 1;
			}
		}
		else if (strcmp(argv[i], "--10X") == 0 || strcmp(argv[i], "--10x") == 0){
			check_input_0_or_1(argv[i + 1]);
			if (atoi(argv[i + 1])){
				SINGLEREADS = 1;
				NEW_FORMAT = 1;
				DATA_TYPE = 2;
			}
		}else if (strcmp(argv[i], "--realign_variants") == 0){
			check_input_0_or_1(argv[i + 1]);
			if (atoi(argv[i + 1])){
				REALIGN_VARIANTS = 1;
			}
		}else if (strcmp(argv[i], "--pacbio") == 0 || strcmp(argv[i], "--SMRT") == 0 || strcmp(argv[i],"--pb") ==0){
			check_input_0_or_1(argv[i + 1]);
			if (atoi(argv[i + 1])){
				REALIGN_VARIANTS = 1; PACBIO =1; MINQ = 4;
				SUM_ALL_ALIGN = 1; 
			}

			// scores based on https://www.researchgate.net/figure/230618348_fig1_Characterization-of-Pacific-Biosciences-dataa-Base-error-mode-rate-for-deletions
			MATCH = log10(1.0 - (0.01 + 0.071 + 0.037)); //log10(1.0 - (0.01 + 0.12 + 0.02));
			MISMATCH = log10(0.01);
			INSERTION_OPEN = log10(0.071); //log10(0.12);
			DELETION_OPEN = log10(0.037); //log10(0.02);
			// estimated these empirically from bam file
			INSERTION_EXTEND = log10(0.26);
			DELETION_EXTEND = log10(0.12);


		}else if (strcmp(argv[i], "--ont") == 0 || strcmp(argv[i], "--ONT") == 0){
			check_input_0_or_1(argv[i + 1]);
			if (atoi(argv[i + 1])){
				REALIGN_VARIANTS = 1;  MINQ = 4;
			}
			// scores based on http://www.nature.com/nmeth/journal/v12/n4/abs/nmeth.3290.html
			MATCH = log10(1.0 - (0.051 + 0.049 + 0.078));
			MISMATCH = log10(0.051);
			INSERTION_OPEN = log10(0.049);
			INSERTION_EXTEND = log10(0.25); // this number has no basis in anything
			DELETION_OPEN = log10(0.078);
			DELETION_EXTEND = log10(0.25); // this number also has no basis in anything
			SUM_ALL_ALIGN = 1; 

		}else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "--v") == 0){
			check_input_0_or_1(argv[i + 1]);
			VERBOSE = atoi(argv[i + 1]);
		}else if (strcmp(argv[i], "--new_format") == 0 || strcmp(argv[i], "--nf") == 0){
			check_input_0_or_1(argv[i + 1]);
			NEW_FORMAT = atoi(argv[i + 1]);
		}else if (strcmp(argv[i], "--maxIS") == 0)
			MAX_IS = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--minIS") == 0)
			MIN_IS = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--fullprint") == 0)
			PRINT_COMPACT = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--PEonly") == 0){
			check_input_0_or_1(argv[i + 1]);
			PEONLY = atoi(argv[i + 1]); // discard single end mapped reads
		}else if (strcmp(argv[i], "--indels") == 0){
			check_input_0_or_1(argv[i + 1]);
			PARSEINDELS = atoi(argv[i + 1]); // allow indels in hairs
		}else if (strcmp(argv[i], "--pflag") == 0){
			check_input_0_or_1(argv[i + 1]);
			IFLAG = atoi(argv[i + 1]); // allow indels in hairs
		}else if (strcmp(argv[i], "--qvoffset") == 0) QVoffset = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0)
			fragment_file = fopen(argv[i + 1], "w");
		else if (strcmp(argv[i], "--logfile") == 0 || strcmp(argv[i], "--log") == 0)
			logfile = fopen(argv[i + 1], "w");
		else if (strcmp(argv[i], "--singlereads") == 0){
			check_input_0_or_1(argv[i + 1]);
			SINGLEREADS = atoi(argv[i + 1]);
		}else if (strcmp(argv[i], "--maxfragments") == 0) MAXFRAG = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--hom") == 0) 
		{
			check_input_0_or_1(argv[i + 1]);
			HOMOZYGOUS = atoi(argv[i + 1]);
			if (HOMOZYGOUS ==1) SINGLEREADS = 1; 
		}
		else if (strcmp(argv[i], "--mbq") == 0) MINQ = atoi(argv[i + 1]);
		else if (strcmp(argv[i], "--noquality") == 0){
			check_input_0_or_1(argv[i + 1]);
			MISSING_QV = atoi(argv[i + 1]);
		}else if (strcmp(argv[i], "--triallelic") == 0){
			check_input_0_or_1(argv[i + 1]);
			TRI_ALLELIC = atoi(argv[i + 1]);
		}else if (strcmp(argv[i], "--fosmids") == 0 || strcmp(argv[i], "--fosmid") == 0){
			check_input_0_or_1(argv[i + 1]);
			LONG_READS = 1;
		}else if (strcmp(argv[i], "--sumall") == 0) { 
			SUM_ALL_ALIGN = atoi(argv[i+1]); 
			if (SUM_ALL_ALIGN >=1) fprintf(stderr, "\nusing sum of all alignments for scoring \n");
		}else if (strcmp(argv[i], "--ep") == 0) { 
			ESTIMATE_PARAMS = atoi(argv[i+1]); 
		}else{
			fprintf(stderr, "\nERROR: Invalid Option \"%s\" specified.\n",argv[i]);
			exit(1);
		}
	}
	if (REALIGN_VARIANTS && strcmp(fastafile, "None") == 0) {
		fprintf(stderr, "\nERROR: In order to realign variants (including --pacbio and --ont options), reference fasta file must be provided with --ref option.\n");
		exit(1);
	}
	if (MINQ < 4) {
		fprintf(stderr, "\nERROR: MINQ must be at least 4.\n");
		exit(1);
	}
	if (bamfiles > 0 && strcmp(variantfile, "None") != 0) {
		bamfilelist = (char**) malloc(sizeof (char*)*bamfiles);
		for (i = 0; i < bamfiles; i++) bamfilelist[i] = (char*) malloc(1024);
		bamfiles = 0;
		for (i = 1; i < argc; i += 2) {
			if (strcmp(argv[i], "--bam") == 0 || strcmp(argv[i], "--bamfile") == 0) strcpy(bamfilelist[bamfiles++], argv[i + 1]);
		}
		fprintf(stderr, "\nExtracting haplotype informative reads from bamfiles %s minQV %d minMQ %d maxIS %d \n\n", bamfilelist[0], MINQ, MIN_MQ, MAX_IS);
	} else {
		print_options();
		return -1;
	}

	HASHTABLE ht;
	ht.htsize = 7919;
	init_hashtable(&ht); // chromosome names are inserted into hashtable from VCF file, ideally this should be done using BAM file header to avoid missing some contigs/chroms
	VARIANT* varlist;
	int chromosomes = 0;

	variants = count_variants(variantfile, sampleid, &samplecol);
	if (variants < 0) return -1;
	varlist = (VARIANT*) malloc(sizeof (VARIANT) * variants);
	chromosomes = read_variantfile(variantfile, varlist, &ht, &hetvariants, samplecol);
	// variants is set to hetvariants only, but this is not correct since
	VARIANTS = variants;
	// there are two options, we include all variants in the chromvars datastructure but only use heterozygous variants for outputting HAIRS
	// variant-id should correspond to line-number in VCF file since that will be used for printing out variants in Hapcut

	//	fprintf(stderr,"read %d variants from file %s chromosomes %d\n",snps,argv[1],chromosomes);
	CHROMVARS* chromvars = (CHROMVARS*) malloc(sizeof (CHROMVARS) * chromosomes);
	build_intervalmap(chromvars, chromosomes, varlist, VARIANTS);

	// read reference fasta file for INDELS, reads entire genome in one shot but saves memory by only keep contigs that are relevant for VCF parsing
	REFLIST* reflist = (REFLIST*) malloc(sizeof (REFLIST));
	reflist->ns = 0; reflist->names = NULL; reflist->lengths = NULL; reflist->sequences = NULL; reflist->current = -1;
	if (strcmp(fastafile, "None") != 0) {
		if (read_fastaheader(fastafile, reflist) > 0) {
			reflist->sequences = calloc(reflist->ns, sizeof (char*)); //(char**)malloc(sizeof(char*)*reflist->ns);
			for (i = 0; i < reflist->ns; i++) {
				reflist->sequences[i] = calloc(reflist->lengths[i] + 1, sizeof (char));
				int chrom = getindex(&ht, reflist->names[i]); reflist->used[i] = 1;
				if (chrom >=0) fprintf(stderr,"found match for reference contig %s in VCF file index \n",reflist->names[i]); 
				else reflist->used[i] = 0; // memory for this chromosome will be freed 
				if (i < 5) fprintf(stderr, "contig %s length %d\n", reflist->names[i], reflist->lengths[i]);
			}
			read_fasta(fastafile, reflist);
		}
	}
	if (readsorted == 0 && bamfiles > 0) {
		for (i = 0; i < bamfiles; i++) {
			int parse_ok = 0;
			if (LONG_READS == 0) parse_ok = parse_bamfile_sorted(bamfilelist[i], &ht, chromvars, varlist, reflist,regions);
			//else parse_ok = parse_bamfile_fosmid(bamfilelist[i], &ht, chromvars, varlist, reflist); // fosmid pool bam file
			if (parse_ok != 0) return parse_ok;
		}
	}

	if (logfile != NULL) fclose(logfile);
	if (fragment_file != NULL && fragment_file != stdout) fclose(fragment_file);

	for (i=0;i<reflist->ns;i++){
		free(reflist->names[i]);
		if (reflist->used[i] ==1) free(reflist->sequences[i]);
	}
	if (reflist->ns > 0) { 
		free(reflist->names); 
		free(reflist->sequences); 
		free(reflist->lengths); 
		free(reflist->used);
	}
	//free(reflist->offsets);
	//int xor = pow(2,16)-1; int c=0;
	for (i=0;i<variants;i++){
		//fprintf(stderr,"variant %d %s %d cov %d %s %s ",i+1,varlist[i].genotype,varlist[i].position,varlist[i].depth,varlist[i].RA,varlist[i].AA);
		//fprintf(stderr,"REF(strand) %d:%d ALT %d:%d\n",varlist[i].A1>>16,varlist[i].A1 & xor,varlist[i].A2>>16,varlist[i].A2 & xor);
		free(varlist[i].genotype); free(varlist[i].RA); 	free(varlist[i].AA);free(varlist[i].chrom);
		if (varlist[i].heterozygous == '1'){
			free(varlist[i].allele1); free(varlist[i].allele2);
		}
	}

	for (i=0;i<chromosomes;i++){
		free(chromvars[i].intervalmap);
	}
	free(chromvars);
	free(sampleid); free(varlist); free(reflist);
	if (bamfiles > 0 && strcmp(variantfile,"None") !=0){
		for (i=0;i<bamfiles;i++)
			free(bamfilelist[i]);
		free(bamfilelist);
	}

	return 0;
}
