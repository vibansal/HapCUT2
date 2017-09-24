// author: Vikas Bansal, vbansal@scripps, 2011-2012
//  program to extract haplotype informative reads from sorted BAM file, input requirements: bamfile and variantfile
//  Jan 13 2012, changed to read directly from BAM file
//  paired-end overlapping reads need to be handled properly
//  add module for RNA-seq data as well
//  add flag for secondary alignments and PCR duplicates (picard mark duplicates)
//  input format for variants should be VCF from now
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
#include "sam.h"
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
//int QVoffset = 33; declared in samread.h
FILE* logfile;
int PFLAG = 1;
int PRINT_FRAGMENTS = 1;
char* GROUPNAME; // for fragments from different pools, SRRxxx
FILE* fragment_file;
int TRI_ALLELIC = 0;
int VERBOSE = 0;

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
int DATA_TYPE = 0;

// NEW FORMAT FOR HIC
// same as old format but with more columns.
// column 3 is the data type (0 for normal, 1 for HiC)
// column 4 is the first SNP index of mate 2 (-1 if no mate 2)
// column 5 is the absolute insert size
int NEW_FORMAT = 0;

//int get_chrom_name(struct alignedread* read,HASHTABLE* ht,REFLIST* reflist);

#include "parsebamread.c"
#include "realignbamread.c"
#include "fosmidbam_hairs.c" // code for parsing fosmid pooled sequence data

//disabled sam file reading
//#include "samhairs.c" // has two functions that handle sam file parsing

void print_options();
int parse_bamfile_sorted(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist);

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
    fprintf(stderr, "--VCF <FILENAME> : variant file with genotypes for a single individual in VCF format\n");
    fprintf(stderr, "--variants : variant file in hapCUT format (use this option or --VCF option but not both), this option will be phased out in future releases\n");
    fprintf(stderr, "--maxIS <INT> : maximum insert size for a paired-end read to be considered as a single fragment for phasing, default 1000\n");
    fprintf(stderr, "--minIS <INT> : minimum insert size for a paired-end read to be considered as single fragment for phasing, default 0\n");
    fprintf(stderr, "--PEonly <0/1> : do not use single end reads, default is 0 (use all reads)\n");
    fprintf(stderr, "--indels <0/1> : extract reads spanning INDELS, default is 0, variants need to specified in VCF format to use this option\n");
    fprintf(stderr, "--noquality <INTEGER> : if the bam file does not have quality string, this value will be used as the uniform quality value, default 0 \n");
    //fprintf(stderr,"--triallelic <0/1> : print information about , default 0 \n");
    fprintf(stderr, "--ref <FILENAME> : reference sequence file (in fasta format), optional but required for indels, should be indexed using samtools\n");
    fprintf(stderr, "--out <FILENAME> : output filename for haplotype fragments, if not provided, fragments will be output to stdout\n\n");
    //fprintf(stderr,"--out : output file for haplotype informative fragments (hairs)\n\n");
}

void check_input_0_or_1(char* x){
    if (!(strcmp(x, "0") == 0 || strcmp(x, "1") == 0)){
        fprintf(stderr, "\nERROR: Invalid input \"%s\" for <0/1> option flag.\n",x);
        exit(1);
    }
}

// extract haplotype informative reads from sorted bam file //
// need to discard reads that are marked as duplicates using flag //

int parse_bamfile_sorted(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist) {
    fprintf(stderr, "reading sorted bamfile %s \n", bamfile);
    int reads = 0;
    struct alignedread* read = (struct alignedread*) malloc(sizeof (struct alignedread));

    if (REALIGN_VARIANTS){
        fcigarlist = calloc(sizeof(int),400000);
    }
    int i = 0;
    int chrom = 0; //int sl=0;
    // int v1,v2;
    int prevchrom = -1;
    int prevtid = -1;

    FRAGMENT* flist = (FRAGMENT*) malloc(sizeof (FRAGMENT) * MAXFRAG);
    int fragments = 0;
    int prevfragments = 0;
    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*10000);

    samfile_t *fp;
    if ((fp = samopen(bamfile, "rb", 0)) == 0) {
        fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
        return -1;
    }
    bam1_t *b = bam_init1();

    while (samread(fp, b) >= 0) {
        fetch_func(b, fp, read);
        if ((read->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) || read->mquality < MIN_MQ) {
            free_readmemory(read);
            continue;
        }
        // find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
        if (read->tid != prevtid) {
            chrom = getindex(ht, read->chrom); // doing this for every read, can replace this by string comparison ..april 4 2012
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

        //if (read->tid == read->mtid) // use mateposition to calculate insert size, march 12 2013, wrong since we need to consider the readlength/cigar
        //{
        //    read->IS = read->mateposition - read->position;
        //}

        fragment.absIS = (read->IS < 0) ? -1 * read->IS : read->IS;
        // add check to see if the mate and its read are on same chromosome, bug for contigs, july 16 2012
        if ((read->flag & 8) || fragment.absIS > MAX_IS || fragment.absIS < MIN_IS || read->IS == 0 || !(read->flag & 1) || read->tid != read->mtid) // single read
        {
            fragment.variants = 0; // v1 =0; v2=0;
            if (chrom >= 0 && PEONLY == 0) {

                fragment.id = read->readid;
                fragment.barcode = read->barcode;
                if (REALIGN_VARIANTS){
                    realign_and_extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
                }else{
                    extract_variants_read(read,ht,chromvars,varlist,0,&fragment,chrom,reflist);
                }
                if (fragment.variants >= 2 || (SINGLEREADS == 1 && fragment.variants >= 1)) {
                    // instead of printing fragment, we could change this to update genotype likelihoods
                    print_fragment(&fragment, varlist, fragment_file);
                }
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
        // BUG here when the fragment list cannot be cleaned due to long mate-pair fragments (accumulated for large IS)
        // fragments >= 100000 and we will clean it repeatedly...
        // need to fix this june 4 2012.... even for long mate-pairs this could be a problem...
        if ((fragments - prevfragments >= 100000) || fragments >= MAXFRAG - 10000 || (chrom != prevchrom && prevchrom != -1 && fragments > 0)) // chrom of current read is not the same as previous read's chromosome...
        {
            if (PFLAG == 1) fprintf(stderr, "cleaning buffer: current chrom %s %d fragments %d\n", read->chrom, read->position, fragments);
            // BUG HERE when trying to clean empty fragment list (fragments ==0)
            if (fragments > 0) clean_fragmentlist(flist, &fragments, varlist, chrom, read->position, prevchrom);
            prevfragments = fragments;
            //fprintf(stderr,"remaining %d\n",fragments);
        }

        reads += 1;
        if (reads % 2000000 == 0){
            fprintf(stderr, "processed %d reads", reads);
            if (DATA_TYPE != 2) fprintf(stderr, ", paired end fragments %d", fragments);
            fprintf(stderr, "\n");
        }
        prevchrom = chrom;
        prevtid = read->tid;
        free_readmemory(read);
    }
    if (fragments > 0) {
        fprintf(stderr, "final cleanup of fragment list: %d current chrom %s %d \n", fragments, read->chrom, read->position);
        clean_fragmentlist(flist, &fragments, varlist, -1, read->position, prevchrom);
    }
    bam_destroy1(b);

    free(flist); free(read); free(fragment.alist);
    if (REALIGN_VARIANTS){
        free(fcigarlist);
    }
    return 0;
}

int main(int argc, char** argv) {
    char samfile[1024];
    char bamfile[1024];
    char variantfile[1024];
    char fastafile[1024];
    strcpy(samfile, "None");
    strcpy(bamfile, "None");
    strcpy(variantfile, "None");
    strcpy(fastafile, "None");
    GROUPNAME = NULL;
    int readsorted = 0;
    char* sampleid = (char*) malloc(1024);
    sampleid[0] = '-';
    sampleid[1] = '\0';
    int samplecol = 10; // default if there is a single sample in the VCF file
    int i = 0, j = 0, variants = 0, hetvariants = 0;
    char** bamfilelist = NULL;
    int bamfiles = 0;

    logfile = NULL;
    fragment_file = stdout; // write fragments to this file if it is present

    if (argc % 2 != 1){
        fprintf(stderr, "\nERROR: Invalid number of arguments specified.\n");
        exit(1);
    }

    for (i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "--bam") == 0 || strcmp(argv[i], "--bamfile") == 0) bamfiles++;
        else if (strcmp(argv[i], "--variants") == 0) strcpy(variantfile, argv[i + 1]);
        else if (strcmp(argv[i], "--reffile") == 0 || strcmp(argv[i], "--ref") == 0) strcpy(fastafile, argv[i + 1]);
        else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
            strcpy(variantfile, argv[i + 1]);
            VCFformat = 1;
        } else if (strcmp(argv[i], "--sorted") == 0){
            check_input_0_or_1(argv[i + 1]);
            readsorted = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--mbq") == 0) MINQ = atoi(argv[i + 1]);
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
        }else if (strcmp(argv[i], "--pacbio") == 0){
            check_input_0_or_1(argv[i + 1]);
            if (atoi(argv[i + 1])){
                REALIGN_VARIANTS = 1;
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
                REALIGN_VARIANTS = 1;
            }
            // scores based on http://www.nature.com/nmeth/journal/v12/n4/abs/nmeth.3290.html
            MATCH = log10(1.0 - (0.051 + 0.049 + 0.078));
            MISMATCH = log10(0.051);
            INSERTION_OPEN = log10(0.049);
            INSERTION_EXTEND = log10(0.25); // this number has no basis in anything
            DELETION_OPEN = log10(0.078);
            DELETION_EXTEND = log10(0.25); // this number also has no basis in anything

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
        else if (strcmp(argv[i], "--noquality") == 0){
            check_input_0_or_1(argv[i + 1]);
            MISSING_QV = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--triallelic") == 0){
            check_input_0_or_1(argv[i + 1]);
            TRI_ALLELIC = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--fosmids") == 0 || strcmp(argv[i], "--fosmid") == 0){
            check_input_0_or_1(argv[i + 1]);
            LONG_READS = 1;
        }else if (strcmp(argv[i], "--groupname") == 0) {
            GROUPNAME = (char*) malloc(1024);
            strcpy(GROUPNAME, argv[i + 1]);
        }else{
            fprintf(stderr, "\nERROR: Invalid Option \"%s\" specified.\n",argv[i]);
            exit(1);
        }
    }
    if (REALIGN_VARIANTS && strcmp(fastafile, "None") == 0) {
        fprintf(stderr, "\nERROR: In order to realign variants (including --pacbio and --ont options), reference fasta file must be provided with --ref option.\n");
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
    init_hashtable(&ht);
    VARIANT* varlist;
    int chromosomes = 0;

    if (VCFformat == 1) {
        variants = count_variants(variantfile, sampleid, &samplecol);
        if (variants < 0) return -1;
        varlist = (VARIANT*) malloc(sizeof (VARIANT) * variants);
        chromosomes = read_variantfile(variantfile, varlist, &ht, &hetvariants, samplecol);
    } else {
        variants = count_variants_oldformat(variantfile);
        if (variants < 0) return -1;
        varlist = (VARIANT*) malloc(sizeof (VARIANT) * variants);
        chromosomes = read_variantfile_oldformat(variantfile, varlist, &ht, variants);
    }
    // variants is set to hetvariants only, but this is not correct since
    VARIANTS = variants;
    // there are two options, we include all variants in the chromvars datastructure but only use heterozygous variants for outputting HAIRS
    // variant-id should correspond to line-number in VCF file since that will be used for printing out variants in Hapcut

    //	fprintf(stderr,"read %d variants from file %s chromosomes %d\n",snps,argv[1],chromosomes);
    CHROMVARS* chromvars = (CHROMVARS*) malloc(sizeof (CHROMVARS) * chromosomes);
    build_intervalmap(chromvars, chromosomes, varlist, VARIANTS);

    // read reference fasta file for INDELS, currently reads entire genome in one go, need to modify to read chromosome by chromosome
    REFLIST* reflist = (REFLIST*) malloc(sizeof (REFLIST));
    reflist->ns = 0;
    reflist->names = NULL;
    reflist->lengths = NULL;
    reflist->sequences = NULL;
    reflist->current = -1;
    if (strcmp(fastafile, "None") != 0) {
        if (read_fastaheader(fastafile, reflist) > 0) {
            reflist->sequences = calloc(reflist->ns, sizeof (char*)); //(char**)malloc(sizeof(char*)*reflist->ns);
            for (i = 0; i < reflist->ns; i++) {
                reflist->sequences[i] = calloc(reflist->lengths[i] + 1, sizeof (char));
                if (i < 5) fprintf(stderr, "contig %s length %d\n", reflist->names[i], reflist->lengths[i]);
            }
            read_fasta(fastafile, reflist);
        }
    }
    if (readsorted == 0 && bamfiles > 0) {
        for (i = 0; i < bamfiles; i++) {
			int parse_ok = 0;
            if (LONG_READS == 0) parse_ok = parse_bamfile_sorted(bamfilelist[i], &ht, chromvars, varlist, reflist);
            else parse_ok = parse_bamfile_fosmid(bamfilelist[i], &ht, chromvars, varlist, reflist); // fosmid pool bam file
			if (parse_ok != 0) return parse_ok;
        }
    }

    if (logfile != NULL) fclose(logfile);
    if (fragment_file != NULL && fragment_file != stdout) fclose(fragment_file);

    for (i=0;i<reflist->ns;i++){
		free(reflist->names[i]);
		free(reflist->sequences[i]);
	}
	free(reflist->names);
	free(reflist->sequences);
	free(reflist->lengths);
	//free(reflist->offsets);

	for (i=0;i<variants;i++){
		free(varlist[i].genotype);
		free(varlist[i].allele1);
		free(varlist[i].allele2);
		free(varlist[i].RA);
		free(varlist[i].AA);
		free(varlist[i].chrom);
	}

	for (i=0;i<chromosomes;i++){
		free(chromvars[i].intervalmap);
	}
	free(chromvars);

    for (i=0;i<ht.htsize;i++){
        if (ht.blist[i] != NULL){
            for(j=0;j<ht.bucketlengths[i];j++){
                free(ht.blist[i]->key);
            }
            free(ht.blist[i]);
        }
    }
    free(ht.blist);
    free(ht.bucketlengths);

	free(sampleid); free(varlist); free(reflist);
	if (bamfiles > 0 && strcmp(variantfile,"None") !=0){
		for (i=0;i<bamfiles;i++)
			free(bamfilelist[i]);
		free(bamfilelist);
	}

    return 0;
}
