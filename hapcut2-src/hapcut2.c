// CODE STARTED SEPT 10 2007 4pm //  april 8 2008 this code used for producing results in ECCB 2008 paper //
// author: VIKAS BANSAL (vbansal@scripps.edu) last modified December 23, 2010

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "common.h"
#include "fragmatrix.h"
#include "pointerheap.h"
#include "readinputfiles.h"
#include "printhaplotypes.c"

// Printing related
int VERBOSE = 0;
int PRINT_FRAGMENT_SCORES = 0; // output the MEC/switch error score of erroneous reads/fragments to a file for evaluation

// Quality-score related parameters
int QVoffset = 33;
int MINQ = 6; // additional base quality filter in hapcut added april 18 2012

// Number of iterations
int MAXITER = 10000;     // maximum number of global iterations
int MAXCUT_ITER = 10000; // maximum number of iterations for max-cut algorithm, if this is proportional to 'N' -> complexity is 'N^2', added march 13 2013
int CONVERGE = 5; // stop iterations on a given block if exceed this many iterations since improvement

// Post-processing related variables
float THRESHOLD = 0.8;
float SPLIT_THRESHOLD = 0.8;
float HOMOZYGOUS_PRIOR = -80; // in log form. assumed to be really unlikely
int CALL_HOMOZYGOUS = 0;
int SPLIT_BLOCKS = 0;
int DISCRETE_PRUNING = 0;
int ERROR_ANALYSIS_MODE = 0;
int SKIP_PRUNE = 0;
int SNVS_BEFORE_INDELS = 0;

int AUTODETECT_LONGREADS = 1;
int LONG_READS = 0; // if this variable is 1, the data contains long read data

// HiC-related global variables
int HIC = 0;
int MAX_HIC_EM_ITER = 1;
int NEW_FRAGFILE_FORMAT = 0;
int HTRANS_BINSIZE = 5000;
int HTRANS_MAXBINS = 10000; // this value will be overwritten at startup
int HTRANS_READ_LOWBOUND = 500;
int HTRANS_MAX_WINDOW = 4000000; // maximum window size for h-trans estimation
char HTRANS_DATA_INFILE[10000];
char HTRANS_DATA_OUTFILE[10000];
int MAX_IS = -1;

#include "find_maxcut.c"   // function compute_good_cut
#include "post_processing.c"  // post-processing functions

int maxcut_haplotyping(char* fragmentfile, char* variantfile, char* outputfile) {
    // IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+

    fprintf_time(stderr, "Calling Max-Likelihood-Cut based haplotype assembly algorithm\n");

    int snps = 0;
    int fragments = 0, iter = 0, components = 0;
    int i = 0, k = 0;
    int* slist;
    int flag = 0;
    float bestscore = 0, miscalls = 0;
    char buffer[MAXBUF];
    int hic_iter=0;
    struct SNPfrags* snpfrag = NULL;
    struct BLOCK* clist;
    char* HAP1;
    float HIC_LL_SCORE = -80;
    float OLD_HIC_LL_SCORE = -80;
    int converged_count=0, split_count, new_components, component;

    int new_fragments = 0;
    struct fragment* new_Flist;

    // READ FRAGMENT MATRIX
    struct fragment* Flist;
    FILE* ff = fopen(fragmentfile, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open fragment file %s\n", fragmentfile);
        exit(0);
    }
    fragments = 0;
    while (fgets(buffer, MAXBUF, ff) != NULL){
        if (!((buffer[0] == '0')&&(buffer[1] == ' ')))
            fragments++;
    }
    fclose(ff);
    Flist     = (struct fragment*) malloc(sizeof (struct fragment)* fragments);

    flag = read_fragment_matrix(fragmentfile, Flist, fragments);

    if (MAX_IS != -1){
        // we are going to filter out some insert sizes
        new_fragments = 0;
        new_Flist = (struct fragment*) malloc(sizeof (struct fragment)* fragments);
        for(i = 0; i < fragments; i++){

            if (Flist[i].isize < MAX_IS){
                new_Flist[new_fragments] = Flist[i];
                new_fragments++;
            }
        }

        Flist = new_Flist;
        fragments = new_fragments;
    }

    if (flag < 0) {
        fprintf_time(stderr, "unable to read fragment matrix file %s \n", fragmentfile);
        return -1;
    }

    //ADD EDGES BETWEEN SNPS
    snps = count_variants_vcf(variantfile);

    if (snps < 0) {
        fprintf_time(stderr, "unable to read variant file %s \n", variantfile);
        return -1;
    }

    float mean_snps_per_read = 0;
    if (AUTODETECT_LONGREADS){
        for (i = 0; i < fragments; i++){
            mean_snps_per_read += Flist[i].calls;
        }
        mean_snps_per_read /= fragments;
        if (mean_snps_per_read >= 3){
            LONG_READS = 1;
        }else{
            LONG_READS = 0;
        }
    }

    snpfrag = (struct SNPfrags*) malloc(sizeof (struct SNPfrags)*snps);
    update_snpfrags(Flist, fragments, snpfrag, snps, &components);

    // 10/25/2014, edges are only added between adjacent nodes in each fragment and used for determining connected components...
    for (i = 0; i < snps; i++) snpfrag[i].elist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));
    if (LONG_READS ==0){
        add_edges(Flist,fragments,snpfrag,snps,&components);
    }else if (LONG_READS >=1){
        add_edges_fosmids(Flist,fragments,snpfrag,snps,&components);
    }

    for (i = 0; i < snps; i++) snpfrag[i].telist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));

    // this considers only components with at least two nodes
    fprintf_time(stderr, "fragments %d snps %d component(blocks) %d\n", fragments, snps, components);

    // BUILD COMPONENT LIST
    clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*components);
    generate_clist_structure(Flist, fragments, snpfrag, snps, components, clist);

    // READ VCF FILE
    read_vcffile(variantfile, snpfrag, snps);

    // INITIALIZE RANDOM HAPLOTYPES
    HAP1 = (char*) malloc(snps + 1);
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].frags == 0 || (SNVS_BEFORE_INDELS && (strlen(snpfrag[i].allele0) != 1 || strlen(snpfrag[i].allele1) != 1))) {
            HAP1[i] = '-';
        } else if (drand48() < 0.5) {
            HAP1[i] = '0';
        } else {
            HAP1[i] = '1';
        }
    }

    // for each block, we maintain best haplotype solution under MFR criterion
    // compute the component-wise score for 'initHAP' haplotype
    miscalls = 0;
    bestscore = 0;
    for (k = 0; k < components; k++) {
        clist[k].SCORE = 0;
        clist[k].bestSCORE = 0;
        for (i = 0; i < clist[k].frags; i++) {
            update_fragscore(Flist, clist[k].flist[i], HAP1);
            clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
        }
        clist[k].bestSCORE = clist[k].SCORE;
        bestscore += clist[k].bestSCORE;
        miscalls += clist[k].SCORE;
    }

    fprintf_time(stderr, "processed fragment file and variant file: fragments %d variants %d\n", fragments, snps);

    int MAXIS = -1;

    if (HIC){

        // determine the probability of an h-trans interaction for read

        for (i=0; i<fragments;i++){

            Flist[i].htrans_prob = -80;

            if (Flist[i].isize > MAXIS)
                MAXIS = Flist[i].isize;
        }

        HTRANS_MAXBINS = MAXIS/HTRANS_BINSIZE + 1;
    }else{
        HTRANS_MAXBINS = 0;
    }

    // read in file with estimated probabilities of Hi-C h-trans interactions with distance
    if (strcmp(HTRANS_DATA_INFILE, "None") != 0){
        int num_bins        = count_htrans_bins(HTRANS_DATA_INFILE);
        float* htrans_probs = (float*) malloc(sizeof(float) * num_bins);
        read_htrans_file(HTRANS_DATA_INFILE, htrans_probs, num_bins);
        for (i=0; i<fragments;i++){
            Flist[i].htrans_prob = log10(htrans_probs[Flist[i].isize / HTRANS_BINSIZE]);
        }
        free(htrans_probs);
    }

    slist = (int*) malloc(sizeof (int)*snps);

    OLD_HIC_LL_SCORE = bestscore;
    for (hic_iter = 0; hic_iter < MAX_HIC_EM_ITER; hic_iter++){
        if (VERBOSE)
            fprintf_time(stdout, "HIC ITER %d\n", hic_iter);
        for (k = 0; k < components; k++){
            clist[k].iters_since_improvement = 0;
        }
        for (i=0; i<snps; i++){
            snpfrag[i].post_hap = 0;
        }
        // RUN THE MAX_CUT ALGORITHM ITERATIVELY TO IMPROVE LIKELIHOOD

        for (iter = 0; iter < MAXITER; iter++) {
            if (VERBOSE)
                fprintf_time(stdout, "PHASING ITER %d\n", iter);
            converged_count = 0;
            for (k = 0; k < components; k++){
                if(VERBOSE && iter == 0)
                    fprintf_time(stdout, "component %d length %d phased %d %d...%d\n", k, clist[k].length, clist[k].phased, clist[k].offset, clist[k].lastvar);
                if (clist[k].SCORE > 0)
                    converged_count += evaluate_cut_component(Flist, snpfrag, clist, k, slist, HAP1);
                else converged_count++;
            }

            if (converged_count == components) {
                //fprintf(stdout, "Haplotype assembly terminated early because no improvement seen in blocks after %d iterations\n", CONVERGE);
                break;
            }
        }

        // H-TRANS ESTIMATION FOR HIC
        if (MAX_HIC_EM_ITER > 1){

            // Possibly break if we're done improving
            HIC_LL_SCORE = 0;
            for (k = 0; k < components; k++){
                HIC_LL_SCORE += clist[k].bestSCORE;
            }
            if (HIC_LL_SCORE >= OLD_HIC_LL_SCORE){
                break;
            }
            OLD_HIC_LL_SCORE = HIC_LL_SCORE;

            likelihood_pruning(snps, Flist, snpfrag, HAP1, 0); // prune for only very high confidence SNPs
            // estimate the h-trans probabilities for the next round
            estimate_htrans_probs(Flist, fragments, HAP1, snpfrag);
        }
    }

    // BLOCK SPLITTING
    new_components = components;
    if (SPLIT_BLOCKS){
        split_count = 0;
        for (k=0; k<components; k++){
            // attempt to split block
            split_count += split_block(HAP1, clist, k, Flist, snpfrag, &new_components);
        }
        if (split_count > 0){
            // regenerate clist if necessary
            free(clist);
            clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*new_components);
            generate_clist_structure(Flist, fragments, snpfrag, snps, new_components, clist);
        }
        components = new_components;
    }else if(ERROR_ANALYSIS_MODE && !HIC){
        for (k=0; k<components; k++){
            // run split_block but don't actually split, just get posterior probabilities
            split_block(HAP1, clist, k, Flist, snpfrag, &new_components);
        }
    }

    // PRUNE SNPS
    if (!SKIP_PRUNE){
        discrete_pruning(snps, fragments, Flist, snpfrag, HAP1);
        likelihood_pruning(snps, Flist, snpfrag, HAP1, CALL_HOMOZYGOUS);
    }
    // PRINT OUTPUT FILE
    fprintf_time(stderr, "OUTPUTTING PRUNED HAPLOTYPE ASSEMBLY TO FILE %s\n", outputfile);
    print_hapfile(clist, components, HAP1, Flist, fragments, snpfrag, variantfile, miscalls, outputfile);

    // FREE UP MEMORY
    for (i = 0; i < snps; i++) free(snpfrag[i].elist);
    for (i = 0; i < snps; i++) free(snpfrag[i].telist);
    component = 0;
    for (i = 0; i < snps; i++) {
        free(snpfrag[i].flist);
        free(snpfrag[i].alist);
        free(snpfrag[i].jlist);
        free(snpfrag[i].klist);

        if (snpfrag[i].component == i && snpfrag[i].csize > 1) // root node of component
        {
            free(clist[component].slist);
            component++;
        }
    }

    for (i = 0; i < components; i++) free(clist[i].flist);
    free(snpfrag);
    free(clist);
    free(Flist);

    return 0;
}

void check_input_0_or_1(char* x){
    if (!(strcmp(x, "0") == 0 || strcmp(x, "1") == 0)){
        fprintf(stderr, "\nERROR: Invalid input \"%s\" for <0/1> option flag.\n",x);
        exit(1);
    }
}

int main(int argc, char** argv) {
    // input arguments are initial fragment file, variant file with variant information and alleles for each variant
    // number of iterations total, when to output the solution, file to output solution .....
    int i = 0;
    int flag = 0;
    char fragfile[10000];
    char varfile[10000];
    char VCFfile[10000];
    char hapfile[10000];
    strcpy(fragfile, "None");
    strcpy(varfile, "None");
    strcpy(hapfile, "None");
    strcpy(HTRANS_DATA_INFILE, "None");
    strcpy(HTRANS_DATA_OUTFILE, "None");

    if (argc % 2 != 1){
        fprintf(stderr, "\nERROR: Invalid number of arguments specified.\n");
        exit(1);
    }

    for (i = 1; i < argc; i += 2) {
        if (argc < 6) break;

        // BASIC OPTIONS
        if (strcmp(argv[i], "--fragments") == 0 || strcmp(argv[i], "--f") == 0) {
            strcpy(fragfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
            strcpy(VCFfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "--out") == 0|| strcmp(argv[i], "--o") == 0) {
            strcpy(hapfile, argv[i + 1]);
            flag++;
        }else if ((strcmp(argv[i], "--converge") == 0) || (strcmp(argv[i], "--c") == 0)) {
            CONVERGE = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "--v") == 0){
            check_input_0_or_1(argv[i + 1]);
            VERBOSE = atoi(argv[i + 1]);
        }
        // READ-TECHNOLOGY OPTIONS
        else if (strcmp(argv[i], "--HiC") == 0 || strcmp(argv[i], "--hic") == 0){
            check_input_0_or_1(argv[i + 1]);
            if (atoi(argv[i + 1])){
                MAX_HIC_EM_ITER = 100; //atoi(argv[i + 1]);
                NEW_FRAGFILE_FORMAT = 1;
                HIC = 1;
            }
        }else if (strcmp(argv[i], "--long_reads") == 0 || strcmp(argv[i], "--lr") == 0){
            check_input_0_or_1(argv[i + 1]);
            LONG_READS = atoi(argv[i + 1]);
            AUTODETECT_LONGREADS = 0;
        }else if (strcmp(argv[i], "--QV_offset") == 0 || strcmp(argv[i], "--qv_offset") == 0 || strcmp(argv[i], "--qo") == 0){
            QVoffset = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--hic_htrans_file") == 0 || strcmp(argv[i], "--hf") == 0){
            NEW_FRAGFILE_FORMAT = 1;
            strcpy(HTRANS_DATA_INFILE, argv[i + 1]);
            HIC = 1;
        }
        // HAPLOTYPE POST-PROCESSING OPTIONS
        else if (strcmp(argv[i], "--threshold") == 0 || strcmp(argv[i], "--t") == 0){
            THRESHOLD = 1.0 - unphred(atof(argv[i + 1]));
        //}else if (strcmp(argv[i], "--split_blocks") == 0 || strcmp(argv[i], "--sb") == 0){
            //SPLIT_BLOCKS = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--split_threshold") == 0 || strcmp(argv[i], "--st") == 0){
            SPLIT_THRESHOLD = 1.0 - unphred(atof(argv[i + 1]));
        }else if (strcmp(argv[i], "--call_homozygous") == 0 || strcmp(argv[i], "--ch") == 0){
            check_input_0_or_1(argv[i + 1]);
            CALL_HOMOZYGOUS = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--discrete_pruning") == 0 || strcmp(argv[i], "--dp") == 0){
            check_input_0_or_1(argv[i + 1]);
            DISCRETE_PRUNING = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--error_analysis_mode") == 0 || strcmp(argv[i], "--ea") == 0){
            check_input_0_or_1(argv[i + 1]);
            ERROR_ANALYSIS_MODE = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--SNVs_before_indels") == 0 || strcmp(argv[i], "--si") == 0){
            check_input_0_or_1(argv[i + 1]);
            SNVS_BEFORE_INDELS = atoi(argv[i + 1]);
        }
        // ADVANCED OPTIONS
        else if (strcmp(argv[i], "--nf") == 0 || strcmp(argv[i], "--new_format") == 0){
            check_input_0_or_1(argv[i + 1]);
            NEW_FRAGFILE_FORMAT = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--max_iter") == 0 || strcmp(argv[i], "--mi") == 0){
            MAXITER = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--maxcut_iter") == 0 || strcmp(argv[i], "--mc") == 0) {
            MAXCUT_ITER = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--htrans_read_lowbound") == 0 || strcmp(argv[i], "--hrl") == 0){
            HTRANS_READ_LOWBOUND = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--htrans_max_window") == 0 || strcmp(argv[i], "--hmw") == 0){
            HTRANS_MAX_WINDOW = atoi(argv[i + 1]);
        }
        // HIDDEN OPTIONS
        else if (strcmp(argv[i], "--htrans_data_outfile") == 0){
            strcpy(HTRANS_DATA_OUTFILE, argv[i + 1]);
        }else if (strcmp(argv[i], "--printscores") == 0 || strcmp(argv[i], "--scores") == 0){
            PRINT_FRAGMENT_SCORES = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--mbq") == 0){
            MINQ = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--skip_prune") == 0 || strcmp(argv[i], "--sp") == 0){
            check_input_0_or_1(argv[i + 1]);
            SKIP_PRUNE = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--max_IS") == 0 || strcmp(argv[i], "--mi") == 0){
            MAX_IS = atoi(argv[i + 1]);
        }else{
            fprintf(stderr, "\nERROR: Invalid Option \"%s\" specified.\n",argv[i]);
            exit(1);
        }
    }

    if (ERROR_ANALYSIS_MODE && HIC){
        fprintf_time(stderr,"WARNING: Switch error quality scores are not intended for use with Hi-C data. Scores will be left blank.\n");
    }

    if (flag != 3) // three essential arguments are not supplied
    {
        print_hapcut_options();
        return 0;
    }

	fprintf(stderr, "\n\n");
    fprintf_time(stderr, "fragment file: %s\n", fragfile);
    fprintf_time(stderr, "variantfile (VCF format):%s\n", VCFfile);
    fprintf_time(stderr, "haplotypes will be output to file: %s\n", hapfile);
    fprintf_time(stderr, "solution convergence cutoff: %d\n", CONVERGE);
    fprintf_time(stderr, "QVoffset: %d\n\n", QVoffset);
    maxcut_haplotyping(fragfile, VCFfile, hapfile);
    return 0;
}
