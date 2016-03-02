// CODE STARTED SEPT 10 2007 4pm //  april 8 2008 this code used for producing results in ECCB 2008 paper //
// author: VIKAS BANSAL (vbansal@scripps.edu) last modified December 23, 2010 

// time ./HAPCUT --fragments /home/vbansal-scripps/Haplotype-assembly/Dec8hairs/NA18507.allfour.matrix.SORTED --variants NA18507.chr6.hetvariants.inputforhapcut --output NA18507.phased --maxiter 60 > na18507.out 

/******** THINGS TO DO ***************************
0. format of fragment file changed: last string for each fragment now has the quality string (offset 33) for whole fragment // changed Feb 20 2011
1. to speed up computation, ignore components for which the MEC score has been stable for 5-10 iterations // TODO
2. implement weighted hapcut version where we take the quality values into account  // implemented this Feb 20 2011 
3. general model for errors: false variants, sequencing errors (bases) and chimeric paired-end reads.....
4. output MEC score with each component in solution  // DONE 
5. output to VCF format and maybe read in data from VCF file .....partly done 
6. how to determine bad SNPs/variants that cause excessive MEC score // working on this 
7. how to assign confidence to haplotypes... individual or whole 
 **********************************************/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "common.h"
#include "fragmatrix.h"
#include "pointerheap.h"
//#include "annealing.h"
#include "readinputfiles.h"
#include "removefalsehets.c"
#include "printhaplotypes.c"
//#define MAXBUF 10000


int RANDOM_START = 1;
int MINCUTALGO = 1;
int QVoffset = 33;
int VCFformat = 0;
// for testing, minq is 0
int MINQ = 0; // additional base quality filter in hapcut added april 18 2012
int MAXCUT_ITER = 100; // maximum number of iterations for max-cut algorithm, if this is proportional to 'N' -> complexity is 'N^2', added march 13 2013
int FOSMIDS = 1; // if this variable is 1, the data is fosmid long read data 
int HIC = 0;
int SCORING_FUNCTION = 5; // 0 = MEC score, 1 = switches
float THRESHOLD = 0.8;
float SPLIT_THRESHOLD = 0.8;
float HOMOZYGOUS_PRIOR = -80; // in log form. assumed to be really unlikely
int PRINT_FRAGMENT_SCORES = 0; // output the MEC/switch error score of erroneous reads/fragments to a file for evaluation 
int MAX_MEMORY = 8000;
int CONVERGE = 1; // stop iterations on a given component/block if exceed this many iterations since improvement
int NEW_CODE = 1; // likelihood based, max-cut calculated using partial likelihoods
int VERBOSE = 0;
int SPLIT_BLOCKS = 0;
int SPLIT_BLOCKS_MAXCUT = 0;
int REFHAP_HEURISTIC = 0;
int ERROR_ANALYSIS_MODE = 0;
int* iters_since_improvement;
int* iters_since_split;

#include "find_maxcut.c"   // function compute_good_cut 

int maxcut_haplotyping(char* fragmentfile, char* variantfile, int snps, char* outputfile, char* htrans_file, int maxiter_hapcut) {
    // IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+
    fprintf(stderr, "calling MAXCUT based haplotype assembly algorithm\n");
    int fragments = 0, iter = 0, components = 0;
    int i = 0, j = 0, k = 0, t = 0, component;
    int* slist;
    int flag = 0;
    float bestscore_mec = 0, calls = 0, miscalls = 0, ll = 0;
    char buffer[MAXBUF];

    /****************************** READ FRAGMENT MATRIX*************************************************/
    struct fragment* Flist;
    FILE* ff = fopen(fragmentfile, "r");
    if (ff == NULL) {
        fprintf(stderr, "couldn't open fragment file %s\n", fragmentfile);
        exit(0);
    }
    fragments = 0;
    while (fgets(buffer, MAXBUF, ff) != NULL) fragments++;
    fclose(ff);
    Flist = (struct fragment*) malloc(sizeof (struct fragment)*fragments);
    flag = read_fragment_matrix(fragmentfile, Flist, fragments);
    if (flag < 0) {
        fprintf(stderr, "unable to read fragment matrix file %s \n", fragmentfile);
        return -1;
    }

    if (VCFformat == 0) snps = count_variants(variantfile);
    else snps = count_variants_vcf(variantfile);
    if (snps < 0) {
        fprintf(stderr, "unable to read variant file %s \n", variantfile);
        return -1;
    }
    fprintf(stderr, "processed fragment file and variant file: fragments %d variants %d\n", fragments, snps);
    /****************************** READ FRAGMENT MATRIX*************************************************/

    struct SNPfrags* snpfrag = (struct SNPfrags*) malloc(sizeof (struct SNPfrags)*snps);
    update_snpfrags(Flist, fragments, snpfrag, snps, &components);
    double MEM_ALLOC = 0;
    for (i = 0; i < snps; i++) MEM_ALLOC += snpfrag[i].edges * 0.002;
    MEM_ALLOC *= 0.008;
    fprintf(stderr, "%f MB memory needs to be allocated for graph edges\n", MEM_ALLOC); // size of struct edge is 16/1000 bytes

    /*	
    //components =  determine_connected_components(Flist,fragments,snpfrag,snps);
    if (MEM_ALLOC >= MAX_MEMORY) 
    {
            fprintf(stderr,"\nstoring the HAPCUT graph structure requires more than %d MB of memory:\n 1. increase the maximum memory available using option \"--maxmem 12000\" where the memory is specified in megabytes OR \n 2. run the program with the options \"--longreads 1 \" to reduce the number of edges stored \n\n",MAX_MEMORY); 
            return -1;
    } 

    // too much memory allocated here for fosmid based data... 
    for (i=0;i<snps;i++) snpfrag[i].telist = (struct edge*)malloc(sizeof(struct edge)*snpfrag[i].edges);  

     */

    // 10/25/2014, edges are only added between adjacent nodes in each fragment and used for determining connected components... 
    // elist data structure is not used in hapcut algorithm anymore... 
    for (i = 0; i < snps; i++) snpfrag[i].elist = (struct edge*) malloc(sizeof (struct edge)*snpfrag[i].edges);
    if (FOSMIDS == 0) add_edges(Flist, fragments, snpfrag, snps, &components);
    else if (FOSMIDS >= 1) add_edges_fosmids(Flist, fragments, snpfrag, snps, &components);
    //add_edges_fosmids(Flist,fragments,snpfrag,snps,&components);
    for (i = 0; i < snps; i++) snpfrag[i].telist = (struct edge*) malloc(sizeof (struct edge)*snpfrag[i].edges);
    //add_edges(Flist,fragments,snpfrag,snps,&components);

    // this considers only components with at least two nodes
    fprintf(stderr, "fragments %d snps %d component(blocks) %d\n", fragments, snps, components);

    struct BLOCK* clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*components);
    component = 0;
    generate_clist_structure(Flist, fragments, snpfrag, snps, components, clist);
    
    /*****************************************************************************************************/
    // read in file with estimated probabilities of Hi-C h-trans interactions with distance
    int num_bins        = count_htrans_bins(htrans_file);
    int* htrans_bins    = (int*) malloc(sizeof(int) * num_bins);
    float* htrans_probs = (float*) malloc(sizeof(float) * num_bins);
    read_htrans_file(htrans_file, htrans_bins, htrans_probs, num_bins);
    /*****************************************************************************************************/

    char* HAP1 = (char*) malloc(snps + 1);
    char* besthap_mec = (char*) malloc(snps + 1);
    char* HAP2 = (char*) malloc(snps + 1);
    struct tm *ts1;
    char buf[80];
    time_t now;
    int split_occured;
    slist = (int*) malloc(sizeof (int)*snps);
    char fn[1000];

    if (VCFformat == 0) read_variantfile(variantfile, snpfrag, snps);
    else read_vcffile(variantfile, snpfrag, snps);
    
    /*****************************************************************************************************/
    int count=0;
    if (RANDOM_START == 1) {
        fprintf(stdout, "starting from a completely random solution SOLUTION \n\n");
        for (i = 0; i < snps; i++) {
            if (snpfrag[i].frags == 0) {
                HAP1[i] = '-';
                HAP2[i] = '-';
            } else {
                count++;
                if (drand48() < 0.5) {
                    HAP1[i] = '0';
                    HAP2[i] = '1';
                } else {
                    HAP1[i] = '1';
                    HAP2[i] = '0';
                }
            }
        }
    }
    
    float avg_depth, total = 0;
    for (i=0; i < snps; i++){
        total += snpfrag[i].frags;
    }
    avg_depth = total / ((float) snps);
    fprintf(stdout, "Avg fragment depth: %f\n",avg_depth);
    
    //frag_cluster_initialize(Flist, fragments, snpfrag, HAP1, snps, clist, components);
    //improve_hap(HAP1,clist,components, snps, fragments, Flist, snpfrag);
    
    for (i = 0; i < snps; i++) {
        besthap_mec[i] = HAP1[i];
    }

    // for each block, we maintain best haplotype solution under MFR criterion 
    // compute the component-wise score for 'initHAP' haplotype 
    miscalls = 0;
    bestscore_mec = 0;
    for (k = 0; k < components; k++) {
        clist[k].split = 0;
        clist[k].MEC = 0;
        clist[k].bestMEC = 0;
        clist[k].calls = 0;
        clist[k].LL = 0;
        for (i = 0; i < clist[k].frags; i++) {
            update_fragscore(Flist, clist[k].flist[i], HAP1);
            clist[k].MEC += Flist[clist[k].flist[i]].currscore;
            clist[k].LL += Flist[clist[k].flist[i]].ll;
            clist[k].calls += Flist[clist[k].flist[i]].calls;

        }
        clist[k].bestMEC = clist[k].MEC;
        bestscore_mec += clist[k].bestMEC;
        miscalls += clist[k].MEC;
        clist[k].bestLL = clist[k].LL;
    }

    //	annealing_haplotyping(Flist,fragments,snpfrag,snps,maxiter,HAP1,HAP2,clist,components,slist); return 1;
    //	annealing_haplotyping_full(Flist,fragments,snpfrag,snps,maxiter,HAP1,HAP2,0); return 1;

    int trueMEC = 0, converged_count=0, split_count, new_components;
    
    // counter arrays where iters_since_whatever[i] refers to component at offset i
    // maintained OUTSIDE of the clist structure since it gets regenerated frequently
    iters_since_improvement = (int*) calloc(snps,sizeof(int));
    iters_since_split = (int*) calloc(snps,sizeof(int));

    /************************** RUN THE MAX_CUT ALGORITHM ITERATIVELY TO IMPROVE MEC SCORE*********************************/
    for (iter = 0; iter < maxiter_hapcut; iter++) {
        trueMEC = mecscore(Flist, fragments, HAP1, &ll, &calls, &miscalls);
        if (SCORING_FUNCTION == 5) miscalls = trueMEC;
        time(&now);
        ts1 = localtime(&now);
        strftime(buf, sizeof (buf), "%a %Y-%m-%d %H:%M:%S %Z", ts1);
        if(VERBOSE) fprintf(stdout, "iter %d current haplotype MEC %f calls %d LL %f %s \n", iter, miscalls, (int) calls, ll, buf);
        //fprintf(stderr, "iter %d current haplotype MEC %f calls %d LL %f %s \n", iter, miscalls, (int) calls, ll, buf);
        if ((iter % 5 == 0 && iter > 0)) {
            // new code added april 7 2012
            for (k = 0; k < components; k++) find_bestvariant_segment(Flist, fragments, snpfrag, clist, k, HAP1, HAP2);

            sprintf(fn, "%s", outputfile); // newfile for every update to score....
            //sprintf(fn,"%s.%f",outputfile,miscalls);   // newfile for every update to score....
            fprintf(stdout, "OUTPUTTING HAPLOTYPE ASSEMBLY TO FILE %s\n", fn);
            //fprintf(stderr, "OUTPUTTING HAPLOTYPE ASSEMBLY TO FILE %s\n", fn);
            //if (VCFformat ==1) print_haplotypes_vcf(clist,components,HAP1,Flist,fragments,snpfrag,snps,fn);
            print_hapfile(clist, components, HAP1, Flist, fragments, snpfrag, variantfile, miscalls, fn);

            // do this only if some option is specified 
            if (PRINT_FRAGMENT_SCORES == 1) print_fragmentmatrix_MEC(Flist, fragments, HAP1, outputfile);

        }
        converged_count = 0; split_occured = 0;
        new_components = components;
        for (k = 0; k < components; k++) // COMPUTATION OF TREE FOR EACH COMPONENT 
        {
            if(VERBOSE && iter == 0) fprintf(stdout, "component %d length %d phased %d %d...%d \n", k, clist[k].length, clist[k].phased, clist[k].offset, clist[k].lastvar);
            // call function for each component only if MEC > 0 april 17 2012
            if (clist[k].MEC > 0) converged_count += evaluate_cut_component(Flist, snpfrag, clist, k, slist, HAP1, iter, &new_components);
            else converged_count++;
            if (clist[k].split) split_occured = 1;
        }
        if (split_occured){
            // regenerate clist because there are newly split blocks
            free(clist);
            clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*new_components);
            generate_clist_structure(Flist, fragments, snpfrag, snps, new_components, clist);
            components = new_components;
        }
        if (converged_count == components) {
            fprintf(stdout, "Haplotype assembly terminated early because no improvement seen in blocks after %d iterations\n", CONVERGE);
            break;
        }
    }
    
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
    }else if(ERROR_ANALYSIS_MODE){
        for (k=0; k<components; k++){
            // run split_block but don't actually split, just get posterior probabilities
            split_block(HAP1, clist, k, Flist, snpfrag, &new_components);
        }
    }
    components = new_components;

    for (k = 0; k < components; k++) find_bestvariant_segment(Flist, fragments, snpfrag, clist, k, HAP1, HAP2);

    refhap_heuristic(snps, fragments, Flist, snpfrag, HAP1);
    prune_snps(snps, Flist, snpfrag, HAP1);
    
    fprintf(stdout, "OUTPUTTING PRUNED HAPLOTYPE ASSEMBLY TO FILE %s\n", outputfile);
    //if (VCFformat ==1) print_haplotypes_vcf(clist,components,HAP1,Flist,fragments,snpfrag,snps,fn);
    print_hapfile(clist, components, HAP1, Flist, fragments, snpfrag, variantfile, miscalls, outputfile);

    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    // input arguments are initial fragment file, variant file with variant information and alleles for each variant 
    // number of iterations total, when to output the solution, file to output solution .....
    time_t ts;
    time(&ts);
    srand48((long int) ts);
    if (MINCUTALGO == 2) RANDOM_START = 0;
    int i = 0, j = 0;
    int flag = 0;
    char fragfile[10000];
    char varfile[10000];
    char VCFfile[10000];
    char hapfile[10000];
    char htrans_file[10000];
    int maxiter = 100;
    strcpy(fragfile, "None");
    strcpy(varfile, "None");
    strcpy(hapfile, "None");
    strcpy(htrans_file, "None");
    for (i = 1; i < argc; i += 2) {
        if (argc < 6) break;
        if (strcmp(argv[i], "--fragments") == 0 || strcmp(argv[i], "--frags") == 0) {
            strcpy(fragfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--variants") == 0) {
            j = strlen(argv[i + 1]); // check if it is a VCF file ending with .vcf or .VCF 
            if (j > 3 && ((argv[i + 1][j - 1] == 'f' && argv[i + 1][j - 2] == 'c' && argv[i + 1][j - 3] == 'v') || (argv[i + 1][j - 1] == 'F' && argv[i + 1][j - 2] == 'C' && argv[i + 1][j - 3] == 'V'))) {
                strcpy(VCFfile, argv[i + 1]);
                VCFformat = 1;
                flag++;
            } else {
                fprintf(stderr, "please provide variant file in VCF format using --VCF option, old variant format is no longer supported\n\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
            strcpy(VCFfile, argv[i + 1]);
            VCFformat = 1;
            flag++;
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "--out") == 0) {
            strcpy(hapfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--threshold") == 0 || strcmp(argv[i], "--t") == 0) THRESHOLD = atof(argv[i + 1]);
        else if (strcmp(argv[i], "--split_threshold") == 0 || strcmp(argv[i], "--st") == 0) SPLIT_THRESHOLD = atof(argv[i + 1]);
        else if (strcmp(argv[i], "--maxiter") == 0) maxiter = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--longreads") == 0 || strcmp(argv[i], "--lr") == 0) // long reads pacbio 
        {
            FOSMIDS = atoi(argv[i + 1]);
        } else if ((strcmp(argv[i], "--converge") == 0) || (strcmp(argv[i], "--c") == 0)) {
            CONVERGE = atoi(argv[i + 1]);
        } else if (strcmp(argv[i], "--MEC") == 0) {
            NEW_CODE = !(atoi(argv[i + 1]));
            if (!NEW_CODE){
                SCORING_FUNCTION = 0; // for new code, scoring function is likelihood based and mec of fragment = -1*LL
                //fprintf(stderr, "SCORING FUNCTION var set to 0\n", NEW_CODE);
            }
        }else if (strcmp(argv[i], "--fosmid") == 0 || strcmp(argv[i], "--fosmids") == 0) FOSMIDS = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--HiC_htrans_file") == 0 || strcmp(argv[i], "--htrans") == 0){
            strcpy(htrans_file, argv[i + 1]);
            HIC = 1;
        }
        else if (strcmp(argv[i], "--sf") == 0 || strcmp(argv[i], "--switches") == 0) SCORING_FUNCTION = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--printscores") == 0 || strcmp(argv[i], "--scores") == 0) PRINT_FRAGMENT_SCORES = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--maxcutiter") == 0) {
            MAXCUT_ITER = atoi(argv[i + 1]);
            fprintf(stderr, "max iterations for max-cut calculations is %d \n", MAXCUT_ITER);
        }else if (strcmp(argv[i], "--QVoffset") == 0 || strcmp(argv[i], "--qvoffset") == 0) QVoffset = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--maxmem") == 0) MAX_MEMORY = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--mbq") == 0) MINQ = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--splitblocks") == 0) SPLIT_BLOCKS = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--splitblocks_maxcut") == 0) SPLIT_BLOCKS_MAXCUT = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "--v") == 0) VERBOSE = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--refhap_heuristic") == 0 || strcmp(argv[i], "--rh") == 0) REFHAP_HEURISTIC = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--error_analysis_mode") == 0 || strcmp(argv[i], "--ea") == 0) ERROR_ANALYSIS_MODE = atoi(argv[i + 1]);
    }
    if (flag != 3) // three essential arguments are not supplied 
    {
        print_hapcut_options();
        return 0;
    } else {
        
        if (NEW_CODE){
            fprintf(stdout, "USING LOG-LIKELIHOOD BASED HAPCUT\n");
        }else{
            fprintf(stdout, "USING MEC BASED HAPCUT\n");
        }
        
        if (VCFformat == 1) {
            fprintf(stderr, "\n\nfragment file: %s\nvariantfile (VCF format):%s\nhaplotypes will be output to file: %s\niterations of maxcut algorithm: %d\nQVoffset: %d\n\n", fragfile, VCFfile, hapfile, maxiter, QVoffset);
            maxcut_haplotyping(fragfile, VCFfile, 0, hapfile, htrans_file, maxiter);
        }else{
            fprintf(stderr, "\n\nfragment file: %s\nvariantfile (variant format):%s\nhaplotypes will be output to file: %s\niterations of maxcut algorithm: %d\nQVoffset: %d\n\n", fragfile, varfile, hapfile, maxiter, QVoffset);
            maxcut_haplotyping(fragfile, varfile, 0, hapfile, htrans_file, maxiter);
        }
    }
    return 1;
}




