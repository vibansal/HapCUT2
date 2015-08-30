// Simple Haplotyping tool forked from Hapcut code

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include "simple.h"
#include "readinputfiles.h"
#include "printhaplotypes.c"
#define MAXBUF 10000

int QVoffset = 33;
float HOMOZYGOUS_PRIOR = -80; // assumed to be really unlikely
int MAX_MEMORY = 8000;

/***********************************************************************************************************/
// this function updates the data structure snpfrag which links the heterozyous SNPs and the haplotype fragments
// // it generates a list of fragments (flist) that affect each SNP 

void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps) {
    int i = 0, j = 0, k = 0;

    // find the first fragment whose endpoint lies at snp 'i' or beyond
    for (i = 0; i < snps; i++) {
        snpfrag[i].frags = 0;
        snpfrag[i].ff = -1;
    }
    for (i = 0; i < fragments; i++) {
        j = Flist[i].list[0].offset;
        k = Flist[i].list[Flist[i].blocks - 1].len + Flist[i].list[Flist[i].blocks - 1].offset;
    }

    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) snpfrag[Flist[i].list[j].offset + k].frags++;
        }
    }
    for (i = 0; i < snps; i++) {
        snpfrag[i].flist = (int*) malloc(sizeof (int)*snpfrag[i].frags);
        snpfrag[i].alist = (char*) malloc(snpfrag[i].frags);
    }

    for (i = 0; i < snps; i++) {
        snpfrag[i].component = -1;
        snpfrag[i].csize = 1;
        snpfrag[i].frags = 0;
    }

    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                snpfrag[Flist[i].list[j].offset + k].flist[snpfrag[Flist[i].list[j].offset + k].frags] = i;
                snpfrag[Flist[i].list[j].offset + k].alist[snpfrag[Flist[i].list[j].offset + k].frags] = Flist[i].list[j].hap[k];
                snpfrag[Flist[i].list[j].offset + k].frags++;
            }
        }
    }
}

// new functions for connected components using disjoint union data structure feb 14 2013 

int find_parent(struct SNPfrags* snpfrag, int x) {
    // all nodes on path to root will be labeled with the same id
    if (snpfrag[x].component != x) snpfrag[x].component = find_parent(snpfrag, snpfrag[x].component);
    return snpfrag[x].component;
}

// disjoint union data structure where for each node we maintain information about its parent
// following parent to parent pointers -> information about connected component (i)

int determine_connected_components(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps) {
    int i = 0, j = 0, k = 0, x = 0, y = 0, parent_x = 0, parent_y = 0;
    for (i = 0; i < snps; i++) {
        snpfrag[i].component = i; // initialize every node to a singleton set
        snpfrag[i].rank = 0;
        snpfrag[i].csize = 1;
    }

    // each fragment is a set of nodes that should end up in the same connected component 
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                if (j == 0 && k == 0) continue;
                x = Flist[i].list[0].offset;
                y = Flist[i].list[j].offset + k; // the pair of SNPs
                parent_x = find_parent(snpfrag, x);
                parent_y = find_parent(snpfrag, y);
                if (parent_x == parent_y) continue;
                else if (snpfrag[parent_x].rank < snpfrag[parent_y].rank) {
                    snpfrag[parent_x].component = parent_y;
                    snpfrag[parent_y].csize += snpfrag[parent_x].csize;
                    snpfrag[parent_x].csize = 0;
                } else if (snpfrag[parent_x].rank > snpfrag[parent_y].rank) {
                    snpfrag[parent_y].component = parent_x;
                    snpfrag[parent_x].csize += snpfrag[parent_y].csize;
                    snpfrag[parent_y].csize = 0;
                } else {
                    snpfrag[parent_y].component = parent_x;
                    snpfrag[parent_x].rank++;
                    snpfrag[parent_x].csize += snpfrag[parent_y].csize;
                    snpfrag[parent_y].csize = 0;
                }
            }
        }
    }

    int components = 0, singletons = 0;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1)components++;
        else if (snpfrag[i].component == i){
            snpfrag[i].component = -1; // -1 signifies singleton
            singletons++;
        }
    }
    fprintf(stdout, " no of non-trivial connected components in graph is %d singletons %d mean size %0.1f\n", components, singletons, (float) (snps - singletons) / components);
    return components;
}
// compute score of fragment
// don't mutate Flist or other structures.
// return score as return value
// homozygous: 0-based index of a homozygous position. -1 if no homozygous pos
// switch_ix: 0-based index of the switch error being tested, -1 if none

float simple_fragscore(struct fragment* Flist, int f, char* h, int homozygous, int switch_ix) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    float good = 0, bad = 0, ll;
    int snp_ix, switched;

    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs

            // conditions to skip this base
            if (h[snp_ix] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}

            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
            prob2 = Flist[f].list[j].p1[k];

            if (h[snp_ix] == Flist[f].list[j].hap[k]) good += prob1;
            else bad += prob1;

            // this is likelihood based calculation
            switched = (switch_ix != -1 && snp_ix >= switch_ix);
            if (snp_ix != homozygous) { // not homozygous
                if ((h[snp_ix] == Flist[f].list[j].hap[k]) != switched) { // true if match, or not match but switched
                    p0 += prob2;
                    p1 += prob;
                } else {
                    p0 += prob;
                    p1 += prob2;
                }
            } else { // homozygous at this postion
                if ((h[snp_ix] == Flist[f].list[j].hap[k]) != switched) { // true if match, or not match but switched
                    p0 += prob2; // both hap1 and hap2 match
                    p1 += prob2;
                } else {
                    p0 += prob;
                    p1 += prob;
                }
            }
        }
    }
    if (p0 > p1) ll = (p0 + log10(1 + pow(10, p1 - p0)));
    else ll = (p1 + log10(1 + pow(10, p0 - p1)));

    return ll;
}

// given a and b, where
// a = log10(x)
// b = log10(y)
// returns log10(x+y)

float addlogs(float a, float b) {
    if (a > b)
        return (a + log10(1 + pow(10, b - a)));
    else
        return (b + log10(1 + pow(10, a - b)));
}

// given a and b, where
// a = log10(x)
// b = log10(y)
// returns log10(x-y)

float subtractlogs(float a, float b) {
    if (a > b)
        return (a + log10(1 - pow(10, b - a)));
    else
        return (b + log10(1 - pow(10, a - b)));
}

// returns an array length [snps] called 'pruned' that indicates which SNPs were pruned:
// 0 indicates not pruned
// 1 indicates pruned (leave unphased in output)
// 2 indicates called homozygous 00
// 3 indicates called homozygous 11
// the fraction to prune is determined by prune_threshold if THRESHOLD_TYPE is 1.
// else, the default behavior is that prune_threshold is the posterior probability cutoff.

void prune_snps(char* pruned, float prune_threshold, float homozygous_threshold, float switch_threshold, int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP) {
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, P_data_Hsw, total, temp;
    float log_hom_prior;
    float log_het_prior;
    struct hap_prob* hap_probs = malloc(snps * sizeof (struct hap_prob));

    for (i = 0; i < snps; i++) {

        // get prior probabilities for homozygous and heterozygous genotypes
        log_hom_prior = snpfrag[i].homozygous_prior;
        log_het_prior = subtractlogs(log10(0.5), log_hom_prior);

        // we only care about positions that are haplotyped
        if (!(HAP[i] == '1' || HAP[i] == '0')) {
            hap_probs[i].post_hap = 1;
            hap_probs[i].post_hapf = 0;
            hap_probs[i].post_00 = 0;
            hap_probs[i].post_11 = 0;
            hap_probs[i].post_sw = 0;
            continue;
        }

        // hold on to original haplotype values
        temp1 = HAP[i];
        // set probabilities to zero
        P_data_H = 0;
        P_data_Hf = 0;
        P_data_H00 = 0;
        P_data_H11 = 0;
        P_data_Hsw = 0; // switch error at i

        //looping over fragments overlapping i and sum up read probabilities
        for (j = 0; j < snpfrag[i].frags; j++) {

            f = snpfrag[i].flist[j];
            // normal haplotypes
            P_data_H += simple_fragscore(Flist, f, HAP, -1, -1);
            // haplotypes with switch error starting at i
            P_data_Hsw += simple_fragscore(Flist, f, HAP, -1, i);
            // haplotypes with i flipped
            if (HAP[i] == '1')
                HAP[i] = '0';
            else if (HAP[i] == '0')
                HAP[i] = '1';
            P_data_Hf += simple_fragscore(Flist, f, HAP, -1, -1);
            // haplotypes with i homozygous 00
            HAP[i] = '0';
            P_data_H00 += simple_fragscore(Flist, f, HAP, i, -1);
            // haplotypes with i homozygous 11
            HAP[i] = '1';
            P_data_H11 += simple_fragscore(Flist, f, HAP, i, -1);
            //return haplotype to original value
            HAP[i] = temp1;
        }

        // denominator of posterior probabilities;
        // sum of all 4 data probabilities times their priors
        total = addlogs(
                addlogs((log_het_prior + P_data_H), (log_het_prior + P_data_Hf)),
                addlogs((log_hom_prior + P_data_H00), (log_hom_prior + P_data_H11)));

        hap_probs[i].post_hap = log_het_prior + P_data_H - total;
        hap_probs[i].post_hapf = log_het_prior + P_data_Hf - total;
        hap_probs[i].post_00 = log_hom_prior + P_data_H00 - total;
        hap_probs[i].post_11 = log_hom_prior + P_data_H11 - total;
        hap_probs[i].post_sw = P_data_Hsw - addlogs(P_data_H, P_data_Hsw);
    }

    // figure out if block should be split at this position, or switched
    for (i = 0; i < snps; i++) {

        // we only care about positions that are haplotyped
        if (!(HAP[i] == '1' || HAP[i] == '0')) continue;

        //fprintf(stderr, "%f  %f\n", 1-pow(10,hap_probs[i].post_sw),pow(10,hap_probs[i].post_sw));
        if (hap_probs[i].post_hap < log10(prune_threshold)) {
            pruned[i] = 1;
            continue;
        }

        // change the status of SNPs that are above the homozygous threshold
        // doing this in a separate loop just to be more decoupled
        // flip the haplotype at this position if necessary
        if (hap_probs[i].post_hapf > hap_probs[i].post_hap) {
            temp = hap_probs[i].post_hap;
            hap_probs[i].post_hap = hap_probs[i].post_hapf;
            hap_probs[i].post_hapf = temp;
            if (HAP[i] == '1')
                HAP[i] = '0';
            else if (HAP[i] == '0')
                HAP[i] = '1';
        }

        if (hap_probs[i].post_00 > log10(homozygous_threshold)) {
            pruned[i] = 2; // 2 specifies 00 homozygous
            continue;
        } else if (hap_probs[i].post_11 > log10(homozygous_threshold)) { // 3 specifies 11 homozygous
            pruned[i] = 3;
            continue;
        }

        if (hap_probs[i].post_sw > log10(switch_threshold)) {
            // block should be switched here
            for (j = i; j < snps; j++) {

                // need to switch this position
                if (HAP[j] == '1')
                    HAP[j] = '0';
                else if (HAP[j] == '0')
                    HAP[j] = '1';

            }
            continue;
        }else if (subtractlogs(0, hap_probs[i].post_sw) < log10(switch_threshold)) {
            // probability that there wasn't a switch error is below threshold
            // block should be split here
            //pruned[i] = 4; // 4 specifies block split
        }

    }

    free(hap_probs);
}

float improve_hap(char* HAP, struct component* clist, int components, int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, int maxiter) {
    int iter;
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_Hsw, temp, new_ll, old_ll;
    float* post_hap = malloc(snps * sizeof (float));
    float* post_hapf = malloc(snps * sizeof (float));
    float* post_sw = malloc(snps * sizeof (float));
    
    // get initial log likelihood
    new_ll = 0;
    for (i = 0; i < fragments; i++) {
        new_ll += simple_fragscore(Flist, i, HAP, -1, -1);
    }

    
    for (iter = 0; iter < maxiter; iter++){
        for (i = 0; i < snps; i++) {

            // skip ahead if positions aren't haplotyped
            if (!(HAP[i] == '1' || HAP[i] == '0')) {
                continue;
            }

            // hold on to original haplotype values
            temp1 = HAP[i];
            // set probabilities to zero
            P_data_H = 0;
            P_data_Hf = 0;
            P_data_Hsw = 0; // switch error at i

            //looping over fragments overlapping i and sum up read probabilities
            for (j = 0; j < snpfrag[i].frags; j++) {

                f = snpfrag[i].flist[j];
                // normal haplotypes
                P_data_H += simple_fragscore(Flist, f, HAP, -1, -1);
                // haplotypes with switch error starting at i
                P_data_Hsw += simple_fragscore(Flist, f, HAP, -1, i);
                // haplotypes with i flipped
                if (HAP[i] == '1')
                    HAP[i] = '0';
                else if (HAP[i] == '0')
                    HAP[i] = '1';
                P_data_Hf += simple_fragscore(Flist, f, HAP, -1, -1);
                //return haplotype to original value
                HAP[i] = temp1;
            }

            // denominator of posterior probabilities;
            // sum of all 4 data probabilities times their priors

            post_hap[i] = P_data_H - addlogs(P_data_H, P_data_Hf);
            post_hapf[i] = subtractlogs(0, post_hap[i]);
            post_sw[i] = P_data_Hsw - addlogs(P_data_H, P_data_Hsw);
        }

        for (i = 0; i < snps; i++) {
            // skip ahead if positions aren't haplotyped
            if (!(HAP[i] == '1' || HAP[i] == '0')) {
                continue;
            }
            // change the status of SNPs that are above the homozygous threshold
            // doing this in a separate loop just to be more decoupled

            // flip the haplotype at this position if necessary

            if (post_hapf[i] > post_hap[i]) {
                //temp = post_hap[i];
                //post_hap[i] = post_hapf[i];
                //post_hapf[i] = temp[i];
                if (HAP[i] == '1')
                    HAP[i] = '0';
                else if (HAP[i] == '0')
                    HAP[i] = '1';
            }

            if (post_sw[i] > log10(0.8)) {
                // block should be switched here
                //fprintf(stderr, "PHASED: %d SLIST LEN: %d\n", clist[snpfrag[i].component].phased,sizeof(clist[snpfrag[i].component].slist)/sizeof(clist[snpfrag[i].component].slist[0]));
                for (j = i; j < snps; j++) {

                    // need to switch this position
                    if (HAP[j] == '1')
                        HAP[j] = '0';
                    else if (HAP[j] == '0')
                        HAP[j] = '1';

                }
                continue;
            }
        }
        
        old_ll = new_ll;
        new_ll = 0;
        for (i = 0; i < fragments; i++) {
            new_ll += simple_fragscore(Flist, i, HAP, -1, -1);
        }
        if (new_ll <= old_ll) break;
    }
    
    free(post_hap);
    free(post_hapf);
    free(post_sw);
    return iter;
}

int simple_haplotyping(char* fragmentfile, char* variantfile, int snps, char* outputfile, int maxiter, float prune_threshold, float homozygous_threshold, float switch_threshold) {
    // IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+

    int fragments, components, i, j, k, flag = 0, c, f,last_covered;

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
    Flist = (struct fragment*) calloc(fragments, sizeof (struct fragment));
    flag = read_fragment_matrix(fragmentfile, Flist, fragments); // include HAP so we can set covered positions while reading in matrix
    if (flag < 0) {
        fprintf(stderr, "unable to read fragment matrix file %s \n", fragmentfile);
        return -1;
    }

    snps = count_variants_vcf(variantfile);
    if (snps < 0) {
        fprintf(stderr, "unable to read variant file %s \n", variantfile);
        return -1;
    }
    fprintf(stderr, "processed fragment file and variant file: fragments %d variants %d\n", fragments, snps);

    /****************************** READ FRAGMENT MATRIX AND VCF*************************************************/

    struct SNPfrags* snpfrag = (struct SNPfrags*) malloc(sizeof (struct SNPfrags)*snps);
    update_snpfrags(Flist, fragments, snpfrag, snps);
    read_vcffile(variantfile, snpfrag, snps);
    
    components = determine_connected_components(Flist, fragments, snpfrag, snps);
    char* done = calloc(snps, sizeof(char));
    struct component* clist = (struct component*) malloc(components*sizeof(struct component));
    c = 0; 
    for (i = 0; i < snps; i++){
        if(!done[i] && snpfrag[i].component != -1){
            clist[c].slist = malloc(snpfrag[i].csize*sizeof(int));
            clist[c].size = snpfrag[i].csize;
            f = 0;
            for (j = 0; j < snps; j++){
                if (snpfrag[j].component == snpfrag[i].component){
                    clist[c].slist[f] = j;
                    done[j] = 1;
                    f++;
                }
            }
            c++;
        }
    }
    free(done);
    /*****************************************************************************************************/
    // initialize haplotype
    char* HAP = (char*) malloc(snps + 1);
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == -1){
            HAP[i] = '-';
        }else{
            if (drand48() < 0.5)
                HAP[i] = '0';
            else
                HAP[i] = '1';
        }
    }
    HAP[snps] = '\0';
    
    last_covered = -1;
    // "layer" on fragments to minimize the number of switch errors that need to be corrected.
    for (i = 0; i < fragments; i++) {
        if (Flist[i].list[0].offset < last_covered) continue;
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                if (snpfrag[Flist[i].list[j].offset + k].component == -1) continue;
                HAP[Flist[i].list[j].offset + k] = Flist[i].list[j].hap[k];
                last_covered = Flist[i].list[j].offset + k;
            }
        }
    }
    // get initial ll score

    // now attempt to improve it
    int iters = improve_hap(HAP, clist, components, snps, fragments, Flist, snpfrag, maxiter);
    fprintf(stdout, "Solution converged after %d iterations\n", iters);

    fprintf(stdout, "OUTPUTTING HAPLOTYPE ASSEMBLY TO FILE %s\n", outputfile);
    char* pruned = calloc(snps, sizeof (char)); // bit array that indicates pruned SNPs, initialize to 0 with calloc
    print_hapfile(clist, components, snps, HAP, Flist, fragments, snpfrag, variantfile, outputfile, pruned);
    prune_snps(pruned, prune_threshold, homozygous_threshold, switch_threshold, snps, Flist, snpfrag, HAP);
    char outputfile_pruned[1000];
    strcpy(outputfile_pruned, outputfile);
    strcat(outputfile_pruned, ".pruned");
    fprintf(stdout, "OUTPUTTING PRUNED HAPLOTYPE ASSEMBLY TO FILE %s\n", outputfile_pruned);
    print_hapfile(clist, components, snps,  HAP, Flist, fragments, snpfrag, variantfile, outputfile_pruned, pruned);

    free(pruned);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    // input arguments are initial fragment file, variant file with variant information and alleles for each variant 
    // number of iterations total, when to output the solution, file to output solution .....
    time_t ts;
    time(&ts);
    srand48((long int) ts);
    int i = 0;
    int flag = 0;
    float prune_threshold = 0.8, homozygous_threshold = 0.8, switch_threshold = 0.8;
    char fragfile[10000];
    char VCFfile[10000];
    char hapfile[10000];
    int maxiter = 10000;
    strcpy(fragfile, "None");
    strcpy(hapfile, "None");
    for (i = 1; i < argc; i += 2) {
        if (argc < 6) break;
        if (strcmp(argv[i], "--fragments") == 0 || strcmp(argv[i], "--frags") == 0) {
            strcpy(fragfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
            strcpy(VCFfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "--out") == 0) {
            strcpy(hapfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--prune_threshold") == 0 || strcmp(argv[i], "--pt") == 0) {
            prune_threshold = atof(argv[i + 1]);
        } else if (strcmp(argv[i], "--homozygous_threshold") == 0 || strcmp(argv[i], "--ht") == 0) {
            homozygous_threshold = atof(argv[i + 1]);
        } else if (strcmp(argv[i], "--switch_threshold") == 0 || strcmp(argv[i], "--st") == 0) {
            switch_threshold = atof(argv[i + 1]);
        } else if (strcmp(argv[i], "--maxiter") == 0) maxiter = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--QVoffset") == 0 || strcmp(argv[i], "--qvoffset") == 0) QVoffset = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--maxmem") == 0) MAX_MEMORY = atoi(argv[i + 1]);

    }
    if (flag != 3) // three essential arguments are not supplied 
    {
        print_hapcut_options();
        return 0;
    } else {
        fprintf(stderr, "\n\nfragment file: %s\nvariantfile (VCF format):%s\nhaplotypes will be output to file: %s\niterations of maxcut algorithm: %d\nQVoffset: %d\n\n", fragfile, VCFfile, hapfile, maxiter, QVoffset);
        simple_haplotyping(fragfile, VCFfile, 0, hapfile, maxiter, prune_threshold, homozygous_threshold, switch_threshold);
    }
    return 0;
}




