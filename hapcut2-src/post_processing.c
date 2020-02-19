#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

extern float THRESHOLD;
extern int DISCRETE_PRUNING;
extern int ERROR_ANALYSIS_MODE;
extern float HOMOZYGOUS_PRIOR;
extern float SPLIT_THRESHOLD;
extern int SPLIT_BLOCKS;
extern int SNVS_BEFORE_INDELS;

void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous);
int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr);

void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous) {
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, total, log_hom_prior, log_het_prior;
    float post_hap, post_hapf, post_00, post_11,ll;
    float PGLL[5]; 

    for (i = 0; i < snps; i++) {

        // "unmask" indel variant
        if (SNVS_BEFORE_INDELS && snpfrag[i].is_indel == '1') HAP1[i] = '0';

        // we only care about positions that are haplotyped
        if (!(HAP1[i] == '1' || HAP1[i] == '0')) continue;

        // hold on to original haplotype values
        temp1 = HAP1[i];
        // set probabilities to zero
        P_data_H = 0;
        P_data_Hf = 0;
        P_data_H00 = 0;
        P_data_H11 = 0;
	for (j=0;j<5;j++) PGLL[j] = 0; 

        //looping over fragments overlapping i and sum up read probabilities
        for (j = 0; j < snpfrag[i].frags; j++) {

            f = snpfrag[i].flist[j];
            P_data_H += fragment_ll1(Flist,f, HAP1, -1, -1);
            flip(HAP1[i]);
            P_data_Hf += fragment_ll1(Flist,f, HAP1, -1, -1);

            if (call_homozygous){
                // haplotypes with i homozygous 00
                HAP1[i] = '0';
                P_data_H00 += fragment_ll1(Flist,f, HAP1, i, -1);
                // haplotypes with i homozygous 11
                HAP1[i] = '1';
                P_data_H11 += fragment_ll1(Flist,f, HAP1, i, -1);
            }
            HAP1[i] = '0'; PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1); // added 03/12/2018
            HAP1[i] = '1'; PGLL[3] += fragment_ll1(Flist,f, HAP1, i, -1); // homozygous genotypes 00,11
            //return haplotype to original value
            HAP1[i] = temp1;
	    ll = fragment_ll1(Flist,f, HAP1, i, 10); // variant 'i' ignored for likelihood calculation if last parameter = 10
	    if (ll < 0) PGLL[4] +=ll;
        }
	PGLL[1] = P_data_H; PGLL[2] = P_data_Hf;


        if (call_homozygous){

            // get prior probabilities for homozygous and heterozygous genotypes
            log_hom_prior = snpfrag[i].homozygous_prior;
            log_het_prior = subtractlogs(log10(0.5), log_hom_prior);

            // denominator of posterior probabilities;
            // sum of all 4 data probabilities times their priors
            total = addlogsEXACT(
                    addlogsEXACT((log_het_prior + P_data_H), (log_het_prior + P_data_Hf)),
                    addlogsEXACT((log_hom_prior + P_data_H00), (log_hom_prior + P_data_H11)));

            post_hap = log_het_prior + P_data_H - total;
            post_hapf = log_het_prior + P_data_Hf - total;
            post_00 = log_hom_prior + P_data_H00 - total;
            post_11 = log_hom_prior + P_data_H11 - total;

            // change the status of SNPs that are above/below threshold
            if (post_00 > log10(0.5)){
                snpfrag[i].genotypes[0] = '0';
                snpfrag[i].genotypes[2] = '0';
                snpfrag[i].post_hap     = post_00;
            }else if (post_11 > log10(0.5)){
                snpfrag[i].genotypes[0] = '1';
                snpfrag[i].genotypes[2] = '1';
                snpfrag[i].post_hap     = post_11;
            }else if (post_hapf > log10(0.5)){
                flip(HAP1[i]);                // SNP should be flipped
                snpfrag[i].post_hap = post_hapf;
            }else{
                snpfrag[i].post_hap = post_hap;
            }
        } else {

            // get prior probabilities for homozygous and heterozygous genotypes
            log_het_prior = log10(0.5);

            total = addlogsEXACT((log_het_prior + P_data_H), (log_het_prior + P_data_Hf));

            post_hap = log_het_prior + P_data_H - total;
            post_hapf = log_het_prior + P_data_Hf - total;

            // change the status of SNPs that are above/below threshold
            if (post_hapf > log10(0.5)){
                flip(HAP1[i]);                // SNP should be flipped
                snpfrag[i].post_hap = post_hapf;
            }else{
                snpfrag[i].post_hap = post_hap;
            }
        }
    }
}

int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr) {
    int i, j, f, s;
    float P_data_H, P_data_Hsw, post_sw;
    int split_occured = 0;
    if (!ERROR_ANALYSIS_MODE){
        for (s = 0; s < clist[k].phased; s++)
            snpfrag[clist[k].slist[s]].csize = 0;
    }
    int curr_component = clist[k].slist[0];
    for (s = 1; s < clist[k].phased; s++) {
        i = clist[k].slist[s]; // i is the current SNP index being considered

        // set probabilities to zero
        P_data_H = 0;
        P_data_Hsw = 0; // switch at i

        //looping over fragments in component c and sum up read probabilities
        for (j = 0; j < clist[k].frags; j++) {
            // normal haplotype
            f = clist[k].flist[j];
            P_data_H += fragment_ll1(Flist,f, HAP, -1, -1);
            // haplotype with switch error starting at i
            P_data_Hsw += fragment_ll1(Flist,f, HAP, -1, i);
        }
        // posterior probability of no switch error
        post_sw = P_data_Hsw - addlogsEXACT(P_data_H, P_data_Hsw);
        //post_sw_total = addlogsEXACT(post_sw_total, post_sw);
        snpfrag[i].post_notsw = subtractlogs(0,post_sw);
        //snpfrag[i].post_notsw_total = subtractlogs(0,post_sw_total);
        // flip the haplotype at this position if necessary
        if (!ERROR_ANALYSIS_MODE){

            if (snpfrag[i].post_notsw < log10(SPLIT_THRESHOLD)) {
                // block should be split here
                curr_component = clist[k].slist[s];
                split_occured = 1;
                //post_sw_total = FLT_MIN;
            }

            snpfrag[i].component = curr_component;
            snpfrag[curr_component].csize++;
        }
    }

    if (!ERROR_ANALYSIS_MODE){
        // subtract the component we started with, which may or may not exist anymore
        (*components_ptr)--;
        for (s = 0; s < clist[k].phased; s++) {
            i = clist[k].slist[s]; // i is the current SNP index being considered
            if (snpfrag[i].csize > 1){
                // this is the head SNP of a component size 2 or greater
                (*components_ptr)++;
            }else{
                snpfrag[i].csize = snpfrag[snpfrag[i].component].csize;
                if (snpfrag[i].csize <= 1) snpfrag[i].component = -1;
            }
        }
    }

    return split_occured;
}
