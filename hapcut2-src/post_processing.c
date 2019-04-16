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


void calculate_hetLL(struct fragment* Flist, int fragments,struct SNPfrags* snpfrag, int snps,char* HAP1);
void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous);
int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr);

// not used currentyl
void calculate_hetLL(struct fragment* Flist, int fragments,struct SNPfrags* snpfrag, int snps,char* HAP1)
{
	float log_half = log10(0.5);
	int i=0,j=0,k=0,snp_ix=0,f=0;
    	for (i = 0; i < snps; i++) snpfrag[i].hetLL = 0.0;
	for (f=0;f<fragments;f++)
	{
		for (j = 0; j < Flist[f].blocks; j++) 
		{
		    for (k = 0; k < Flist[f].list[j].len; k++) 
		    {
			snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs
			if (HAP1[snp_ix] == '-') continue;
			snpfrag[snp_ix].hetLL += log_half;
		    }
		}
	}
}

void unphased_optim(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) 
{   
    float log_half = log10(0.5);
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, total,ll,max,delta;
    fprintf(stderr,"unphased genotype likelihood calculation for data %d \n",snps);

    for (i = 0; i < snps; i++) {
        // "unmask" indel variant
        if (SNVS_BEFORE_INDELS && (strlen(snpfrag[i].allele0) != 1 || strlen(snpfrag[i].allele1) != 1)) HAP1[i] = '0';
	//if (HAP1[i] == '-') fprintf(stdout,"SNP %d %d %s %s %s %d\n",i,snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1,snpfrag[i].genotypes,snpfrag[i].frags);
        // we only care about positions that are haplotyped
        //if (!(HAP1[i] == '1' || HAP1[i] == '0')) continue;

        temp1 = HAP1[i];
        P_data_H = 0; P_data_Hf = 0; ll =0;
	for (j=0;j<5;j++) snpfrag[i].PGLL[j] = 0;
	snpfrag[i].hetLL = 0.0; 

        for (j = 0; j < snpfrag[i].frags; j++) {

            f = snpfrag[i].flist[j];
            P_data_H += fragment_ll1(Flist, f, HAP1, -1, -1);
            flip(HAP1[i]);
            P_data_Hf += fragment_ll1(Flist,f, HAP1, -1, -1);

            HAP1[i] = '0'; snpfrag[i].PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1); // added 03/12/2018, genotype = 0|0
            HAP1[i] = '1'; snpfrag[i].PGLL[3] += fragment_ll1(Flist,f, HAP1, i, -1); // genotype = 1|1
            //return haplotype to original value
            HAP1[i] = temp1;
	    ll = fragment_ll1(Flist,f, HAP1, i, 10); // variant 'i' ignored for likelihood calculation if last parameter = 10
	    if (ll < 0) snpfrag[i].PGLL[4] +=ll;  // genotype likelihood where variant is 'unphased' or missing
	    snpfrag[i].hetLL += log_half;
        }
	snpfrag[i].PGLL[1] = P_data_H; snpfrag[i].PGLL[2] = P_data_Hf; // genotypes = 0|1 and 1|0, actually original and flipped

	max = P_data_H; 
	if (P_data_Hf > max) max = P_data_Hf; 
	delta = snpfrag[i].hetLL+snpfrag[i].PGLL[4]-max;

	if (HAP1[i] == '0') temp1 = '1'; else temp1 = '0';
	fprintf(stdout,"SNP %d %d %s %s ",i,snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1);
	fprintf(stdout," %c|%c cov %d hetLL %0.2f %0.2f ",HAP1[i],temp1,snpfrag[i].frags,snpfrag[i].hetLL,snpfrag[i].hetLL+snpfrag[i].PGLL[4]);
	fprintf(stdout,"hets= %0.2f %0.2f ",snpfrag[i].PGLL[1],snpfrag[i].PGLL[2]);
	fprintf(stdout,"homs= %0.2f %0.2f delta %0.2f ",snpfrag[i].PGLL[0],snpfrag[i].PGLL[3],delta);
	if (snpfrag[i].PGLL[0] > max) fprintf(stdout,"ref-hom %0.2f ",snpfrag[i].PGLL[0]-max);
	if (snpfrag[i].PGLL[3] > max) fprintf(stdout,"alt-hom %0.2f ",snpfrag[i].PGLL[3]-max);
	if (P_data_Hf > P_data_H) fprintf(stdout,"flip-phase ");
	else if (delta > 0) fprintf(stdout,"unphased ");
	fprintf(stdout,"\n");

    }
}

void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous) {
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, total, log_hom_prior, log_het_prior;
    float post_hap, post_hapf, post_00, post_11,ll;

    for (i = 0; i < snps; i++) {

        // "unmask" indel variant
        if (SNVS_BEFORE_INDELS && (strlen(snpfrag[i].allele0) != 1 || strlen(snpfrag[i].allele1) != 1)) HAP1[i] = '0';

        // we only care about positions that are haplotyped
        if (!(HAP1[i] == '1' || HAP1[i] == '0')) continue;

        // hold on to original haplotype values
        temp1 = HAP1[i];
        // set probabilities to zero
        P_data_H = 0;
        P_data_Hf = 0;
        P_data_H00 = 0;
        P_data_H11 = 0;
	for (j=0;j<5;j++) snpfrag[i].PGLL[j] = 0; 

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
            HAP1[i] = '0'; snpfrag[i].PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1); // added 03/12/2018
            HAP1[i] = '1'; snpfrag[i].PGLL[3] += fragment_ll1(Flist,f, HAP1, i, -1); // homozygous genotypes 00,11
            //return haplotype to original value
            HAP1[i] = temp1;
	    ll = fragment_ll1(Flist,f, HAP1, i, 10); // variant 'i' ignored for likelihood calculation if last parameter = 10
	    if (ll < 0) snpfrag[i].PGLL[4] +=ll;
        }
	snpfrag[i].PGLL[1] = P_data_H; snpfrag[i].PGLL[2] = P_data_Hf;


        if (call_homozygous){

            // get prior probabilities for homozygous and heterozygous genotypes
            log_hom_prior = snpfrag[i].homozygous_prior;
            log_het_prior = subtractlogs(log10(0.5), log_hom_prior);

            // denominator of posterior probabilities;
            // sum of all 4 data probabilities times their priors
            total = addlogs(
                    addlogs((log_het_prior + P_data_H), (log_het_prior + P_data_Hf)),
                    addlogs((log_hom_prior + P_data_H00), (log_hom_prior + P_data_H11)));

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

            total = addlogs((log_het_prior + P_data_H), (log_het_prior + P_data_Hf));

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
        post_sw = P_data_Hsw - addlogs(P_data_H, P_data_Hsw);
        //post_sw_total = addlogs(post_sw_total, post_sw);
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
