/*
 Author: Peter Edge
 */

#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

extern float THRESHOLD;
extern int DISCRETE_PRUNING;
extern int ERROR_ANALYSIS_MODE;
extern int HTRANS_MAXBINS;
extern int HTRANS_BINSIZE;
extern float HOMOZYGOUS_PRIOR;
extern float SPLIT_THRESHOLD;
extern int SPLIT_BLOCKS;
extern int HTRANS_READ_LOWBOUND;
extern char HTRANS_DATA_OUTFILE[10000];
extern int HTRANS_MAX_WINDOW;
extern int SNVS_BEFORE_INDELS;

float HIC_EM_THRESHOLD = 0.99; // use a strict-ish threshold for the HiC haplotype SNPs that we'll estimate h-trans from

void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous);
void discrete_pruning(int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1);
int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr);
void split_blocks_old(char* HAP, struct BLOCK* clist, int components, struct fragment* Flist, struct SNPfrags* snpfrag);
//int maxcut_split_blocks(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter);
void improve_hap(char* HAP, struct BLOCK* clist, int components, int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag);
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag);

void likelihood_pruning(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, int call_homozygous) {
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, total, log_hom_prior, log_het_prior;
    float post_hap, post_hapf, post_00, post_11;

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

        //looping over fragments overlapping i and sum up read probabilities
        for (j = 0; j < snpfrag[i].frags; j++) {

            f = snpfrag[i].flist[j];
            // normal haplotypes

            P_data_H += fragment_ll(Flist, f, HAP1, -1, -1);
            // haplotypes with i flipped
            flip(HAP1[i]);
            P_data_Hf += fragment_ll(Flist, f, HAP1, -1, -1);

            if (call_homozygous){
                // haplotypes with i homozygous 00
                HAP1[i] = '0';
                P_data_H00 += fragment_ll(Flist, f, HAP1, i, -1);
                // haplotypes with i homozygous 11
                HAP1[i] = '1';
                P_data_H11 += fragment_ll(Flist, f, HAP1, i, -1);
            }

            //return haplotype to original value
            HAP1[i] = temp1;
        }



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

// the discrete heuristic for pruning SNPs (introduced by RefHap, used by ProbHap)
void discrete_pruning(int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) {
    int i, j, f, k, snp_ix, frag_assign;
    float prob, prob2, p0, p1;
    int* good = (int*) malloc(snps*sizeof(int)); // good[i] is the number of frag bases consistent with phasing at SNP i
    int* bad  = (int*) malloc(snps*sizeof(int)); // bad [i] is the number of frag bases inconsistent with phasing at SNP i

    // compute the optimal assignment of each fragment based on log likelihood
    for (f = 0; f < fragments; f++){
        p0 = 0; p1 = 0;
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs

                // conditions to skip this base
                if (HAP1[snp_ix] == '-') continue;

                prob = QVoffset - (int) Flist[f].list[j].qv[k];
                prob /= 10;
                prob2 = Flist[f].list[j].p1[k];
                // this is likelihood based calculation
                if (HAP1[snp_ix] == Flist[f].list[j].hap[k]) {
                    p0 += prob2;
                    p1 += prob;
                } else {
                    p0 += prob;
                    p1 += prob2;
                }
            }
        }

        if (p0 > p1)
            frag_assign = 1;
        else if (p1 > p0)
            frag_assign = 2;
        else if (drand48() < 0.5)
            frag_assign = 1;
        else
            frag_assign = 2;


        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs

                if (HAP1[snp_ix] == '-') continue;

                if ((frag_assign == 1 && HAP1[snp_ix] == Flist[f].list[j].hap[k])
                  ||(frag_assign == 2 && HAP1[snp_ix] != Flist[f].list[j].hap[k]))
                    good[snp_ix]++;
                else
                    bad[snp_ix]++;
            }
        }
    }

    for (i = 0; i < snps; i++){

        snpfrag[i].pruned_discrete_heuristic = (int)(good[i] == bad[i]);

    }

    free(good); free(bad);
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
            P_data_H += fragment_ll(Flist, f, HAP, -1, -1);
            // haplotype with switch error starting at i
            P_data_Hsw += fragment_ll(Flist, f, HAP, -1, i);
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


// estimate probabilities of h-trans to feed back into HapCUT algorithm
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag){

    float* MLE_sum   = calloc(HTRANS_MAXBINS,sizeof(float));
    float* MLE_count = calloc(HTRANS_MAXBINS,sizeof(float));
    float* p_htrans      = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_sum   = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_count = calloc(HTRANS_MAXBINS,sizeof(float));

    int i=0,j=0,k=0,f=0,bin;
    int i_minus = 0, i_plus = 0;
    int e_window_size = HTRANS_BINSIZE; //track the effective window size

    // consistent counts, inconsistent counts
    int i1=-1, i2=-1;
    char a1='-', a2='-', h1='-', h2='-';
    float q1=0,q2=0;
    // count the number of consistent vs inconsistent for a given insert size bin
    for (f = 0; f < fragments; f++){

        // consider mate pairs only
        if (Flist[f].mate2_ix == -1 || Flist[f].isize == -1){
            continue;
        }

        // insert size bin
        bin = Flist[f].isize / HTRANS_BINSIZE;

        if (bin < 0){
            fprintf_time(stderr,"ERROR: bin less than 0");
            exit(1);
        }
        if (bin >= HTRANS_MAXBINS){
            fprintf_time(stderr,"ERROR: bin greater than HTRANS_MAXBINS");
            exit(1);
        }

        // keep things very simple by only sampling 1-snp mates
        if (Flist[f].calls < 2){
            continue;
        }

        i1 = Flist[f].list[0].offset;
        a1 = Flist[f].list[0].hap[0];
        q1 = Flist[f].list[0].pv[0];

        i2 = Flist[f].mate2_ix;
        for (j=0; j<Flist[f].blocks; j++){
            for (k=0; k<Flist[f].list[j].len; k++){
                if (Flist[f].list[j].offset+k == Flist[f].mate2_ix){
                    a2 = Flist[f].list[j].hap[k];
                    q2 = Flist[f].list[j].pv[k];
                    break;
                }
            }
        }

        h1 = HAP[i1];
        h2 = HAP[i2];

        if (h1 == '-' || h2 == '-'
         || snpfrag[i1].post_hap < log10(HIC_EM_THRESHOLD)|| snpfrag[i2].post_hap < log10(HIC_EM_THRESHOLD)
         || (snpfrag[i1].bcomp != snpfrag[i2].bcomp)
         || (snpfrag[i1].bcomp == -1)
         || (snpfrag[i2].bcomp == -1)
         || (snpfrag[i1].frags == 1)
         || (snpfrag[i2].frags == 1)){
            continue;
        }

        MLE_count[bin]++;

        assert(i1 < Flist[f].mate2_ix);
        assert(i2 == Flist[f].mate2_ix);
        assert (h1 == '1' || h1 == '0');
        assert (h2 == '1' || h2 == '0');
        assert (a1 == '1' || a1 == '0');
        assert (a2 == '1' || a2 == '0');

        if ((a1 == a2) == (h1 == h2)){ // phase match
            MLE_sum[bin] += ((1-q1)*q2 + (1-q2)*q1);
        }else{                         // phase mismatch
            MLE_sum[bin] += ((1-q1)*(1-q2) + q1*q2);
        }
    }


    for (i = 0; i < HTRANS_MAXBINS; i++){
        adj_MLE_count[i] = MLE_count[i];
        adj_MLE_sum[i] = MLE_sum[i];
        i_minus = i;
        i_plus = i;
        e_window_size = HTRANS_BINSIZE; //track the effective window size
        for (j = 0; j< 100000; j++){
            if (adj_MLE_count[i] >= HTRANS_READ_LOWBOUND) break;
            i_minus--;
            i_plus++;
            if (i_minus >= 0){
                adj_MLE_count[i] += MLE_count[i_minus];
                adj_MLE_sum[i] += MLE_sum[i_minus];
                e_window_size += HTRANS_BINSIZE;
            }
            if(i_plus < HTRANS_MAXBINS){
                adj_MLE_count[i] += MLE_count[i_plus];
                adj_MLE_sum[i] += MLE_sum[i_plus];
                e_window_size += HTRANS_BINSIZE;
            }
            if (e_window_size >= HTRANS_MAX_WINDOW) break; // cut off window expansion if it's larger than some amount
        }
    }

    // compute the MLE for each bin
    for (i = 0; i < HTRANS_MAXBINS; i++){
        if (adj_MLE_count[i] == 0 || adj_MLE_sum[i] == 0)
            p_htrans[i] = -80;
        else
            p_htrans[i] = log10(adj_MLE_sum[i] / adj_MLE_count[i]);
    }

    // assign the probabilities to fragments based in insert size
    for (f = 0; f < fragments; f++){
        if (Flist[f].mate2_ix == -1){
            Flist[f].htrans_prob = -80;
            continue;
        }
        Flist[f].htrans_prob = p_htrans[Flist[f].isize / HTRANS_BINSIZE];
    }

    if (strcmp(HTRANS_DATA_OUTFILE,"None") != 0){
        FILE* fp;
        fp = fopen(HTRANS_DATA_OUTFILE, "w");
        // compute the MLE for each bin
        for (i = 0; i < HTRANS_MAXBINS; i++){
            fprintf(fp,"%d-%d\t%f\n",i*HTRANS_BINSIZE,(i+1)*HTRANS_BINSIZE,
                    pow(10,p_htrans[i]));
        }
        fclose(fp);
    }

    free(MLE_sum);
    free(MLE_count);
    free(adj_MLE_sum);
    free(adj_MLE_count);
    free(p_htrans);
    return 0;
}
