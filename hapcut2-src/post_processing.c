/*
 Author: Peter Edge
 */

#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

extern int REFHAP_HEURISTIC;
extern int ERROR_ANALYSIS_MODE;
extern int HTRANS_MAXBINS;
extern int HTRANS_BINSIZE;
extern float HOMOZYGOUS_PRIOR;
extern float SPLIT_THRESHOLD;
extern int SPLIT_BLOCKS;
extern int HTRANS_MLE_COUNT_LOWBOUND;
extern char HTRANS_DATA_OUTFILE[10000];
extern int MAX_WINDOW_SIZE;
extern int HIC_STRICT_FILTER;
// sets snpfrag[i].prune_status:
// 0 indicates not pruned
// 1 indicates pruned (leave unphased in output)
// 2 indicates called homozygous 00
// 3 indicates called homozygous 11
// the posterior probability cutoff is defined by THRESHOLD

void prune_snps(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, float threshold);
void refhap_heuristic(int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1);
int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr);
void split_blocks_old(char* HAP, struct BLOCK* clist, int components, struct fragment* Flist, struct SNPfrags* snpfrag);
//int maxcut_split_blocks(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter);
void improve_hap(char* HAP, struct BLOCK* clist, int components, int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag);
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag, float* MLE_sum, float* MLE_count);
int combine_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag, float* MLE_sum, float* MLE_count);

void prune_snps(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1, float threshold) {
    int i, j, f;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, total, log_hom_prior, log_het_prior;
    float post_hap, post_hapf, post_00, post_11;

    for (i = 0; i < snps; i++) {

        // get prior probabilities for homozygous and heterozygous genotypes
        log_hom_prior = snpfrag[i].homozygous_prior;
        log_het_prior = subtractlogs(log10(0.5), log_hom_prior);

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
            // haplotypes with i homozygous 00
            HAP1[i] = '0';
            P_data_H00 += fragment_ll(Flist, f, HAP1, i, -1);
            // haplotypes with i homozygous 11
            HAP1[i] = '1';
            P_data_H11 += fragment_ll(Flist, f, HAP1, i, -1);
            //return haplotype to original value
            HAP1[i] = temp1;
        }

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
        if (post_00 > log10(threshold)){
            snpfrag[i].prune_status = 2; // 2 specifies 00 homozygous
            snpfrag[i].post_hap = post_00;
        }else if (post_11 > log10(threshold)){
            snpfrag[i].prune_status = 3; // 3 specifies 11 homozygous
            snpfrag[i].post_hap = post_11;
        }else if (post_hapf > log10(threshold)){
            flip(HAP1[i]);                // SNP should be flipped
            snpfrag[i].post_hap = post_hapf;
        }else if (post_hap < log10(threshold)){
            if (!REFHAP_HEURISTIC && !ERROR_ANALYSIS_MODE)
                snpfrag[i].prune_status = 1; // remove the SNP entirely
            snpfrag[i].post_hap = post_hap;
        }
    }
}

// the refhap heuristic for pruning SNPs
void refhap_heuristic(int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) {
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

        if (good[i] == bad[i]){
            snpfrag[i].pruned_refhap_heuristic = 1; // this isn't used to prune, just recorded.

            if (REFHAP_HEURISTIC)
                snpfrag[i].prune_status = 1;
        }else{
            snpfrag[i].pruned_refhap_heuristic = 0;
            if (REFHAP_HEURISTIC)
                snpfrag[i].prune_status = 0;
        }
    }

    free(good); free(bad);
}

// create an array of a random permutation of integers 1..size
int* create_randperm(int size){
    int* array = (int*) malloc(size * sizeof(int));
    int i, swap_ix, temp;
    for (i = 0; i < size; i++){
        array[i] = i;
    }
    for (i = size-1; i > 0; i--){
        swap_ix = (int)((double)rand() / ((double)RAND_MAX + 1) * i); //credit to http://c-faq.com/lib/randrange.html
        temp = array[i];
        array[i] = array[swap_ix];
        array[swap_ix] = temp;
    }
    return array;
}

void improve_hap(char* HAP, struct BLOCK* clist, int components, int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag) {

    int iter;
    int i, j, f, c, s, a;
    char temp1;
    float P_data_H, P_data_Hf, P_data_Hsw, new_ll, old_ll;
    int* randperm;
    // get initial log likelihood
    new_ll = 0;
    for (i = 0; i < fragments; i++) {
        new_ll += fragment_ll(Flist, i, HAP, -1, -1);
    }
    // iterate until convergence or 100 iters
    for (iter = 0; iter < 100; iter++) {

        for (c = 0; c < components; c++) {
            // we want to consider snps in random order
            randperm = create_randperm(clist[c].phased);
            for (s = 0; s < clist[c].phased; s++) {
                i = clist[c].slist[randperm[s]]; // i is the current SNP index being considered.

                // skip ahead if positions aren't haplotyped
                if (!(HAP[i] == '1' || HAP[i] == '0')) {
                    continue;
                }

                // hold on to original haplotype values
                temp1 = HAP[i];
                // set probabilities to zero
                P_data_H = 0; P_data_Hf = 0; P_data_Hsw = 0; // switch error at i

                //looping over fragments overlapping i and sum up read probabilities
                for (j = 0; j < clist[c].frags; j++) {

                    //f = snpfrag[i].flist[j];
                    f = clist[c].flist[j];
                    // normal haplotypes
                    P_data_H += fragment_ll(Flist, f, HAP, -1, -1);
                    // haplotype with switch error starting at i
                    P_data_Hsw += fragment_ll(Flist, f, HAP, -1, i);
                    // haplotype with i flipped
                    flip(HAP[i]);
                    P_data_Hf += fragment_ll(Flist, f, HAP, -1, -1);
                    //return haplotype to original value
                    HAP[i] = temp1;
                }

                // flip the haplotype at this position if necessary
                if (P_data_Hsw > P_data_H) {
                    // block should be switched here
                    if (randperm[s] > clist[c].phased/2){
                        // it is less work to move forward and flip SNPs (starting with s)
                        for (a = randperm[s]; a < clist[c].phased; a++) {
                            j = clist[c].slist[a];
                            // need to switch this position
                            flip(HAP[i]);
                        }
                    }else{
                        // it is less work to move backwards and flip all previous SNPs in block
                        for (a = randperm[s]-1; a >= 0; a--) {
                            j = clist[c].slist[a];
                            // need to switch this position
                            flip(HAP[i]);
                        }
                    }
                }else if (P_data_Hf > P_data_H) {
                    flip(HAP[i]);
                }
            }
            free(randperm);
        }

        old_ll = new_ll;
        new_ll = 0;
        for (i = 0; i < fragments; i++) {
            new_ll += fragment_ll(Flist, i, HAP, -1, -1);
        }

        if (new_ll <= old_ll) break;
    }
}

// ideally we'd have a snpfrag field for which frags cross the snp, gaps or not.
// then it would be valid to simply compute the posterior from those alone
// and not the whole component.
void split_blocks_old(char* HAP, struct BLOCK* clist, int components, struct fragment* Flist, struct SNPfrags* snpfrag) {
    int i, j, f, c, s, blocks_since_split=0;
    float P_data_H, P_data_Hsw, post_sw;

    for (c = 0; c < components; c++) {
        //post_sw_total = FLT_MIN;
        blocks_since_split=0;
        for (s = 0; s < clist[c].phased; s++) {

            i = clist[c].slist[s]; // i is the current SNP index being considered
            if (snpfrag[i].prune_status == 1 || HAP[i] == '-' || blocks_since_split < 3){
                blocks_since_split++;
                continue;
            }// set probabilities to zero
            P_data_H = 0;
            P_data_Hsw = 0; // switch at i

            //looping over fragments in component c and sum up read probabilities
            for (j = 0; j < clist[c].frags; j++) {
                // normal haplotype
                f = clist[c].flist[j];
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
                    snpfrag[i].split = 1;
                    blocks_since_split = 0;
                    //post_sw_total = FLT_MIN;
                }else{
                    blocks_since_split++;
                }
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
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag, float* MLE_sum, float* MLE_count){
    int j,k,f,bin,count, matches, joined, block;
    // consistent counts, inconsistent counts
    int i1=-1, i2=-1;
    char a1='-', a2='-', h1='-', h2='-';
    float q1=0,q2=0;
    // count the number of consistent vs inconsistent for a given insert size bin
    for (f = 0; f < fragments; f++){

        if (!Flist[f].use_for_htrans_est) continue;

        // consider mate pairs only
        if (Flist[f].mate2_ix == -1 || Flist[f].isize == -1){
            continue;
        }

        // insert size bin
        bin = Flist[f].isize / HTRANS_BINSIZE;

        if (bin < 0){
            fprintf(stderr,"ERROR: bin less than 0");
            exit(1);
        }
        if (bin >= HTRANS_MAXBINS){
            fprintf(stderr,"ERROR: bin greater than HTRANS_MAXBINS");
            exit(1);
        }

        // mark reads that aren't supported by the phase
        count = 0;
        matches = 0;
        joined = 1;
        block = -1;
        for (j=0; j<Flist[f].blocks; j++){
            if (!joined) break;
            for (k=0; k<Flist[f].list[j].len; k++){
                if (!((Flist[f].list[j].hap[k] == '1' || Flist[f].list[j].hap[k] == '0')
                    &&(HAP[Flist[f].list[j].offset+k] == '1' || HAP[Flist[f].list[j].offset+k] == '0')))
                    continue;
                if (block == -1){
                    block = snpfrag[Flist[f].list[j].offset+k].bcomp;
                }else if (block != snpfrag[Flist[f].list[j].offset+k].bcomp){
                    joined = 0;
                    break;
                }

                count++;
                if (Flist[f].list[j].hap[k] == HAP[Flist[f].list[j].offset+k])
                    matches++;
            }
        }
        if ((!joined) || (count >= 2 && !(matches == 0 || matches == count))){
            Flist[f].hic_strict_filtered = 1;
        }

        // keep things very simple by only sampling 1-snp mates
        if (Flist[f].calls != 2){
            continue;
        }

        i1 = Flist[f].list[0].offset;
        a1 = Flist[f].list[0].hap[0];
        q1 = Flist[f].list[0].pv[0];

        if (Flist[f].blocks == 1){
            i2 = Flist[f].list[0].offset+1;
            a2 = Flist[f].list[0].hap[1];
            q2 = Flist[f].list[0].pv[1];
        }else if (Flist[f].blocks == 2){
            i2 = Flist[f].list[1].offset;
            a2 = Flist[f].list[1].hap[0];
            q2 = Flist[f].list[1].pv[0];
        }else{
            fprintf(stderr,"ERROR: Inconsistent block structure in estimate_htrans_probs");
        }

        h1 = HAP[i1];
        h2 = HAP[i2];

        if (h1 == '-' || h2 == '-'
         || snpfrag[i1].prune_status == 1 || snpfrag[i2].prune_status == 1
         || (snpfrag[i1].bcomp != snpfrag[i2].bcomp)
         || (snpfrag[i1].bcomp == -1)
         || (snpfrag[i2].bcomp == -1)){
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
    return 0;
}

// estimate probabilities of h-trans to feed back into HapCUT algorithm
int combine_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag, float* MLE_sum, float* MLE_count){
    float* p_htrans      = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_sum   = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_count = calloc(HTRANS_MAXBINS,sizeof(float));
    int i_minus = 0, i_plus = 0, i=0, j=0, f=0;
    int e_window_size = HTRANS_BINSIZE; //track the effective window size

    for (i = 0; i < HTRANS_MAXBINS; i++){
        adj_MLE_count[i] = MLE_count[i];
        adj_MLE_sum[i] = MLE_sum[i];
        i_minus = i;
        i_plus = i;
        e_window_size = HTRANS_BINSIZE; //track the effective window size
        for (j = 0; j< 100000; j++){
            if (adj_MLE_count[i] >= HTRANS_MLE_COUNT_LOWBOUND) break;
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
            if (e_window_size >= MAX_WINDOW_SIZE) break; // cut off window expansion if it's larger than some amount
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

    free(adj_MLE_sum);
    free(adj_MLE_count);
    free(p_htrans);
    return 0;
}
