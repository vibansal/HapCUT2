
#define DEBUG 0

#include "maxcut_functions.c"
#include "maxcut_lr.c"
#include "common.h"
#include <float.h>

extern int MAX_ITER;
extern int CONVERGE;
extern float HOMOZYGOUS_PRIOR;
extern float THRESHOLD;
extern float SPLIT_THRESHOLD;
extern int VERBOSE;
extern int SPLIT_BLOCKS;
extern int SPLIT_BLOCKS_MAXCUT;
extern int* iters_since_improvement;
extern int* iters_since_split;
extern int REFHAP_HEURISTIC;
extern int ERROR_ANALYSIS_MODE;
extern int HTRANS_MAXBINS;
extern int HTRANS_BINSIZE;
/* edge weights
   weight 334 373 hap 10 alleles 01 W 3.048192
   weight 334 375 hap 11 alleles 00 W 3.523493
   weight 335 298 hap 00 alleles 01 W -2.702092
   weight 335 301 hap 00 alleles 01 W -3.098797
   we are trying to find negative weight cuts, lower the better....

 */

//int NEW_CODE = 1; // likelihood based, max-cut calculated using partial likelihoods

/****** CODE TO FIND MAX CUT FOR EACH COMPONENT *************/
int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter, int* components_ptr);
float compute_goodcut(struct SNPfrags* snpfrag, char* hap, int* slist, struct BLOCK* component, struct fragment* Flist, int algo);
void prune_snps(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1);
void refhap_heuristic(int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1);
int split_block(char* HAP, struct BLOCK* clist, int k, struct fragment* Flist, struct SNPfrags* snpfrag, int* components_ptr);
void split_blocks_old(char* HAP, struct BLOCK* clist, int components, struct fragment* Flist, struct SNPfrags* snpfrag);
//int maxcut_split_blocks(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter);
void improve_hap(char* HAP, struct BLOCK* clist, int components, int snps, int fragments, struct fragment* Flist, struct SNPfrags* snpfrag);
//float addlogs(float a, float b);

void single_variant_flips_LL(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1) {
    int t = 0, i = 0, f = 0;
    double delta = 0;
    float mec_ll;
    float chim_ll;

    for (t = 0; t < clist[k].phased; t++) {
        for (i = 0; i < snpfrag[slist[t]].frags; i++) {
            f = snpfrag[slist[t]].flist[i];
            calculate_fragscore(Flist, f, HAP1, &mec_ll, &chim_ll);
            delta += mec_ll;
        }

        // flip the single SNP allele between the two haplotypes
        if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
        else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';

        for (i = 0; i < snpfrag[slist[t]].frags; i++) {
            f = snpfrag[slist[t]].flist[i];
            calculate_fragscore(Flist, f, HAP1, &mec_ll, &chim_ll);
            delta -= mec_ll;
        }

        if (delta > 0) // newMEC is not better than original MEC so we need to revert back to oldMEC for each fragment
        {
            if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
            else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';
        } else {
            clist[k].MEC += delta;
        }
    }
}


// MEC for likelihood function is -1*LL -> preserves comparisons

void single_variant_flips(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1) {
    int t = 0, i = 0, f = 0;
    float newscore = clist[k].bestMEC;

    for (t = 0; t < clist[k].phased; t++) {
        if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
        else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';

        for (i = 0; i < snpfrag[slist[t]].frags; i++) {
            f = snpfrag[slist[t]].flist[i];
            newscore -= Flist[f].currscore;
            update_fragscore(Flist, f, HAP1);
            newscore += Flist[f].currscore;
        }

        if (newscore > clist[k].bestMEC) // newMEC is not better than original MEC so we need to revert back to oldMEC for each fragment
        {
            if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
            else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';
            for (i = 0; i < snpfrag[slist[t]].frags; i++) {
                f = snpfrag[slist[t]].flist[i];
                update_fragscore(Flist, f, HAP1);
            }
        } else {
            clist[k].MEC = newscore;
            clist[k].bestMEC = newscore;
        }
    }
}

/**************** DETERMINISTIC MAX_CUT MEC IMPLEMENTATION *********************************************************//////

// function called from hapcut.c for each component...

int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter, int* components_ptr) {
    clist[k].split = 0;

    if (iters_since_improvement[clist[k].offset] > CONVERGE && (!SPLIT_BLOCKS_MAXCUT || iters_since_split[clist[k].offset] > CONVERGE)) {
        return 1; // return 1 only if converged
    }

    int i = 0, j = 0, t = 0, first_in, first_out, count1, count2;
    float cutvalue, post;
    /*
       i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++)
       {
       if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; }
       }*/
    // slist is copied from clist[k].slist
    for (t = 0; t < clist[k].phased; t++) slist[t] = clist[k].slist[t]; // new way to determine slist

    // not required but we do it to avoid errors
    clist[k].MEC = 0;
    clist[k].LL = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].MEC += Flist[clist[k].flist[i]].currscore;
        clist[k].LL += Flist[clist[k].flist[i]].ll;
    }
    clist[k].bestMEC = clist[k].MEC;
    // evaluate the impact of flipping of each SNP on MEC score, loop needed to improve MEC score
    single_variant_flips(Flist, snpfrag, clist, k, slist, HAP1);

    clist[k].MEC = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].MEC += Flist[clist[k].flist[i]].currscore;
    }
    clist[k].bestMEC = clist[k].MEC;
    //if (fabsf(newscore-clist[k].MEC) >= 0.001) fprintf(stdout,"old %f %f MECcheck\n",newscore,clist[k].MEC);

    cutvalue = 10;

    if (clist[k].MEC > 0) cutvalue = compute_goodcut(snpfrag, HAP1, slist, &clist[k], Flist, MINCUTALGO);
    // flip the subset of columns in slist with positive value
    //if (cutvalue <= 3 || MINCUTALGO == 2) { //getchar();
    for (i = 0; i < clist[k].phased; i++) {
        if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
        else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
    }

    clist[k].bestMEC = clist[k].MEC;
    clist[k].MEC = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].MEC += Flist[clist[k].flist[i]].currscore;
    }

    post = (-1.0*clist[k].bestMEC)-(addlogs((-1.0*clist[k].bestMEC), (-1.0*clist[k].MEC)));

    if ((iters_since_improvement[clist[k].offset] <= CONVERGE && clist[k].MEC >= clist[k].bestMEC)
        ||(SPLIT_BLOCKS_MAXCUT && post >= log10(SPLIT_THRESHOLD))){// revert to old haplotype
        // we aren't splitting blocks, or cut is above threshold
        // flip back the SNPs in the cut
        iters_since_improvement[clist[k].offset]++;
        for (i = 0; i < clist[k].phased; i++) {
            if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
            else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
        }
        clist[k].MEC = 0;
        for (i = 0; i < clist[k].frags; i++) {
            update_fragscore(Flist, clist[k].flist[i], HAP1);
            clist[k].MEC += Flist[clist[k].flist[i]].currscore;
        }
    }else if (SPLIT_BLOCKS_MAXCUT && iters_since_improvement[clist[k].offset] > CONVERGE && post < log10(SPLIT_THRESHOLD)){
        // solution has converged so we are trying to split blocks
        // this cut is under the threshold so we cut it out as a separate block
        clist[k].split = 1;
        iters_since_split[clist[k].offset] = 0;

        first_in = -1;  // first element in the cut
        first_out = -1; // first element not in the cut
        count1 = 0; count2 = 0; // counts for size of each side of cut

        for (j = 0; j < clist[k].phased; j++){

            if (slist[j] > 0){
                // j is in the cut
                if (first_in == -1){
                    first_in = slist[j];
                    snpfrag[first_in].csize = 0;
                }
                snpfrag[slist[j]].component = first_in;
                snpfrag[first_in].csize ++;
                count1++;
            }else{
                // j is not in the cut
                if (first_out == -1){
                    first_out = clist[k].slist[j];
                    snpfrag[first_out].csize = 0;
                }
                snpfrag[clist[k].slist[j]].component = first_out;
                snpfrag[first_out].csize ++;
                count2++;
            }
        }

        iters_since_improvement[first_in] = CONVERGE+1;
        iters_since_improvement[first_out] = CONVERGE+1;

        for (j = 0; j < clist[k].phased; j++){
            snpfrag[clist[k].slist[j]].csize = snpfrag[snpfrag[clist[k].slist[j]].component].csize;
            if (snpfrag[clist[k].slist[j]].csize <= 1){
                snpfrag[clist[k].slist[j]].component = -1;
            }
        }

        if (count1 > 1 && count2 > 1)
            (*components_ptr)++;
        else if (count1 == 1 && count2 == 1)
            (*components_ptr)--;

    }else{
        clist[k].bestMEC = clist[k].MEC; // update current haplotype
        iters_since_improvement[clist[k].offset] = 0;
    }

    if (VERBOSE && iter > 0 && clist[k].MEC > 0) fprintf(stdout, "component %d offset %d length %d phased %d  calls %d MEC %0.1f cutvalue %f bestMEC %0.2f\n", k, clist[k].offset, clist[k].length, clist[k].phased, clist[k].calls, clist[k].MEC, cutvalue, clist[k].bestMEC);

    if (SPLIT_BLOCKS_MAXCUT && iters_since_improvement[clist[k].offset] > CONVERGE && !clist[k].split)
        iters_since_split[clist[k].offset]++;

    // return 0 to specify that this component hasn't converged.
    return 0;
}

/********* THIS IS THE MAIN FUNCTION FOR HAPLOTYPE ASSEMBLY USING ITERATIVE MAX CUT computations **************/

float compute_goodcut(struct SNPfrags* snpfrag, char* hap, int* slist, struct BLOCK* component, struct fragment* Flist, int algo) {
    // given a haplotype 'hap' and a fragment matrix, find a cut with positive score
    int totaledges = 0, i = 0, j = 0, k = 0, l = 0, f = 0;
    int wf = 0; //if (drand48() < 0.5) wf=1;
    float W = 0;
    int N = component->phased;

    /* CODE TO set up the read-haplotype consistency graph */
    for (i = 0; i < N; i++) {
        snpfrag[slist[i]].tedges = 0;
        k = -1;
        // edges contain duplicates in sorted order, but tedges is unique count of edges
        for (j = 0; j < snpfrag[slist[i]].edges; j++) {
            if (k != snpfrag[slist[i]].elist[j].snp) {
                snpfrag[slist[i]].tedges++;
                k = snpfrag[slist[i]].elist[j].snp;
            }
        }
    }
    for (i = 0; i < N; i++) {
        snpfrag[slist[i]].tedges = 0;
        k = -1;
        for (j = 0; j < snpfrag[slist[i]].edges; j++) {
            if (k != snpfrag[slist[i]].elist[j].snp) {
                snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].snp = snpfrag[slist[i]].elist[j].snp;
                k = snpfrag[slist[i]].elist[j].snp;
                W = (float) edge_weight(hap, slist[i], k, snpfrag[slist[i]].elist[j].p, Flist, snpfrag[slist[i]].elist[j].frag);
                if (wf == 0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls - 1; //(fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1);
                snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = W;
                snpfrag[slist[i]].tedges++;
                totaledges++;
            } else if (k == snpfrag[slist[i]].elist[j].snp) {
                W = (float) edge_weight(hap, slist[i], k, snpfrag[slist[i]].elist[j].p, Flist, snpfrag[slist[i]].elist[j].frag);
                if (wf == 0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls - 1; //(fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1);
                snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges - 1].w += W;
            }
        }
    }
    /* CODE TO find 'K' biggest edges in MEC graph, negative weight edges in graph  */
    int K = 5;
    int smallest = 0;
    float smallw = 1000;
    if (totaledges / 2 < K) K = totaledges / 2;
    EDGE* edgelist = (EDGE*) malloc(sizeof (EDGE) * K);
    j = 0;
    i = 0;
    k = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < snpfrag[slist[i]].tedges; j++) {
            if (k < K) {
                edgelist[k].s = slist[i];
                edgelist[k].t = snpfrag[slist[i]].telist[j].snp;
                edgelist[k].w = snpfrag[slist[i]].telist[j].w;
                if (edgelist[k].w < smallw) {
                    smallest = k;
                    smallw = edgelist[k].w;
                }
                k++;
            } else {
                if (snpfrag[slist[i]].telist[j].w > smallw) {
                    edgelist[smallest].s = slist[i];
                    edgelist[smallest].t = snpfrag[slist[i]].telist[j].snp;
                    edgelist[smallest].w = snpfrag[slist[i]].telist[j].w;
                    smallw = 1000;
                    for (l = 0; l < K; l++) {
                        if (edgelist[l].w < smallw) {
                            smallest = l;
                            smallw = edgelist[l].w;
                        }
                    }
                }
            }
        }
    }

    /* CODE TO set up the read-haplotype consistency graph */


    // edge contraction algorithm: merge vertices until only two nodes left or total edge weight of graph is negative
    int startnode = (int) (drand48() * N);
    if (startnode == N) startnode--;
    int secondnode = -1; // root of 2nd cluster initially not there
    // chose a positive edge to initialize the two clusters and run this algorithm $O(m)$ times for each block
    // a negative weight cut should have at least one negative edge or if there is no negative weight edge, the edge with lowest weight

    for (i = 0; i < N; i++) {
        snpfrag[slist[i]].revmap = i;
    }
    int V = N;
    float curr_cut = 0, best_cut = 10000;
    int snp_add;
    int c1 = 0, c2 = 0;
    char* bestmincut;
    //	int size_small,best_small=0,secondlast=0,last=0;
    int iter = 0, maxiter = N / 10;
    if (N / 10 < 1) maxiter = 1;
    if (maxiter >= MAXCUT_ITER && MAXCUT_ITER >= 1) maxiter = MAXCUT_ITER; // added march 13 2013

    int fixheap = 0;
    PHEAP pheap;
    pinitheap(&pheap, N); // heap for maxcut

    /*****************************Maintain two clusters and add each vertex to one of these two ******************/
    bestmincut = (char*) malloc(N);
    for (i = 0; i < N; i++) bestmincut[i] = '0';
    //for (iter=0;iter<totaledges*(int)(log2(totaledges));iter++)
    for (iter = 0; iter < maxiter + K; iter++) {
        pheap.length = N - 2;
        V = N - 2;
        if (iter < K) {
            startnode = edgelist[iter].s;
            secondnode = edgelist[iter].t;
            if (DEBUG) fprintf(stdout, " edge sel %d %d %f \n", startnode, secondnode, edgelist[iter].w);
        } else {
            if (NEW_CODE == 0 || (NEW_CODE == 1 && drand48() < 0.5)) {
                i = (int) (drand48() * totaledges - 0.0001);
                j = 0;
                while (i >= snpfrag[slist[j]].tedges) {
                    i -= snpfrag[slist[j]].tedges;
                    j++;
                }
                startnode = slist[j];
                secondnode = snpfrag[slist[j]].telist[i].snp;
                if (snpfrag[slist[j]].telist[i].w >= 1) continue;
            } else {
                // find node with high MEC score, initialize as startnode
                j = (int) (drand48() * N);
                if (j >= N) j = N - 1;
                startnode = slist[j];
                secondnode = -1;
                pheap.length = N - 1;
                V = N - 1;
            }
        }

        for (i = 0; i < N; i++) snpfrag[slist[i]].parent = slist[i];
        // new code added for heap based calculation
        for (i = 0; i < N; i++) snpfrag[slist[i]].score = 0;
        j = 0; // heap only has N-2 elements (startnode and secondnode are not there)
        for (i = 0; i < N; i++) {
            if (slist[i] != startnode && slist[i] != secondnode) {
                pheap.elements[j] = i;
                snpfrag[slist[i]].heaploc = j;
                j++;
            }
        }
        //		for (i=0;i<N;i++) fprintf(stdout,"heaploc %d %d %d-%d\n",slist[i],snpfrag[slist[i]].heaploc,startnode,secondnode);

        for (i = 0; i < component->frags; i++) {
            f = component->flist[i];
            Flist[f].scores[0] = 0.0;
            Flist[f].scores[1] = 0.0;
            Flist[f].scores[2] = 0.0;
            Flist[f].scores[3] = 0.0;

            Flist[f].htscores[0] = 0.0;
            Flist[f].htscores[1] = 0.0;
            Flist[f].htscores[2] = 0.0;
            Flist[f].htscores[3] = 0.0;

            Flist[f].init = '1';
        }
        if (NEW_CODE == 1) // long reads, likelihood based
        {
            init_fragment_scores(snpfrag, Flist, hap, startnode, secondnode);
        } else {
            init_neighbor_scores(snpfrag, startnode, Flist, hap, 1);
            init_neighbor_scores(snpfrag, secondnode, Flist, hap, -1);
        }

        pbuildmaxheap(&pheap, snpfrag, slist);
        //V = N-2;
        while (V > 0) // more than two clusters, this loop is O(N^2)
        {
            snp_add = pheap.elements[0];
            premovemax(&pheap, snpfrag, slist);
            fixheap = 0;
            //if (N < 30) fprintf(stdout,"standard best score %f snp %d %d V %d\n",snpfrag[slist[snp_add]].score,snp_add,slist[snp_add],V);
            if (snpfrag[slist[snp_add]].score > 0) snpfrag[slist[snp_add]].parent = startnode;
            else if (snpfrag[slist[snp_add]].score < 0) {
                if (secondnode < 0) {
                    secondnode = slist[snp_add];
                    //fprintf(stderr,"secondnode found %d %f V %d N %d\n",secondnode,snpfrag[slist[snp_add]].score,V,N);
                }

                snpfrag[slist[snp_add]].parent = secondnode;
            } else if (secondnode < 0) secondnode = slist[snp_add];
            else // score is 0
            {
                if (drand48() < 0.5) snpfrag[slist[snp_add]].parent = startnode;
                else snpfrag[slist[snp_add]].parent = secondnode;
            }
            V--;

            if (NEW_CODE == 1) {
                update_fragment_scores(snpfrag, Flist, hap, startnode, secondnode, slist[snp_add], &pheap, slist);
                for (i = 0; i < N; i++) {
                    if (DEBUG) fprintf(stdout, "score %d %f hap %c \n", slist[i], snpfrag[slist[i]].score, hap[slist[i]]);
                }
                if (DEBUG) fprintf(stdout, "init frag-scores %d...%d new node added %d parent %d\n\n", startnode, secondnode, slist[snp_add], snpfrag[slist[snp_add]].parent);
            } else update_neighbor_scores(snpfrag, slist[snp_add], startnode, secondnode, Flist, hap, &pheap, slist);

            if (fixheap == 1) pbuildmaxheap(&pheap, snpfrag, slist);
        }
        if (secondnode == -1) continue; // cut is empty, so we should ignore this cut

        // compute score of the cut computed above
        for (i = 0; i < N; i++) {
            if (snpfrag[slist[i]].parent == startnode) snpfrag[slist[i]].parent = 0;
            else snpfrag[slist[i]].parent = 1;
        }
        c1 = 0;
        c2 = 0;
        for (i = 0; i < N; i++) {
            if (snpfrag[slist[i]].parent == 0) c1++;
            else c2++;
        }
        //curr_cut = -1*cut_MEC(snpfrag,Flist,hap,slist,N,component);
        //fprintf(stderr,"component %d MEC %f\n",N,curr_cut); //fprintf(stderr,"component %d cut value%f\n",N,curr_cut);
        if (c1 == 0 || c2 == 0) {
            fprintf(stdout, " cut size is 0 red \n");
            exit(0);
        }

        if (NEW_CODE == 1) {
            curr_cut = cut_score(Flist, snpfrag, component, hap);
            // cut score returns difference between likelihood of current haplotype and new haplotype => smaller it is, better the cut
            if (DEBUG) fprintf(stdout, "cut size %d %d %f best %f\n", c1, c2, curr_cut, best_cut);
        } else {
            curr_cut = find_cutvalue(snpfrag, Flist, hap, slist, N);
            if (DEBUG) fprintf(stdout, "cut size %d %d %f best %f\n", c1, c2, curr_cut, best_cut);
        }

        // for new likelihood based cut score, the score of the cut should always be less than 0 since it is difference of the log-likelihoods of old and new haplotypes
        if (curr_cut < best_cut) // negative weight cut is better...
        {
            best_cut = curr_cut;
            for (i = 0; i < N; i++) {
                if (snpfrag[slist[i]].parent == 1) bestmincut[i] = '1';
                else bestmincut[i] = '0';
            }
        }
        //exit(0);

    }
    for (i = 0; i < N; i++) {
        if (bestmincut[i] == '1') slist[i] = -1 * slist[i] - 1;
    }
    free(bestmincut);
    free(pheap.elements);
    free(edgelist);
    return best_cut;
}

// sets snpfrag[i].prune_status:
// 0 indicates not pruned
// 1 indicates pruned (leave unphased in output)
// 2 indicates called homozygous 00
// 3 indicates called homozygous 11
// the posterior probability cutoff is defined by THRESHOLD

void prune_snps(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) {
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
        if (post_00 > log10(THRESHOLD)){
            snpfrag[i].prune_status = 2; // 2 specifies 00 homozygous
            snpfrag[i].post_hap = post_00;
        }else if (post_11 > log10(THRESHOLD)){
            snpfrag[i].prune_status = 3; // 3 specifies 11 homozygous
            snpfrag[i].post_hap = post_11;
        }else if (post_hapf > log10(THRESHOLD)){
            flip(HAP1[i]);                // SNP should be flipped
            snpfrag[i].post_hap = post_hapf;
        }else if (post_hap < log10(THRESHOLD)){
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
        snpfrag[i].prune_status = 0;

        if (good[i] == bad[i]){
            snpfrag[i].pruned_refhap_heuristic = 1; // this isn't used to prune, just recorded.

            if (REFHAP_HEURISTIC)
                snpfrag[i].prune_status = 1;
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
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag, int snps){
    int i,j,k,snp_ix,mate1match,mate2match, on_mate1, bin, skip;
    float Q =0;
    int Qcount=0;
    // consistent counts, inconsistent counts
    float* p_htrans = calloc(HTRANS_MAXBINS,sizeof(float));
    int* cons       = calloc(HTRANS_MAXBINS,sizeof(int));
    int* incons     = calloc(HTRANS_MAXBINS,sizeof(int));

    // count the number of consistent vs inconsistent for a given insert size bin
    for (i = 0; i < fragments; i++){
        skip = 0;
        // keep things very simple by only sampling 1-snp mates
        if (Flist[i].calls != 2 || Flist[i].mate2_ix == -1){
            continue;
        }

        mate1match = -1;
        mate2match = -1;

        for (j = 0; j < Flist[i].blocks; j++){
            if (skip) break;

            for (k = 0; k < Flist[i].list[j].len; k++) {

                snp_ix = Flist[i].list[j].offset + k;
                on_mate1 = (snp_ix < Flist[i].mate2_ix); // which mate are we on?

                if (HAP[snp_ix] == '-' || snpfrag[snp_ix].prune_status == 1){
                    skip = 1;
                    break;
                }

                // for computing average quality score
                Q += Flist[i].list[j].pv[k];
                Qcount ++;

                if (on_mate1){
                    // on mate 1
                    if(HAP[snp_ix] == Flist[i].list[j].hap[k])
                        mate1match = 1;
                    else
                        mate1match = 0;
                }else{
                    // on mate 2
                    if(HAP[snp_ix] == Flist[i].list[j].hap[k])
                        mate2match = 1;
                    else
                        mate2match = 0;
                }

                bin = Flist[i].isize / HTRANS_BINSIZE;
                if (mate1match == mate2match){
                    cons[bin]++;
                }else{
                    incons[bin]++;
                }
            }
        }
    }


    // compute the average quality score Q, and use it to estimate P(flip) ignoring h-trans

    Q /= Qcount;
    float log_denom = 2*log10(2*Q-1); // denominator for htrans estimate
    float numval    = 2*Q*(1-Q);
    float emp; // empirically estimate mate flip probability

    for (i = 0; i < HTRANS_MAXBINS; i++){
        if (cons[i] + incons[i] == 0)
            continue;

        emp = (incons[i]/(cons[i] + incons[i]));
        p_htrans[i] = log10(emp - numval) - log_denom;
    }

    // assign the probabilities to fragments based in insert size
    for (i = 0; i < fragments; i++){
        // no mate
        if (Flist[i].mate2_ix == -1){
            Flist[i].htrans_prob = -80;
            continue;
        }

        Flist[i].htrans_prob = p_htrans[Flist[i].isize/HTRANS_BINSIZE];

    }

    free(cons);
    free(incons);
    free(p_htrans);

    return 0;
}
