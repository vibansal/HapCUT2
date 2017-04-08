
#define DEBUG 0

#include "maxcut_lr.c"
#include "common.h"
#include "like_scores.c" // additional functions for likelihood based calculation
#include <float.h>
#include <assert.h>     /* assert */

extern int MAX_ITER;
extern int CONVERGE;
extern int VERBOSE;

/****** CODE TO FIND MAX CUT FOR EACH COMPONENT *************/
int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1);
float compute_goodcut(struct SNPfrags* snpfrag, char* hap, int* slist, struct BLOCK* component, struct fragment* Flist);

// likelihood function is -1*LL 

void single_variant_flips(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1) {
    int t = 0, i = 0, f = 0;
    float newscore = clist[k].bestSCORE;

    for (t = 0; t < clist[k].phased; t++) {
        if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
        else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';

        for (i = 0; i < snpfrag[slist[t]].frags; i++) {
            f = snpfrag[slist[t]].flist[i];
            newscore -= Flist[f].currscore;
            update_fragscore(Flist, f, HAP1);
            newscore += Flist[f].currscore;
        }

        if (newscore > clist[k].bestSCORE) // newMEC is not better than original MEC so we need to revert back to oldMEC for each fragment
        {
            if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';
            else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';
            for (i = 0; i < snpfrag[slist[t]].frags; i++) {
                f = snpfrag[slist[t]].flist[i];
                update_fragscore(Flist, f, HAP1);
            }
        } else {
            clist[k].SCORE = newscore;
            clist[k].bestSCORE = newscore;
        }
    }
}


float edge_weight(char* hap, int i, int j, char* p, struct fragment* Flist, int f) {
    float q1 = 1, q2 = 1;
    int k = 0, l = 0;
    // new code added so that edges are weighted by quality scores, running time is linear in length of fragment !! 08/15/13 | reduce this
    for (k = 0; k < Flist[f].blocks; k++) {
        for (l = 0; l < Flist[f].list[k].len; l++) {
            if (Flist[f].list[k].offset + l == i) q1 = Flist[f].list[k].pv[l];
            else if (Flist[f].list[k].offset + l == j) q2 = Flist[f].list[k].pv[l];
        }
    }
    float p1 = q1 * q2 + (1 - q1)*(1 - q2);
    float p2 = q1 * (1 - q2) + q2 * (1 - q1);
    if (hap[i] == hap[j] && p[0] == p[1]) return log10(p1 / p2);
    else if (hap[i] != hap[j] && p[0] != p[1]) return log10(p1 / p2);
    else if (hap[i] == hap[j] && p[0] != p[1]) return log10(p2 / p1);
    else if (hap[i] != hap[j] && p[0] == p[1]) return log10(p2 / p1);
    else return 0;
}

/**************** DETERMINISTIC MAX_CUT MEC IMPLEMENTATION *********************************************************//////

// function called from hapcut.c for each component...

int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1) {

    if (clist[k].iters_since_improvement > CONVERGE) {
        return 1; // return 1 only if converged
    }
    int i = 0, j = 0, t = 0, count1, count2;

    float cutvalue;
    /*
       i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++)
       {
       if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; }
       }*/
    // slist is copied from clist[k].slist
    for (t = 0; t < clist[k].phased; t++) slist[t] = clist[k].slist[t]; // new way to determine slist

    // not required but we do it to avoid errors
    clist[k].SCORE = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
    }
    clist[k].bestSCORE = clist[k].SCORE;
    // evaluate the impact of flipping of each SNP on MEC score, loop needed to improve MEC score
    single_variant_flips(Flist, snpfrag, clist, k, slist, HAP1);

    clist[k].SCORE = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
    }
    clist[k].bestSCORE = clist[k].SCORE;

    cutvalue = 10;

    if (clist[k].SCORE > 0) cutvalue = compute_goodcut(snpfrag, HAP1, slist, &clist[k], Flist);

    // flip the subset of columns in slist with positive value
    //if (cutvalue <= 3 || MINCUTALGO == 2) { //getchar();
    for (i = 0; i < clist[k].phased; i++) {
        if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
        else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
    }

    clist[k].bestSCORE = clist[k].SCORE;
    clist[k].SCORE = 0;
    for (i = 0; i < clist[k].frags; i++) {
        update_fragscore(Flist, clist[k].flist[i], HAP1);
        clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
    }

    count1 = 0; count2 = 0; // counts for size of each side of cut
    for (j = 0; j < clist[k].phased; j++){
        if (slist[j] > 0){
            count1++;
        }else{
            count2++;
        }
    }

    if (clist[k].iters_since_improvement <= CONVERGE && clist[k].SCORE >= clist[k].bestSCORE){// revert to old haplotype
        // flip back the SNPs in the cut

        clist[k].iters_since_improvement++;
        for (i = 0; i < clist[k].phased; i++) {
            if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
            else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
        }
        clist[k].SCORE = 0;
        for (i = 0; i < clist[k].frags; i++) {
            update_fragscore(Flist, clist[k].flist[i], HAP1);
            clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
        }
    }else{
        clist[k].bestSCORE = clist[k].SCORE; // update current haplotype
        clist[k].iters_since_improvement = 0;
    }

    if (VERBOSE && clist[k].SCORE > 0) fprintf_time(stdout, "component %d offset %d length %d phased %d LL %0.1f cutvalue %f bestLL %0.2f\n", k, clist[k].offset, clist[k].length, clist[k].phased, -1*clist[k].SCORE, cutvalue, -1*clist[k].bestSCORE);

    // return 0 to specify that this component hasn't converged.
    return 0;
}

/********* THIS IS THE MAIN FUNCTION FOR HAPLOTYPE ASSEMBLY USING ITERATIVE MAX CUT computations **************/

float compute_goodcut(struct SNPfrags* snpfrag, char* hap, int* slist, struct BLOCK* component, struct fragment* Flist) {
    // given a haplotype 'hap' and a fragment matrix, find a cut with positive score
    int totaledges = 0, i = 0, j = 0, k = 0, l = 0, f = 0;
    int wf = 0; //if (drand48() < 0.5) wf=1;
    float W = 0;
    int N = component->phased;
    int iters_since_improved_cut = 0;

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
                if (wf == 0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls - 1;
                snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = W;
                snpfrag[slist[i]].tedges++;
                totaledges++;
            } else if (k == snpfrag[slist[i]].elist[j].snp) {
                W = (float) edge_weight(hap, slist[i], k, snpfrag[slist[i]].elist[j].p, Flist, snpfrag[slist[i]].elist[j].frag);
                if (wf == 0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls - 1;
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
            if (drand48() < 0.5) {
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
        }

        init_fragment_scores(snpfrag, Flist, hap, startnode, secondnode);

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

                    //fprintf_time(stderr,"secondnode found %d %f V %d N %d\n",secondnode,snpfrag[slist[snp_add]].score,V,N);
                }

                snpfrag[slist[snp_add]].parent = secondnode;
            } else if (secondnode < 0) secondnode = slist[snp_add];
            else // score is 0
            {
                if (drand48() < 0.5) snpfrag[slist[snp_add]].parent = startnode;
                else snpfrag[slist[snp_add]].parent = secondnode;
            }
            V--;


            update_fragment_scores(snpfrag, Flist, hap, startnode, secondnode, slist[snp_add], &pheap, slist);

            for (i = 0; i < N; i++) {
                if (DEBUG) fprintf(stdout, "score %d %f hap %c \n", slist[i], snpfrag[slist[i]].score, hap[slist[i]]);
            }
            if (DEBUG) fprintf(stdout, "init frag-scores %d...%d new node added %d parent %d\n\n", startnode, secondnode, slist[snp_add], snpfrag[slist[snp_add]].parent);
            

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

        if (c1 == 0 || c2 == 0) {
            if (DEBUG) fprintf(stdout, " cut size is 0 red \n");
            exit(0);
        }

        curr_cut = cut_score(Flist, snpfrag, component, hap);
        // cut score returns difference between likelihood of current haplotype and new haplotype => smaller it is, better the cut
        if (DEBUG) fprintf(stdout, "cut size %d %d %f best %f\n", c1, c2, curr_cut, best_cut);

        // for new likelihood based cut score, the score of the cut should always be less than 0 since it is difference of the log-likelihoods of old and new haplotypes
        if (curr_cut < best_cut) // negative weight cut is better...
        {
            best_cut = curr_cut;
            for (i = 0; i < N; i++) {
                if (snpfrag[slist[i]].parent == 1) bestmincut[i] = '1';
                else bestmincut[i] = '0';
            }
        }else{
            iters_since_improved_cut++;
        }
        if (iters_since_improved_cut > CONVERGE){
            break;
        }
    }

    for (i = 0; i < N; i++) {
        if (bestmincut[i] == '1') slist[i] = -1 * slist[i] - 1;
    }
    free(bestmincut);
    free(pheap.elements);
    free(edgelist);
    return best_cut;
}
