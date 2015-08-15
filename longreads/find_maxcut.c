
#define DEBUG 0

#include "maxcut_functions.c"
#include "maxcut_lr.c"
#include "common.h"

extern int THRESHOLD_TYPE;
extern int CONVERGE;
extern float HOMOZYGOUS_PRIOR;
/* edge weights 
   weight 334 373 hap 10 alleles 01 W 3.048192 
   weight 334 375 hap 11 alleles 00 W 3.523493 
   weight 335 298 hap 00 alleles 01 W -2.702092 
   weight 335 301 hap 00 alleles 01 W -3.098797 
   we are trying to find negative weight cuts, lower the better....

 */

//int NEW_CODE = 1; // likelihood based, max-cut calculated using partial likelihoods

/****** CODE TO FIND MAX CUT FOR EACH COMPONENT *************/

int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter);
float compute_goodcut(struct SNPfrags* snpfrag, char* hap, int* slist, struct BLOCK* component, struct fragment* Flist, int algo);
void prune_snps(char* pruned, float prune_threshold, float homozygous_threshold, int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1);
float addlogs(float a, float b);

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
    float newscore = clist[k].MEC;

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
                newscore -= Flist[f].currscore;
                update_fragscore(Flist, f, HAP1);
                newscore += Flist[f].currscore;
            }
        } else {
            clist[k].MEC = newscore;
            clist[k].bestMEC = newscore;
        }
    }
}

/**************** DETERMINISTIC MAX_CUT MEC IMPLEMENTATION *********************************************************//////

// function called from hapcut.c for each component...

int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, int iter) {
    int improved = 1; // assume improvement
    if (clist[k].iters_since_improvement > CONVERGE) {
        return 1; // return 1 only if converged
    }

    int i = 0, j = 0, t = 0;
    int f = 0;
    float cutvalue;
    float newscore;
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
    newscore = clist[k].MEC;
    // evaluate the impact of flipping of each SNP on MEC score, loop needed to improve MEC score
    single_variant_flips(Flist, snpfrag, clist, k, slist, HAP1);

    newscore = clist[k].MEC;
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
    if (cutvalue <= 3 || MINCUTALGO == 2) { //getchar();
        //printf("code reached here \n");
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
        if (clist[k].MEC > clist[k].bestMEC) // new haplotype is not better than current haplotype, revert to old haplotype
        {
            improved = 0;
            for (i = 0; i < clist[k].phased; i++) {
                if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
                else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
            }
            clist[k].MEC = 0;
            for (i = 0; i < clist[k].frags; i++) {
                update_fragscore(Flist, clist[k].flist[i], HAP1);
                clist[k].MEC += Flist[clist[k].flist[i]].currscore;
            }
        } else clist[k].bestMEC = clist[k].MEC; // update current haplotype
    }
    if (iter > 0 && clist[k].MEC > 0) fprintf(stdout, "component %d offset %d length %d phased %d  calls %d MEC %0.1f cutvalue %f bestMEC %0.2f\n", k, clist[k].offset, clist[k].length, clist[k].phased, clist[k].calls, clist[k].MEC, cutvalue, clist[k].bestMEC);

    if (improved) {
        clist[k].iters_since_improvement = 0;
    } else {
        clist[k].iters_since_improvement++;
    }
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
    float oldscore;
    int snp_add;
    int moved = 1, c1 = 0, c2 = 0;
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

// compare pruned snps by the log-likelihood of the dataset with them removed.
// this is so we can sort them and prune the worst ones.

int compare_hap_probs(const void *a, const void *b) {
    const struct hap_prob *ia = (const struct hap_prob*) a;
    const struct hap_prob *ib = (const struct hap_prob*) b;

    return ia->post_hap - ib->post_hap;
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

// returns an array length [snps] called 'pruned' that indicates which SNPs were pruned:
// 0 indicates not pruned
// 1 indicates pruned (leave unphased in output)
// 2 indicates called homozygous 00
// 3 indicates called homozygous 11
// the fraction to prune is determined by prune_threshold if THRESHOLD_TYPE is 1.
// else, the default behavior is that prune_threshold is the posterior probability cutoff.

void prune_snps(char* pruned, float prune_threshold, float homozygous_threshold, int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) {
    int i, j, f, num_to_prune;
    char temp1;
    float P_data_H, P_data_Hf, P_data_H00, P_data_H11, total;
    float log_hom_prior = log10(HOMOZYGOUS_PRIOR);
    float log_het_prior = log10(0.5 - HOMOZYGOUS_PRIOR);
    struct hap_prob* hap_probs = malloc(snps * sizeof (struct hap_prob));

    for (i = 0; i < snps; i++) {

        // this is for sorting, later
        hap_probs[i].snp_ix = i;

        // we only care about positions that are haplotyped
        if (!(HAP1[i] == '1' || HAP1[i] == '0')) {
            hap_probs[i].post_hap = 1;
            hap_probs[i].post_hapf = 0;
            hap_probs[i].post_00 = 0;
            hap_probs[i].post_11 = 0;
            continue;
        }

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
            P_data_H += simple_fragscore(Flist, f, HAP1, -1);
            // haplotypes with i flipped
            if (HAP1[i] == '1') HAP1[i] = '0';
            else if (HAP1[i] == '0') HAP1[i] = '1';
            P_data_Hf += simple_fragscore(Flist, f, HAP1, -1);
            // haplotypes with i homozygous 00
            HAP1[i] = '0';
            P_data_H00 += simple_fragscore(Flist, f, HAP1, i);
            // haplotypes with i homozygous 11
            HAP1[i] = '1';
            P_data_H11 += simple_fragscore(Flist, f, HAP1, i);
            //return haplotype to original value
            HAP1[i] = temp1;
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

    }

    // change the status of SNPs that are above the homozygous threshold
    // doing this in a separate loop just to be more decoupled
    for (i = 0; i < snps; i++) {
        // indexing by HAP[i] would do the exact same thing but this way it's robust to qsort() and therefore code structural changes.
        if (hap_probs[i].post_00 > log10(homozygous_threshold))
            pruned[hap_probs[i].snp_ix] = 2; // 2 specifies 00 homozygous
        else if (hap_probs[i].post_11 > log10(homozygous_threshold)) // 3 specifies 11 homozygous
            pruned[hap_probs[i].snp_ix] = 3;
    }

    if (THRESHOLD_TYPE == 1) {
        // prune threshold is a fraction of SNPs (useful mostly for fixing prune rate to compare tools)

        // sort haplotype probs so we can select the lowest values to prune
        qsort(hap_probs, snps, sizeof (struct hap_prob), compare_hap_probs);
        num_to_prune = (int) (prune_threshold * snps);
        // mark pruned SNPs as such in pruned array
        for (i = 0; i < num_to_prune; i++)
            pruned[hap_probs[i].snp_ix] = 1;
    } else {
        // prune threshold is a posterior probability (this is the default)
        for (i = 0; i < snps; i++) {
            if (hap_probs[i].post_hap < log10(prune_threshold))
                pruned[hap_probs[i].snp_ix] = 1;
        }
    }

    free(hap_probs);
}