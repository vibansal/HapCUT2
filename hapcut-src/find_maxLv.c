
#include <stdlib.h>
#include "fragmatrix.h"
int UNIMPROVED_CUTOFF = 10; // cut off iterations for a component if it hasn't improved log likelihood in this many iterations

void evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, int* slist, char* HAP1, char* HAP2, int iter);
float compute_goodcut(struct SNPfrags* snpfrag, char* HAP1, char* HAP2, char* slist, int N, struct fragment* Flist, struct BLOCK* clist, int k);

int evaluate_cut_component(struct fragment* Flist, char* S, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, char* HAP1, char* HAP2, int iter){
    int i, t;
    int N = clist[k].phased;
    
    // if this component hasn't improved in 10 iterations then skip evaluation
    if (clist[k].iters_since_improvement >= UNIMPROVED_CUTOFF){
        return 1; // unlikely to improve
    }
   
    int* slist = clist[k].slist
    
    // compute a 'max-cut' of columns
    clist[k].LL = compute_goodcut(snpfrag, HAP1, HAP2, S, slist, N, Flist, clist, k);
    
    if (clist[k].LL > clist[k].bestLL){
        // found a new best log likelihood
        clist[k].bestLL = clist[k].LL;
    }else{
        for (i = 0; i < N; i++){
            // haplotype wasn't improved
            // flip back haplotypes at members of S
            if (S[i]){
                fliphaps(HAP1, HAP2, i);

            }
        }
    }
    
    return 0; // component was improved
}

float compute_goodcut(struct SNPfrags* snpfrag, char* HAP1, char* HAP2, char* S, int* slist, int N, struct fragment* Flist, struct BLOCK* clist, int k){
    int i, j, max_i;
    int node1, node2, best_node1, best_node2;
    int min_i = -1;
    float edge_ll, best_ll, lv, max_lv;
    best_ll = -1;
    
    for (i = 0; i < N-1; i++) {
        node1 = slist[i]
        node2 = slist[i+1]
        S[node1] == 1;
        S[node2] == 1;
        edge_ll = abs(cut_difference_Lv(i, S, HAP1, HAP2, Flist, clist, k))
        S[node1] == 0;
        S[node2] == 0;
        if (edge_ll > best_ll){
            best_ll = edge_ll;
            best_node1 = node1;
            best_node2 = node2;
        }
    }
    
    S[best_node1] = 1;
    S[best_node2] = 1;

    fliphaps(HAP1, HAP2, best_node1);
    max_i = -1;
    max_lv = -1;
    for (i = 0; i < N; i++){
        if (S[slist[i]]){
            continue;
        }
        lv = cut_difference_Lv(i, S, HAP1, HAP2, Flist, clist, k);
        if(abs(lv) > max_lv){
            max_lv = lv;
            max_i = i;
        }
    }
    
    
    
    
    return 0;
}