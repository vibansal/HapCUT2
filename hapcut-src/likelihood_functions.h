/* 
 * File:   likelihood_functions.h
 * Author: peter
 *
 * Created on August 7, 2015, 5:08 PM
 */

#ifndef LIKELIHOOD_FUNCTIONS_H
#define	LIKELIHOOD_FUNCTIONS_H

#ifdef	__cplusplus
extern "C" {
#endif
 
#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include "common.h"
#include <assert.h>
#include <string.h>
    
int evaluate_cut_component(struct fragment* Flist, struct SNPfrags* snpfrag, struct BLOCK* clist, int k, char* HAP1, char* HAP2, int iter);
float compute_goodcut(struct SNPfrags* snpfrag, char* HAP1, char* HAP2,char* S1, char* S, int* slist, int N, struct fragment* Flist, struct BLOCK* clist, int k);

void flip_haps(char* hap1, char* hap2, int v);
float fragment_likelihood_h(char* hap, struct fragment* Flist, int f, char* S1, int partial, char* S);
float fragment_loglikelihood_H(char* hap1, char* hap2, struct fragment* Flist, int f, char* S1, int partial, char* S);
float data_loglikelihood(char* hap1, char* hap2, struct fragment* Flist, char* S1, int partial, char* S, int* frag_ix, int numfrags);
float cut_difference_Lv(int v, char* S1, char* S, char* hap1, char* hap2, struct fragment* Flist, struct BLOCK* clist, int k);

#ifdef	__cplusplus
}
#endif

#endif	/* LIKELIHOOD_FUNCTIONS_H */

