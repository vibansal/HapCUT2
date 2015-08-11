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
    
#include <assert.h>
    
float fragment_likelihood_h(char* hap, struct fragment* Flist, int f, int partial, int* S);
float fragment_loglikelihood_H(char* hap1, char* hap2, struct fragment* Flist, int f, int partial, int* S);
float data_loglikelihood(char* hap1, char* hap2, struct fragment* Flist, int fragments, int partial, int* S);
float cut_difference_Lv(int v, int* S, char* hap1, char* hap2, struct fragment* Flist, int fragments);

#ifdef	__cplusplus
}
#endif

#endif	/* LIKELIHOOD_FUNCTIONS_H */

