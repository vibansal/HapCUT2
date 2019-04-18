#ifndef _HAPCONTIG_H
#define _HAPCONTIG_H
#include <stdint.h>
#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>

#include "common.h"

// haplotype block/connected component for phasing 

struct BLOCK {
    int offset;
    int length;
    int phased;
    char* haplotype;
    int* flist;
    int frags;
    float SCORE, bestSCORE, lastSCORE;
    int* slist; // ordered list of variants in this connected component
    int lastvar; // index of the first and last variants in this connected component
    int iters_since_improvement;
};


void generate_contigs(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int components, struct BLOCK* clist);

//void update_score(struct BLOCK* contig,struct fragment* Flist,char* HAP1);

int print_contigs(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* outfile);

float calculate_N50(struct BLOCK* clist, int blocks, struct SNPfrags* snpfrag, char* h1); // and other statistics, TBD

#endif
