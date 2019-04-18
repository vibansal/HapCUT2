#ifndef _VARIANTGRAPH_H
#define _VARIANTGRAPH_H
#include <stdint.h>
#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"

struct edge {
    int snp;
    int frag;
    char p[2];
    float w;
};

//void label_node(struct SNPfrags* snpfrag, int node, int comp);
float edge_weight(char* hap, int i, int j, char* p, struct fragment* Flist, int f);
void add_edges(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);
void add_edges_longreads(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);
void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps);
int edge_compare(const void *a, const void *b);


void print_variant(struct SNPfrags* snpfrag,int i,FILE* OUTFILE);

#endif
