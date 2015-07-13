#ifndef _FRAGMATRIX_H
#define _FRAGMATRIX_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"

//#define _GNU_SOURCE

int fraglength(struct fragment* Flist, int f);
float edge_weight(char* hap, int i, int j, char* p, struct fragment* Flist, int f);
int sample_block(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int start, int end, char* h1, int snps);

int compare_haps(struct BLOCK* clist, int components, char* orig, char* h1, struct SNPfrags* snpfrag, int snps);
int compare_Flist_hap(struct SNPfrags* snpfrag, int snps, struct fragment* Flist, int fragments, char* h, int Z, int QV);
int mutate_Flist(struct fragment* Flist, int fragments, double errrate);
int correct_fragment(struct fragment* Flist, int f, char* h);
int mecscore(struct fragment* Flist, int fragments, char* h, float* ll, float* calls, float* miscalls);
float compute_fragscore(struct fragment* Flist, int f, char* h, float* ll);
void update_fragscore(struct fragment* Flist, int f, char* h);
void frag_cluster_initialize(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* h1, int snps, struct BLOCK* clist, int comps);

void label_node(struct SNPfrags* snpfrag, int node, int comp);
void add_edges(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);
void add_edges_fosmids(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);

void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps);
void output_current_solution(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, char* hap, char* best);

int print_hapfile(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fname, int score, char* outfile);
void print_haplotypes_vcf(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, char* outfile);
void print_hapcut_options();

void print_fragmentmatrix_MEC(struct fragment* Flist, int fragments, char* h, char* outfileprefix);


int print_block_frags(struct BLOCK* blist, int block, char* aaron, char* h1, int* bn, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fragfile);

int print_blocks(struct BLOCK* blist, int blocks, char* aaron, char* h1, char* current, int* bn, struct fragment* Flist, int fragments, int mecscore, struct SNPfrags* snpfrag);

int print_block(struct BLOCK* blist, int block, char* aaron, char* h1, char* current, int* bn, struct fragment* Flist, int fragments, int mecscore, struct SNPfrags* snpfrag);

void generate_example_2(int n);

void generate_example(int n, int M);

int phase_score(char* hap, int i, int j, char* p);

int edge_compare(const void *a, const void *b);
int fragment_compare(const void *a, const void *b);
//int read_fragment_matrix(char* fragmentfile,struct fragment* Flist,int fragments );
//int read_variantfile(char* variantfile,struct SNPfrags* snpfrag,int snps);
//int read_haplotypefile(char* hapfile,struct SNPfrags* snpfrag,int snps,char* HAP1,char* initHAP,int* bn);

int determine_connected_components(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps);
void generate_clist_structure(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int components, struct BLOCK* clist);


#endif






