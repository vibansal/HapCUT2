#ifndef _FRAGMENTS_H
#define _FRAGMENTS_H
#include <stdint.h>
#include "variant.h"

// fragment block =  consecutive variants with allele-calls, each fragment is a linked list of such blocks
struct block {
    int offset;
    char* hap;
    short len;
    float* pv;
    char* qv;
    float* p1;
};

struct fragment {
    char* id;
    short blocks;
    struct block* list;
    int component;
    float currscore;
    int calls;
    float ll;
    float scores[4]; // added 03/02/15
    float htscores[4]; // scores assuming a hi-c h-trans interaction added 3/6/16
    int data_type; // data type -- 0:normal, 1:HiC
    float htrans_prob; // probability of an h-trans interaction for this read
    int mate2_ix;     // snp index of second mate; -1 if this fragment has one mate
    int isize;        // approximate insert size
    
    int PS; int PQ; char HP; // haplotype assignments for each fragment, HP= 0/1, PS = integer, PQ = probability that assignment of fragment is correct
};

int fragment_compare(const void *a, const void *b); // sorting fragment list 

void fragment_assignments(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag,char* h,char* outfile);

void calculate_fragscore1(struct fragment* Flist, int f, char* h, float* ll);
void update_fragscore1(struct fragment* Flist, int f, char* h);
float fragment_ll1(struct fragment* Flist, int f, char* h, int homozygous, int switch_ix);

void print_fragment(struct fragment* FRAG,FILE* OUTFILE);

int filter_fragments(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,struct fragment* nFlist);

void free_fragmentlist(struct fragment* Flist, int fragments); // freememory


//void filter_fragmentlist(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag);

#endif
