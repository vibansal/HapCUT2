#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H
#include <stdint.h>

// fragment blocks
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// haplotype block
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

struct edge {
    int snp;
    int frag;
    char p[2];
    float w;
};

typedef struct EDGE {
    int s, t;
    float w;
} EDGE;

struct SNPfrags { // single variant, also holds all information for node in graph 
    int* flist;
    int* jlist; // list of j indexes used to index into Flist[f].list
    int* klist; // list of k indexes used to index into Flist[f].list[j]
    int frags;
    char* alist; // alist is the alleles corresponding to each fragment that covers this SNP
    int component;
    int edges; // those that span the interval between this snp and the next snp
    int csize;
    struct edge* elist;
    int bcomp; // index of clist to which this snp belongs: reverse mapping
    struct edge* telist;
    int tedges; // temporary edge list and number of edges for MIN CUT computation
    int parent;
    float score;
    float htscore; // htrans
    int heaploc;
    char* id;
    char* allele0;
    char* allele1;
    char* chromosome;
    int position;
    // changed on feb 1 2012 to be pointers (char* id, char* chrom)
    char* genotypes; // VCF genotypes 0|1 1|0 or 0/1 added feb 1 2012
    char  ignore; // ignore this variant for all hapcut computations, 12/17/2018
    float post_notsw;
    float post_hap;
    int pruned_discrete_heuristic; // for error analysis mode
    float homozygous_prior; // prior probability of homozygousity. Based on GQ field of VCF.
    float PGLL[5]; // phased genotype likelihoods, 00,01,10,11 and ./. (if variant is ignored) 
    float hetLL; // unphased heterozygous 0/1 likelihood of genotype
};

#endif
