
#ifndef _COMMON_H
#define _COMMON_H
#include <stdint.h>
//Tue May 29 23:13:29 PDT 2007
extern int QVoffset;
extern int VCFformat;
extern int MINQ;
extern int FOSMIDS;
extern int SCORING_FUNCTION;
#define MAXBUF 10000

//fragment block

struct block {
    int offset;
    char* hap;    
    short len;
    float* pv;
    char* qv;
    float* post;
};

struct fragment {
    char* id;
    short blocks;
    struct block* list;
    char clust; // l indicates if it belongs to the haplotype or it's complement
    int component;
    float currscore;
    int calls;
    float ll;
    char sb; // single base fragment not useful for phasing 
};

// haplotype block

struct BLOCK {
    int offset;
    int length;
    int phased;
    char* haplotype;
    int* flist;
    int frags;
    //struct NODE* tree; int nodes; // nodes is the number of nodes in tree  
    float MEC, bestMEC, lastMEC;
    float LL, bestLL; // log likellihood scores 
    int calls;
    int* slist; // ordered list of variants in this connected component
    int lastvar; // index of the first and last variants in this connected component
    // firstvariant is same as offset
};

struct edge {
    //uint32_t snp; uint32_t frag; char p[2]; uint32_t w; 
    int snp;
    int frag;
    char p[2];
    float w;
};

typedef struct EDGE {
    int s, t;
    float w;
} EDGE;

struct SNPfrags {
    int* flist;
    int frags;
    char* alist; // alist is the alleles corresponding to each fragment that covers this SNP
    int component;
    int edges; // those that span the interval between this snp and the next snp 
    int csize;
    int blockno;
    short best_mec;
    struct edge* elist;
    int ff;
    int bcomp; // index of clist to which this snp belongs: reverse mapping  
    struct edge* telist;
    int tedges; // temporary edge list and number of edges for MIN CUT computation 
    int parent;
    float score;
    int revmap;
    int heaploc;
    char* id;
    char* allele0;
    char* allele1;
    char* chromosome;
    int position;
    // changed on feb 1 2012 to be pointers (char* id, char* chrom)
    char* genotypes; // VCF genotypes 0|1 1|0 or 0/1 added feb 1 2012 
    float L00, L01, L10, L11, Lnovar; // change in likelihood if this SNP is made homozygous or removed
    float rMEC;
    int R0, R1; // counts of bases supporting allele0 and allele1

    // added on april 24 2012 for singleton reads
    int A0, A1;
    float G00, G01, G11; // allele counts and genotype likelihoods using singleton reads that cover only one variant
    int rank; // rank of tree at this node
};

#endif
