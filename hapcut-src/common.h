
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

// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) ((a > b) ? (a + log10(1 + pow(10, b - a))) : (b + log10(1 + pow(10, a - b))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) ((a > b) ? (a + log10(1 - pow(10, b - a))) : (b + log10(1 - pow(10, a - b))))

#define flip(allele) if (allele == '1') allele = '0'; else if (allele == '0') allele = '1'

struct block {
    int offset;
    char* hap;
    short len;
    float* pv;
    char* qv;
    float* post;
    float* p1;
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
    int vbits;
    char sb; // single base fragment not useful for phasing 
    float scores[4]; // added 03/02/15
    float htscores[4]; // scores assuming a hi-c h-trans interaction added 3/6/16
    char init;
    int data_type; // data type -- 0:normal, 1:HiC
    float htrans_prob; // probability of an h-trans interaction for this read
    int mate2_ix;     // snp index of second mate; -1 if this fragment has one mate
    int isize;        // approximate insert size
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
    int split;
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
    int* jlist; // list of j indexes used to index into Flist[f].list
    int* klist; // list of k indexes used to index into Flist[f].list[j]
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
    float htscore; // htrans
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
    int prune_status;
    float post_notsw;
    float post_hap;
    int pruned_refhap_heuristic; // for error analysis mode
    char split;

    // added on april 24 2012 for singleton reads
    int A0, A1;
    float G00, G01, G11; // allele counts and genotype likelihoods using singleton reads that cover only one variant
    int rank; // rank of tree at this node
    float homozygous_prior; // prior probability of homozygousity. Based on GQ field of VCF.
};


#endif
