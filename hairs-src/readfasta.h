#ifndef INC_readfasta_H
#define INC_readfasta_H
#include <stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<ctype.h>

typedef struct {
    int chrom;
    int start;
    int end;
    char* annotation;
    // chrom indexes into REFLIST names
} INTERVAL;
// for storing list of intervals read from bedfile

// structure for chromosome/sequences/contigs 

typedef struct {
    int ns; // no of sequences 
    char** names;
    int* lengths;
    unsigned char** sequences; // changed to unsigned char july 5 2012 to avoid warnings in kmertable.c
    uint64_t* offsets; // from fasta index file for each chromosome

    int current;
    int* lookup; // current variable added for indexing into REFLIST// july 20 2011
    char fastafile[1024]; // name of fasta file

    INTERVAL* intervallist;
    int intervals;
    int* first_interval_chrom; // index to first interval for each chromosome in interval list
    int cinterval; // current interval that is closest to the current base being examined for variant calling

} REFLIST;


int read_fastaheader(char* fastafile, REFLIST* reflist);
int read_fasta(char* seqfile, REFLIST* reflist);
int read_chromosome(REFLIST* reflist, int chrom, FILE* fp);
int read_next_chromosome(REFLIST* reflist, int chrom, FILE* fp);
int read_bedfile(char* bedfile, REFLIST* reflist);
REFLIST* init_reflist(char* fastafile, REFLIST* reflist); // initialize reflist 

#endif
