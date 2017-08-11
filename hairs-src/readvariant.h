#ifndef _READVARIANT_H
#define _READVARIANT_H
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "hashtable.h"
#include "readfasta.h"
//#define _GNU_SOURCE
extern FILE* fragment_file; // FILE to which the fragments will be output, if NULL, output to stdout

extern int TRI_ALLELIC;

extern int BSIZE;
extern int PRINT_FRAGMENTS;
//int VARIANTS = 0;

typedef struct {
    char* id;
    char* chrom;
    int position;
    short altalleles;
    char* RA;
    char* AA; // alternate alleles
    double* GLL; // genotype likelihoods added 11/25/13
    char* genotype; // encoded as integers 0 1 2 3 4 5 6 7
    short type;
    // changed this to char* on April 3 2012
    char* allele1;
    char* allele2; // temporary for SNPs
    char heterozygous; // only heterozygous variants will be used for printing out HAIRS
    int depth;
    int A1, A2;
    int H1, H2;
    // total reads covering this variant (haploid/diploid, A1-> reads supporting reference allele (single-read)
    //	float L11,L12,L22; // genotype likelihoods for three possible genotypes
} VARIANT;

// information about the variants on each chromosome

typedef struct {
    int variants;
    int first;
    int last;
    int blocks;
    int* intervalmap;
} CHROMVARS;

/*
typedef struct
{
        char allele; char qv; int varid;  // allele is 0/1    varid is index to varlist[varid] gives all information about the variant
} allele;

typedef struct
{
        char* id; 	int variants; allele* alist;
        int blocks; int paired; int matepos;

} FRAGMENT;

int compare_fragments(const void *a,const void *b);
 */

int count_variants(char* vcffile, char* sampleid, int* samplecol);

int count_variants_oldformat(char* snpfile);
int read_variantfile_oldformat(char* snpfile, VARIANT* varlist, HASHTABLE* ht, int snps);

int parse_variant(VARIANT* variant, char* buffer, int samplecol);

int read_variantfile(char* vcffile, VARIANT* varlist, HASHTABLE* ht, int* hetvariants, int samplecol);

void build_intervalmap(CHROMVARS* chromvars, int chromosomes, VARIANT* varlist, int variants);

int calculate_rightshift(VARIANT* varlist, int ss, REFLIST* reflist);

#endif
