
#ifndef _HAPFRAGMENT_H
#define _HAPFRAGMENT_H
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "readvariant.h"

extern int SINGLEREADS;
typedef struct
{
        char allele; char qv; int varid;  // allele is 0/1    varid is index to varlist[varid] gives all information about the variant
} allele;

typedef struct
{
        char* id;       int variants; allele* alist;
        int blocks; int paired; int matepos;
        
} FRAGMENT;

int compare_fragments(const void *a,const void *b);

int compare_alleles(const void *a,const void *b);

int print_fragment(FRAGMENT* fragment,VARIANT* varlist,FILE* outfile);

// make sure they are in the correct order, i+1 could be < i 
int print_matepair(FRAGMENT* f1, FRAGMENT* f2,VARIANT* varlist,FILE* outfile);

void clean_fragmentlist(FRAGMENT* flist,int* fragments,VARIANT* varlist,int currchrom,int currpos,int prevchrom);

#endif
