
#ifndef INC_samread_H
#define INC_samread_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include "readfasta.h"
#include "sam.h"

extern int QVoffset;
extern int MINQ;
extern int MISSING_QV;

extern char INT_CIGAROP[];


struct alignedread
{
	int readlength;
	char* readid;
	char* chrom; char* matechrom; 
	int matech; char matestrand; 
	int flag; 	int position; int mquality; int mateposition; int IS;
	char* sequence; char* quality;
	char strand;
	int* cigarlist; int cigs;
	int mismatches; int indels; // no of mismatches and no of insertions/deletions
	int alignedbases; int clipped, gapped; 
	int cflag; int span;  
	int tid; int mtid; // matetid
	
	int findex; // index in array of fragments
	int mateindex; // index in array of reads of mate
	int blockid; int cluster; 

};

int fetch_func(const bam1_t *b, void *data,struct alignedread* read);

void free_readmemory(struct alignedread* read);

#endif
