#ifndef HIC_H
#define HIC_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"

// HiC-related global variables
//extern int HIC;
//extern int MAX_HIC_EM_ITER;
extern int HTRANS_BINSIZE;
extern int HTRANS_MAXBINS; // this value will be overwritten at startup
extern int HTRANS_READ_LOWBOUND;
extern int HTRANS_MAX_WINDOW; // maximum window size for h-trans estimation

int count_htrans_bins(char* htrans_file);

int read_htrans_file(char* htrans_file, float* htrans_probs, int num_bins);

void init_HiC(struct fragment* Flist,int fragments,char* htrans_data_file);

// calculate new estimates of htrans error rates and write them to file 
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag,char* htrans_OUTFILE);

#endif
