#ifndef READINPUTFILES_H
#define READINPUTFILES_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"
#include "fragmatrix.h"

int read_fragment_matrix(char* fragmentfile, struct fragment* Flist, int fragments);

// counts # of variants in VCF file, copies from rarevariantassociation testing code (readvcf.c), allows for arbitrary long lines
int count_variants_vcf(char* vcffile);

int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps);

int count_variants(char* variantfile);

int read_variantfile(char* variantfile, struct SNPfrags* snpfrag, int snps);

int read_haplotypefile(char* hapfile, struct SNPfrags* snpfrag, int snps, char* HAP1, char* initHAP, int* bn);

int count_htrans_bins(char* htrans_file);

int read_htrans_file(char* htrans_file, float* htrans_probs, int num_bins);

#endif
