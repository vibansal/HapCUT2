#ifndef READINPUTFILES_H
#define READINPUTFILES_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "simple.h"
#include "common.h"

int read_fragment_matrix(char* fragmentfile, struct fragment* Flist, int fragments);
int count_variants_vcf(char* vcffile);
int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps);
int fragment_compare(const void *a, const void *b);

#endif
