#ifndef READVCF_H
#define READVCF_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"
//#include "htslib/hts.h"
//#include "htslib/vcf.h"


// counts # of variants in VCF file, copies from rarevariantassociation testing code (readvcf.c), allows for arbitrary long lines
int count_variants_vcf(char* vcffile);

int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps);

int count_variants(char* variantfile);

int read_variantfile(char* variantfile, struct SNPfrags* snpfrag, int snps);

#endif
