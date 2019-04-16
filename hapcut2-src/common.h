#ifndef _COMMON_H
#define _COMMON_H
#include <stdint.h>
#include "datastructures.h"

extern int QVoffset;
extern int MINQ;

#define MAXBUF 100000

// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 - pow(10, (b) - (a)))) : ((b) + log10(1.0 - pow(10.0, (a) - (b)))))

#define flip(allele) if (allele == '1') allele = '0'; else if (allele == '0') allele = '1'

int fprintf_time(FILE *stream, const char *format, ...);
float phred(float x);
float unphred(float x);
void check_input_0_or_1(char* x);

#endif
