#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

/*
for Hi-C data, calculate number of unique linkages per variant
also calculate the mean and median base quality value

*/


int detect_long_reads(struct fragment* Flist,int fragments)
{
    int i=0;
    int long_reads=0;
    float mean_snps_per_read = 0;
    if (AUTODETECT_LONGREADS){
        for (i = 0; i < fragments; i++){
            mean_snps_per_read += Flist[i].calls;
        }
        mean_snps_per_read /= fragments;
        if (mean_snps_per_read >= 3){
            long_reads = 1;
        }else{
            long_reads = 0;
        }
    }
    fprintf(stderr,"mean number of variants per read is %0.2f \n",mean_snps_per_read);
    return long_reads;
}

void mean_insert_lengths()
{
}

void mean_fragment_lengths(struct Fragment* Flist,int fragments)
{
	// for linked-read data 
}
