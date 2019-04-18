#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

extern float HOMOZYGOUS_PRIOR;

// similar to SNPfrag but compact
typedef struct
{
	int frags;
	int* flist;
	float* PGLL; // phased genotype likelihoods
	float* GLL;  // unphased genotype likelihoods
	float* priors; // genotype priors
	
} PVAR;

void init_PVAR(PVAR* varlist,int snps,struct fragment* Flist,int fragments)
{
    int i=0,j=0,k=0,h=0,s=0;

    for (i=0;i<snps;i++) 
    {
	varlist[i].frags =0;
	varlist[i].PGLL = calloc(5,sizeof(float));
	varlist[i].GLL = calloc(3,sizeof(float));
	for (j=0;j<3;j++) varlist[i].GLL[j] = 0.0; 
    }
    // calculate number of fragments covering each variant 
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) varlist[Flist[i].list[j].offset + k].frags++;
        }
    }
    for (i = 0; i < snps; i++) {
        varlist[i].flist = (int*) malloc(varlist[i].frags*sizeof (int));
        varlist[i].frags = 0;
    }
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++)
            {
                s = Flist[i].list[j].offset + k;           // index in snp list
                h = varlist[s].frags; 
               	varlist[s].flist[h] = i;
                varlist[s].frags++;
		// we can update genotype likelihoods GLL
	    }
        }
    }
}

// what about filtering low-quality allele calls, this is not being done currently
// this function needs the full fragment list, need to update varlist[i].frags and varlist[i].flist 
void local_optimization(DATA* data) 
{
    int snps = data->snps;
    int fragments = data->full_fragments;
    struct fragment* Flist = data->full_Flist;
    char* HAP1 = data->HAP1;    
    PVAR* varlist = calloc(sizeof(PVAR),snps);
    init_PVAR(varlist,snps,Flist,fragments); 

    float log_half = log10(0.5);
    float ll,max,delta;
    int i, j, f,unphased;
    char temp1;
    fprintf(stderr,"phased genotype likelihood calculation for data %d \n",snps);

    for (i = 0; i < snps; i++) {
	//if (HAP1[i] == '-') print_variant(data->snpfrag,i,stdout);
	for (j=0;j<5;j++) varlist[i].PGLL[j] = 0;
	ll =0.0;

        temp1 = HAP1[i];
        for (j = 0; j < varlist[i].frags; j++) {
            f = varlist[i].flist[j];
            HAP1[i] = '0'; 
	    // fragment_ll is not divided by half
	    varlist[i].PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1) + log_half; //  genotype = 0|0, variant 'i' is homozygous
            varlist[i].PGLL[1] += fragment_ll1(Flist, f, HAP1, -1, -1)+ log_half; // genotype= 0|1
            HAP1[i] = '1'; 
	    varlist[i].PGLL[2] += fragment_ll1(Flist,f, HAP1, -1, -1) + log_half; // genotype = 1|0
	    varlist[i].PGLL[3] += fragment_ll1(Flist,f, HAP1, i, -1) + log_half; // genotype = 1|1, variant 'i' is homozygous

	    varlist[i].PGLL[4] += fragment_ll1(Flist,f, HAP1, i, 10) + 2*log_half; // variant 'i' ignored for likelihood calculation  + heterozygous
	    varlist[i].GLL[1] += log_half; // genotype likelihood where variant is 'unphased' heterozygous
	}
        //return haplotype to original value
        HAP1[i] = temp1;

	max = varlist[i].PGLL[1];  
	if (varlist[i].PGLL[2] > max) max = varlist[i].PGLL[2];
	delta = varlist[i].PGLL[4]-max;

	if (HAP1[i] == '0') temp1 = '1'; 
	else if (HAP1[i] == '1') temp1 = '0';
	else temp1 = '-';
	fprintf(stdout,"SNP %d %9d %s %s (%.3s) ",i,data->snpfrag[i].position,data->snpfrag[i].allele0,data->snpfrag[i].allele1,data->snpfrag[i].genotypes);
	fprintf(stdout," %c|%c cov %3d %3.2f ",HAP1[i],temp1,varlist[i].frags,varlist[i].PGLL[4]);
	fprintf(stdout,"het: %3.2f %3.2f ",varlist[i].PGLL[1],varlist[i].PGLL[2]);
	fprintf(stdout,"hom: %3.2f %3.2f delta %0.2f ",varlist[i].PGLL[0],varlist[i].PGLL[3],delta);
	if (varlist[i].PGLL[0] > max) fprintf(stdout,"ref-hom %0.2f ",varlist[i].PGLL[0]-max);
	if (varlist[i].PGLL[3] > max) fprintf(stdout,"alt-hom %0.2f ",varlist[i].PGLL[3]-max);
	//if (P_data_Hf > P_data_H) fprintf(stdout,"flip-phase ");
	else if (delta > 0) fprintf(stdout,"unphased ");
	fprintf(stdout,"\n");

    }
}

