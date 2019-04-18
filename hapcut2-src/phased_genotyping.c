#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

extern float HOMOZYGOUS_PRIOR;

// not used currently
void calculate_hetLL(struct fragment* Flist, int fragments,struct SNPfrags* snpfrag, int snps)
{
	float log_half = log10(0.5);
	int i=0,j=0,k=0,snp_ix=0,f=0;
    	for (i = 0; i < snps; i++) snpfrag[i].hetLL = 0.0;
	for (f=0;f<fragments;f++)
	{
		for (j = 0; j < Flist[f].blocks; j++) 
		{
		    for (k = 0; k < Flist[f].list[j].len; k++) 
		    {
			snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs
			snpfrag[snp_ix].hetLL += log_half;
		    }
		}
	}
}
// what about filtering low-quality allele calls, this is not being done currently
// this function needs the full fragment list, need to update snpfrag[i].frags and snpfrag[i].flist 
void local_optimization(int snps, struct fragment* Flist, struct SNPfrags* snpfrag, char* HAP1) 
{   
    float log_half = log10(0.5);
    float ll,max,delta;
    int i, j, f,unphased;
    char temp1;
    fprintf(stderr,"phased genotype likelihood calculation for data %d \n",snps);

    for (i = 0; i < snps; i++) {

	if (HAP1[i] == '-') print_variant(snpfrag,i,stdout);

	for (j=0;j<5;j++) snpfrag[i].PGLL[j] = 0;
	snpfrag[i].hetLL = 0.0; ll =0.0;

        temp1 = HAP1[i];
        for (j = 0; j < snpfrag[i].frags; j++) {
            f = snpfrag[i].flist[j];
            HAP1[i] = '0'; 
	    snpfrag[i].PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1); //  genotype = 0|0, variant 'i' is homozygous
            snpfrag[i].PGLL[1] += fragment_ll1(Flist, f, HAP1, -1, -1); // genotype= 0|1
            HAP1[i] = '1'; 
	    snpfrag[i].PGLL[3] += fragment_ll1(Flist,f, HAP1, i, -1); // genotype = 1|1, variant 'i' is homozygous
	    snpfrag[i].PGLL[2] += fragment_ll1(Flist,f, HAP1, -1, -1); // genotype = 1|0

	    ll = fragment_ll1(Flist,f, HAP1, i, 10); // variant 'i' ignored for likelihood calculation 
	    if (ll < 0) snpfrag[i].PGLL[4] +=ll + log_half;  // genotype likelihood where variant is 'unphased' heterozygous
	    else snpfrag[i].PGLL[4] += log_half;
	    snpfrag[i].hetLL += log_half; // simply x 0.5 for each allele 
	}
        //return haplotype to original value
        HAP1[i] = temp1;

	max = snpfrag[i].PGLL[1];  
	if (snpfrag[i].PGLL[2] > max) max = snpfrag[i].PGLL[2];
	delta = snpfrag[i].PGLL[4]-max;

	if (HAP1[i] == '0') temp1 = '1'; else temp1 = '0';
	fprintf(stdout,"SNP %d %d %s %s ",i,snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1);
	fprintf(stdout," %c|%c cov %d %0.2f ",HAP1[i],temp1,snpfrag[i].frags,snpfrag[i].PGLL[4]);
	fprintf(stdout,"hets= %0.2f %0.2f ",snpfrag[i].PGLL[1],snpfrag[i].PGLL[2]);
	fprintf(stdout,"homs= %0.2f %0.2f delta %0.2f ",snpfrag[i].PGLL[0],snpfrag[i].PGLL[3],delta);
	if (snpfrag[i].PGLL[0] > max) fprintf(stdout,"ref-hom %0.2f ",snpfrag[i].PGLL[0]-max);
	if (snpfrag[i].PGLL[3] > max) fprintf(stdout,"alt-hom %0.2f ",snpfrag[i].PGLL[3]-max);
	//if (P_data_Hf > P_data_H) fprintf(stdout,"flip-phase ");
	else if (delta > 0) fprintf(stdout,"unphased ");
	fprintf(stdout,"\n");

    }
}

