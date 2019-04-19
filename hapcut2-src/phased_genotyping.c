#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

// similar to SNPfrag but only for genotype likelihoods
typedef struct
{
	int frags;
	int* flist;
	float* PGLL; // phased genotype likelihoods, 0/0, 0/1, 1/1, 0|1, 1|0 
	float* GLL;  // unphased genotype likelihoods
	float postp;
	char phased; // phased or unphased
	int genotype; // one of 5 possible genotypes, we can add 0/2, 1/2, 2/2 for three alleles if needed
	int bestgeno;
	
} PVAR;

const char* GTYPES[] = {"0/0","0/1","1/1","0|1","1|0"};


void init_PVAR(PVAR* varlist,int snps,struct fragment* Flist,int fragments) // also calculates unphased genotype likelihoods (GLL)
{
    int i=0,j=0,k=0,h=0,s=0;
    float prob=0;
    float log_half = log10(0.5);

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
		if (Flist[i].list[j].hap[k] == '0') 	
		{
			varlist[s].GLL[0] += Flist[i].list[j].p1[k];
			prob = QVoffset - (int) Flist[i].list[j].qv[k]; prob /= 10;
			varlist[s].GLL[2] += prob;

		}
		if (Flist[i].list[j].hap[k] == '1') 
		{
			varlist[s].GLL[2] += Flist[i].list[j].p1[k];
			prob = QVoffset - (int) Flist[i].list[j].qv[k]; prob /= 10;
			varlist[s].GLL[0] += prob;
		}
	    	varlist[s].GLL[1] += log_half; // genotype likelihood where variant is 'unphased' heterozygous
	    }
        }
    }
}

void free_PVAR(PVAR* varlist,int snps)
{
    int i=0;
    for (i=0;i<snps;i++)
    {
        free(varlist[i].PGLL) ;    free(varlist[i].GLL);  free(varlist[i].flist);
    }
    free(varlist);
}

// find most likely genotype and calculate posterior prob. for a variant
int ML_genotype(float PGLL[],float priors[],float* postp)
{
	int i=0,best=0;
	float bestLL=-1000000;
	float sumprobs=0;
	for (i=0;i<5;i++)
	{
		if (PGLL[i]+priors[i] > bestLL) 
		{
			best = i; 
			bestLL = PGLL[i] + priors[i]; 
		}
		if (i==0) sumprobs = PGLL[i] + priors[i];
		else sumprobs = addlogs(PGLL[i]+priors[i],sumprobs);
	}
	*postp = pow(10,PGLL[best]+priors[best]-sumprobs);
	return best;
}

void calculate_PGLL(PVAR* varlist,int i,char* HAP1,struct fragment* Flist)
{
    float log_half = log10(0.5);
	int j=0,f=0; 
    char temp1;
	
	for (j=0;j<5;j++) varlist[i].PGLL[j] = 0;
        temp1 = HAP1[i];
        for (j = 0; j < varlist[i].frags; j++) {
            f = varlist[i].flist[j];
            HAP1[i] = '0'; 
	    // fragment_ll is not divided by half
	    varlist[i].PGLL[0] += fragment_ll1(Flist,f, HAP1, i, -1) + log_half; //  genotype = 0|0, variant 'i' is homozygous
            varlist[i].PGLL[3] += fragment_ll1(Flist, f, HAP1, -1, -1)+ log_half; // genotype= 0|1
            HAP1[i] = '1'; 
	    varlist[i].PGLL[4] += fragment_ll1(Flist,f, HAP1, -1, -1) + log_half; // genotype = 1|0
	    varlist[i].PGLL[2] += fragment_ll1(Flist,f, HAP1, i, -1) + log_half; // genotype = 1|1, variant 'i' is homozygous

	    varlist[i].PGLL[1] += fragment_ll1(Flist,f, HAP1, i, 10) + 2*log_half; // variant 'i' ignored for likelihood calculation  + heterozygous, 0/1
	}
        //return haplotype to original value
        HAP1[i] = temp1;
}

// identify 'SOLID' variants that don't need to be updated 
// what about filtering low-quality allele calls, this is not being done currently
// this function needs the full fragment list with alleles for homozygous variants 
void local_optimization(DATA* data) 
{
    float max,delta;
    int i, j, f;
    char temp1;

    float hetprior = 0.001; 
    float priors[5] = {0.0,0.0,0.0,0.0,0.0};
    priors[1] = hetprior; priors[2] = hetprior/2; priors[0] = 1.0-priors[1]-priors[2];
    priors[3] = hetprior/2; priors[4] = hetprior/2;
    for (i=0;i<5;i++) priors[i]= log10(priors[i]);

    char* HAP1 = data->HAP1;    
    PVAR* varlist = calloc(sizeof(PVAR),data->snps);
    init_PVAR(varlist,data->snps,data->full_Flist,data->full_fragments); 

    fprintf(stderr,"phased genotype likelihood calculation for data %d \n",data->snps);

    for (i = 0; i < data->snps; i++) {
	if (HAP1[i] == '-') 
	{
		varlist[i].phased  = '0'; 
		if (data->snpfrag[i].genotypes[0] == '0' && data->snpfrag[i].genotypes[2] == '0') varlist[i].genotype = 0;
		else if (data->snpfrag[i].genotypes[0] == '1' && data->snpfrag[i].genotypes[2] == '1') varlist[i].genotype = 2;
		else if (data->snpfrag[i].genotypes[0] == '0' && data->snpfrag[i].genotypes[2] == '1') varlist[i].genotype = 1;
		else if (data->snpfrag[i].genotypes[2] == '0' && data->snpfrag[i].genotypes[0] == '1') varlist[i].genotype = 1;
	}
	else 
	{
		varlist[i].phased = '1'; 
		if (HAP1[i] == '0') varlist[i].genotype = 3; 
		else if (HAP1[i] == '1') varlist[i].genotype = 4; 
	}
	calculate_PGLL(varlist,i,data->HAP1,data->full_Flist);

	varlist[i].bestgeno = ML_genotype(varlist[i].PGLL,priors,&varlist[i].postp);

	delta = varlist[i].PGLL[varlist[i].bestgeno]-varlist[i].PGLL[varlist[i].genotype];

	if (HAP1[i] == '0') temp1 = '1'; 
	else if (HAP1[i] == '1') temp1 = '0';
	else temp1 = '-';
	fprintf(stdout,"SNP %d %9d %s %s (%.3s) ",i,data->snpfrag[i].position,data->snpfrag[i].allele0,data->snpfrag[i].allele1,data->snpfrag[i].genotypes);
	fprintf(stdout," %c|%c %3d %3.2f ",HAP1[i],temp1,varlist[i].frags,varlist[i].PGLL[1]);
	fprintf(stdout,"het: %3.2f %3.2f ",varlist[i].PGLL[3],varlist[i].PGLL[4]);
	fprintf(stdout,"hom: %3.2f %3.2f ",varlist[i].PGLL[0],varlist[i].PGLL[2]);
	if(varlist[i].bestgeno != varlist[i].genotype) fprintf(stdout,"change %s->%s %f ",GTYPES[varlist[i].genotype],GTYPES[varlist[i].bestgeno],varlist[i].postp);

	if (delta >= 1 && varlist[i].phased == '1' && varlist[i].bestgeno ==1) fprintf(stdout,"unphase_genotype delta %0.2f",delta);
	else if (delta > 0 && varlist[i].bestgeno != varlist[i].genotype && varlist[i].bestgeno ==0) fprintf(stdout,"non-variant delta %0.2f",delta);
	else if (delta > 1 && varlist[i].bestgeno != varlist[i].genotype) fprintf(stdout,"update_genotype delta %0.2f",delta);
	fprintf(stdout,"\n");
    }
}

