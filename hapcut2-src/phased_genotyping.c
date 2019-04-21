#include "common.h"
#include <float.h>
#include <assert.h>     /* assert */

const char* GTYPES[] = {"0/0","0/1","1/1","0|1","1|0"};

void init_varlist(PVAR* varlist,int snps,struct fragment* Flist,int fragments) // also calculates unphased genotype likelihoods (GLL)
{
    int i=0,j=0,k=0,s=0,h=0;
    float prob=0;
    float log_half = log10(0.5);
    for (i=0;i<snps;i++) 
    {
	varlist[i].frags =0;
	varlist[i].PGLL = calloc(5,sizeof(float));
	varlist[i].GLL = calloc(3,sizeof(float));
	for (j=0;j<3;j++) varlist[i].GLL[j] = 0.0; 
	varlist[i].maxhet = 0.0; // what is the best genotype likelihood for '0/1' genotype, useful for filtering
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
			varlist[s].maxhet += Flist[i].list[j].p1[k];
			prob = QVoffset - (int) Flist[i].list[j].qv[k]; prob /= 10;
			varlist[s].GLL[2] += prob;
			varlist[s].AC0 +=1;

		}
		if (Flist[i].list[j].hap[k] == '1') 
		{
			varlist[s].GLL[2] += Flist[i].list[j].p1[k];
			varlist[s].maxhet += Flist[i].list[j].p1[k];
			prob = QVoffset - (int) Flist[i].list[j].qv[k]; prob /= 10;
			varlist[s].GLL[0] += prob;
			varlist[s].AC1 +=1;
		}
	    	varlist[s].GLL[1] += log_half; // genotype likelihood where variant is 'unphased' heterozygous
	    }
        }
    }
}

void init_genotypes(DATA* data) // set the genotypes based on VCF file
{
    int i=0;
    for (i = 0; i < data->snps; i++) {
	if (data->HAP1[i] == '-') 
	{
		data->varlist[i].phased  = '0'; 
		if (data->snpfrag[i].genotypes[0] == '0' && data->snpfrag[i].genotypes[2] == '0') data->varlist[i].genotype = 0;
		else if (data->snpfrag[i].genotypes[0] == '1' && data->snpfrag[i].genotypes[2] == '1') data->varlist[i].genotype = 2;
		else if (data->snpfrag[i].genotypes[0] == '0' && data->snpfrag[i].genotypes[2] == '1') data->varlist[i].genotype = 1;
		else if (data->snpfrag[i].genotypes[2] == '0' && data->snpfrag[i].genotypes[0] == '1') data->varlist[i].genotype = 1;
	}
	else 
	{
		data->varlist[i].phased = '1'; 
		if (data->HAP1[i] == '0') data->varlist[i].genotype = 3; 
		else if (data->HAP1[i] == '1') data->varlist[i].genotype = 4; 
	}
    }
}

void free_varlist(PVAR* varlist,int snps)
{
    int i=0;
    for (i=0;i<snps;i++)
    {
        free(varlist[i].PGLL) ;    free(varlist[i].GLL);  free(varlist[i].flist);
    }
    free(varlist);
}

// find most likely genotype and calculate posterior prob. for a variant
int best_genotype(float PGLL[],float priors[],float* postp)
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

void calculate_geno_likelihoods(PVAR* varlist,int i,char* HAP1,struct fragment* Flist)
{
    float log_half = log10(0.5);
    int j=0,f=0; 
    char temp;
	
    for (j=0;j<5;j++) varlist[i].PGLL[j] = 0;
    temp = HAP1[i];
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
    HAP1[i] = temp;
}

void print_variant_GLL(struct SNPfrags* snpfrag,PVAR* varlist,char* HAP1,int i)
{
	fprintf(stdout,"VAR %d %9d %s %s (%.8s) %s ",i,snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1,snpfrag[i].genotypes,snpfrag[i].id);
	fprintf(stdout,"%s ",GTYPES[varlist[i].genotype]);
	fprintf(stdout," %c %3d %3.2f ",HAP1[i],varlist[i].frags,varlist[i].PGLL[1]);
	fprintf(stdout,"het: %3.2f %3.2f ",varlist[i].PGLL[3],varlist[i].PGLL[4]);
	fprintf(stdout,"0/0:%3.2f 1/1:%3.2f ",varlist[i].PGLL[0],varlist[i].PGLL[2]);
	fprintf(stdout,"%0.1f,%0.1f,%0.1f %0.1f %d:%d\n",varlist[i].GLL[0],varlist[i].GLL[1],varlist[i].GLL[2],varlist[i].AC0,varlist[i].AC1);
}

// what about filtering low-quality allele calls, this is not being done currently, use MINQ
// this function uses the full fragment list with alleles for homozygous variants 
int local_optimization(DATA* data) 
{
    float max,delta;
    int i, j, f;
    float hetprior = 0.001; 
    float priors[5] = {0.0,0.0,0.0,0.0,0.0};
    // prior for unphased variant should be lower = frequency of such variants...
    priors[1] = hetprior; priors[2] = hetprior/2; priors[0] = 1.0-priors[1]-priors[2]; 
    priors[3] = hetprior/2; priors[4] = hetprior/2;
    for (i=0;i<5;i++) priors[i]= log10(priors[i]);

    PVAR* varlist = data->varlist;
    init_varlist(varlist,data->snps,data->full_Flist,data->full_fragments); 
    init_genotypes(data);

    int changes=0; int iter=0; int flag=0;

    for (iter=0;iter < 4;iter++) {
	changes=0;

    for (i = 0; i < data->snps; i++) {
	if (varlist[i].frags < 1) continue; 
	calculate_geno_likelihoods(varlist,i,data->HAP1,data->full_Flist);
	varlist[i].bestgeno = best_genotype(varlist[i].PGLL,priors,&varlist[i].postp);
	delta = varlist[i].PGLL[varlist[i].bestgeno]+priors[varlist[i].bestgeno]-varlist[i].PGLL[varlist[i].genotype]-priors[varlist[i].genotype];

	flag =0;
	if (varlist[i].bestgeno != varlist[i].genotype && delta > 0) flag =1;
	if (iter ==0 || flag ==1) print_variant_GLL(data->snpfrag,data->varlist,data->HAP1,i);
	if (flag ==0) continue;

	if(varlist[i].bestgeno != varlist[i].genotype) fprintf(stdout,"change %s->%s %f ",GTYPES[varlist[i].genotype],GTYPES[varlist[i].bestgeno],varlist[i].postp);
	if (varlist[i].phased == '1')
	{
		if ( delta > 0 && (varlist[i].bestgeno ==0 || varlist[i].bestgeno == 1 || varlist[i].bestgeno ==2)) // phased to unphased
		{
			varlist[i].genotype = varlist[i].bestgeno;
			data->HAP1[i] = '-'; 
			varlist[i].phased = '0'; 
			fprintf(stdout,"phased2unphase");
			changes +=1;
		}
		else // phased to phased, 0|1 to 1|0 or vice versa
		{
			fprintf(stdout,"flipping-phase to %s",GTYPES[varlist[i].bestgeno]);
			varlist[i].genotype = varlist[i].bestgeno;
			if (data->HAP1[i] == '0') data->HAP1[i] = '1'; 
			else if (data->HAP1[i] == '1') data->HAP1[i] = '0'; 
		}
	}
	else // unphased variant 
	{
		if ( (varlist[i].bestgeno ==0 || varlist[i].bestgeno == 1 || varlist[i].bestgeno ==2)) // unphased to unphased
		{
			varlist[i].genotype = varlist[i].bestgeno;
			fprintf(stdout,"changing unphased genotype");
		}
		else if (varlist[i].postp >= 0.98) // unphased to phased, phased variants should be high confidence, >=0.9x PP
		{
			varlist[i].genotype = varlist[i].bestgeno;
			if (varlist[i].bestgeno == 3) data->HAP1[i] = '0';
			else if (varlist[i].bestgeno == 4) data->HAP1[i] = '1';
			varlist[i].phased = '1';
			fprintf(stdout,"unphased2phased");
			changes +=1;
		}
	}
	fprintf(stdout," iter:%d\n \n",iter+1);
    }
    fprintf(stderr,"phased genotype likelihood updates for data %d changes %d\n",data->snps,changes);
    // should run one round of HapCUT2 to update phased variant haplotype
    if (changes ==0) break;
    }
}
