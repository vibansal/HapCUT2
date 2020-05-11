#include "readinputfiles.h"
#include "common.h"
#include <assert.h>

extern float HOMOZYGOUS_PRIOR;

// counts # of variants in VCF file, allows for arbitrary long lines
int count_variants_vcf(char* vcffile) {
    FILE* fp = fopen(vcffile, "r");
    if (fp == NULL) {
        fprintf_time(stderr, "could not open file %s \n", vcffile);
        return -1;
    }
    int variants = 0;
    int sbuf_size = 100000;
    char buffer[sbuf_size];
    while (fgets(buffer, sbuf_size, fp)) {

        if (buffer[0] == '#') { // only for header lines, keep reading until we get an endline character
            while (buffer[strlen(buffer) - 1] != '\n') 
	    {
		if (fgets(buffer, sbuf_size, fp) == NULL) break;
	    }
            continue;
        }
        if (buffer[strlen(buffer) - 1] == '\n') variants++; // count only when we get an endline character
    }
    fclose(fp);
    fprintf_time(stderr, "read %d variants from %s file \n", variants, vcffile);
    return variants;
}

// read variants from VCF file, this code  assumes that column 10 contains the genotypes for the individual we are interested in phasing

int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps) {
    char* buffer = (char*)malloc(1000000);
    int i = 0, k=0,k1=0,k2=0, var = 0;
    int GQ_ix, GT_ix;
    int l0=0,l1=0;
    FILE* fp = fopen(vcffile, "r");
    char** string_list = (char**)malloc(sizeof(char*)*4096); // if VCF file has more than 4096 columns... 
    char** format_list = (char**)malloc(sizeof(char*)*4096);
    char** geno_list = (char**)malloc(sizeof(char*)*4096);
    //char* format; 
    int header_lines =0;
    while (fgets(buffer, 1000000, fp)) {
        if (buffer[0] == '#') 
	{
		header_lines++;
		continue;
	}
	if (header_lines ==0 && var==0) fprintf(stderr,"VCF file has no header\n");

	k = splitString(buffer,'\t',string_list); 
	if (k < 10) 
	{
		fprintf(stderr,"ERROR reading VCF file: less than 10 columns in VCF file %s\n",vcffile); 
		return -1;
	}
	if (k > 10 && var ==0) 
	{
		fprintf(stderr,"more than 10 columns in VCF file %s:%d, HapCUT2 will use the genotypes in column 10 to phase the sample\n",vcffile,k); 
	}
	l0 = strlen(string_list[3]);  l1 = strlen(string_list[4]);  // length of alleles
	snpfrag[var].chromosome = (char*)malloc(strlen(string_list[0])+1); strcpy(snpfrag[var].chromosome,string_list[0]); 
	snpfrag[var].position = atoi(string_list[1]); 
	snpfrag[var].id = (char*)malloc(strlen(string_list[2])+1); strcpy(snpfrag[var].id,string_list[2]); 
	snpfrag[var].allele0 = (char*)malloc(l0+1); strcpy(snpfrag[var].allele0,string_list[3]); 
	snpfrag[var].allele1 = (char*)malloc(l1+1); strcpy(snpfrag[var].allele1,string_list[4]); 
	snpfrag[var].genotypes = (char*)malloc(strlen(string_list[9])+1); strcpy(snpfrag[var].genotypes,string_list[9]); 
        snpfrag[var].homozygous_prior = HOMOZYGOUS_PRIOR; // default value
	if (l0 != 1 || l1 !=1) snpfrag[var].is_indel = '1';
	else snpfrag[var].is_indel = '0'; 
	
        GT_ix = -1; GQ_ix = -1; // the index of format field for GQ
	k1 = splitString(string_list[8],':',format_list); 
	for (i=0;i<k1;i++)
	{
		if (strcmp(format_list[i],"GT") ==0) GT_ix = i;
		else if (strcmp(format_list[i],"GQ") ==0) GQ_ix = i; 
	}
	k2 = splitString(string_list[9],':',geno_list); 
	if (GQ_ix >= 0)
	{
        	snpfrag[var].homozygous_prior = atof(geno_list[GQ_ix]) / -10.0; // log prior probability of homozygousity if GQ is present
	}
	if (GT_ix == -1 || k1 != k2)  // k1 should be equal to k2
	{
		fprintf(stderr,"ERROR reading VCF file: no GT field found or mismatch between FORMAT and GQ cols %s\n",vcffile); 
		return -1;
	}
	//for (i=0;i<k;i++) fprintf(stdout,"%d %s\t",i,string_list[i]); fprintf(stdout,"\n");
	//for (i=0;i<k1;i++) fprintf(stdout,"%d %s\t",i,format_list[i]); fprintf(stdout,"\n\n");

	for (i=0;i<k1;i++) free(format_list[i]); // free array of split strings
	for (i=0;i<k2;i++) free(geno_list[i]); // free array of split strings
	for (i=0;i<k;i++) free(string_list[i]); // free array of split strings
	var++;
    }
    if (var != snps) fprintf(stderr,"ERROR, variant count does not equal to expected %d %d \n",var,snps); 
    fclose(fp);
    free(string_list); free(format_list); free(geno_list);
    free(buffer);
    return 1;
}


/////////////////////////////// OLD functions /////////////////////////////////////////////////////////////////////////////

int count_variants(char* variantfile) {
    int snps = 0;
    char buffer[MAXBUF];
    FILE* ff = fopen(variantfile, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open variant file %s\n", variantfile);
        exit(0);
    }
    while (fgets(buffer, MAXBUF, ff) != NULL) snps++;
    fclose(ff);
    return snps;
}
