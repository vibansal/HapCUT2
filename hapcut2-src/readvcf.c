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
    char buffer[10000];
    while (fgets(buffer, 10000, fp)) {
        if (buffer[0] == '#') {
            while (buffer[strlen(buffer) - 1] != '\n') fgets(buffer, 10000, fp);
            continue;
        }
        if (buffer[strlen(buffer) - 1] == '\n') variants++;
    }
    fclose(fp);
    fprintf_time(stderr, "read %d variants from %s file \n", variants, vcffile);
    return variants;
}

// 1000 genomes file have large deletions, need to allow longer VCF line
// this function is distinct from the function used to read VCF file for parsing haplotype informative reads, why ???
// read variants from VCF file, this code doesn't check the vcf file and assumes that column 10 contains the genotypes for the individual we are interested in phasing

int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps) {
    char buffer[500000];
    char temp[1000];
    char GQ[100];
    char* gen;
    int i = 0, j = 0, k=0, s = 0, e = 0, var = 0, GQ_ix, format_ix;
    FILE* fp = fopen(vcffile, "r");
    while (fgets(buffer, 500000, fp)) {
        if (buffer[0] == '#') continue;
        i = 0;
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        snpfrag[var].chromosome = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].chromosome[j - s] = buffer[j];
        snpfrag[var].chromosome[j - s] = '\0';
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        for (j = s; j < e; j++) temp[j - s] = buffer[j];
        temp[j - s] = '\0';
        snpfrag[var].position = atoi(temp);

        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        snpfrag[var].id = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].id[j - s] = buffer[j];
        snpfrag[var].id[j - s] = '\0';
        //strcpy(snpfrag[var].id,".");
        //variant->id = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->id[j-s] = buffer[j]; variant->id[j-s] = '\0';

        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        snpfrag[var].allele0 = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].allele0[j - s] = buffer[j];
        snpfrag[var].allele0[j - s] = '\0';


        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        snpfrag[var].allele1 = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].allele1[j - s] = buffer[j];
        snpfrag[var].allele1[j - s] = '\0';
        
	// set INDEL flag
	if (strlen(snpfrag[var].allele0) != 1 || strlen(snpfrag[var].allele1) != 1) snpfrag[var].is_indel = '1';
	else snpfrag[var].is_indel = '0'; 

        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        //variant->quality = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->quality[j-s] = buffer[j]; variant->quality[j-s] = '\0';
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        //variant->filter = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->filter[j-s] = buffer[j]; variant->filter[j-s] = '\0';
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;
        //variant->info = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->info[j-s] = buffer[j]; variant->info[j-s] = '\0';
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t') i++;
        e = i;

        // assert that GT field is the first field
        assert(buffer[s] == 'G' && buffer[s+1] == 'T' && (buffer[s+2] == ':' || buffer[s+2] == '\t'));

        // check format string for presence of GQ
        format_ix = 0;
        GQ_ix = -1; // the index of format field for GQ
        for (j=s;j<e;j++){
            if (buffer[j] == ':')
                format_ix++;
            else if((j+1 < e) && buffer[j] == 'G' && buffer[j+1] == 'Q'){
                GQ_ix = format_ix;
            }
        }

        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n') i++;
        e = i;

        snpfrag[var].genotypes = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].genotypes[j - s] = buffer[j];
        snpfrag[var].genotypes[j - s] = '\0';
        //len = j-s;
        gen = snpfrag[var].genotypes; //  for convenience
        snpfrag[var].homozygous_prior = HOMOZYGOUS_PRIOR; // default value
        int no_gq = 0;
        if (GQ_ix != -1) {

            // get to the index where GQ is located.
            k = 0; // where we are in gen
            for (j=0; j<GQ_ix; j++){
                if (no_gq){
                  break;
                }

                while(gen[k] != ':') {
                    k++;
                    if (gen[k] == '\n' || gen[k] == '\0') {
                        no_gq = 1;
                        break;
                    }
                }
                k++; // step past the ':'

                if (gen[k] == '\n' || gen[k] == '\0') {
                    no_gq = 1;
                    break;
                }
            }

            // reached GQ field. read it in.
            if (!no_gq) {
                j=0;

                while (gen[k] != '\0' && gen[k] != ':') {
                    GQ[j] = gen[k];
                    k++;
                    j++;
                }
                GQ[j] = '\0';

                snpfrag[var].homozygous_prior = atof(GQ) / -10.0; // log prior probability of homozygousity
            }
        }

        var++;
    }
    fclose(fp);
    return 1;
}

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
