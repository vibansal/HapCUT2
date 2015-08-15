#include "readinputfiles.h"

int read_fragment_matrix(char* fragmentfile, struct fragment* Flist, int fragments) {
    int i = 0, j = 0, k = 0, t = 0, t1 = 0;
    int blocks = 0, type = 0, l = 0, biter = 0, offset = 0;
    char buffer[MAXBUF];
    char blockseq[5000];
    char ch;

    FILE* ff = fopen(fragmentfile, "r");
    if (ff == NULL) {
        fprintf(stderr, "couldn't open fragment file \n");
        return -1;
    }
    for (i = 0; i < fragments; i++) {
        //		fprintf(stdout,"%s \n",buffer);
        j = 0;
        ch = fgetc(ff);
        while (ch != '\n') {
            buffer[j] = ch;
            j++;
            ch = fgetc(ff);
        }
        buffer[j] = '\0';
        k = 0;
        t = 0;
        type = 0;
        while (k < j) {
            while (buffer[k] != ' ' && buffer[k] != '\t' && k < j && buffer[k] != '\0') {
                blockseq[t] = buffer[k];
                t++;
                k++;
            }
            k++;
            while ((buffer[k] == ' ' || buffer[k] == '\t') && k < j) k++;
            blockseq[t] = '\0';

            if (type == 0) // read the number of blocks in fragment
            {
                blocks = 0;
                for (l = 0; l < t; l++) {
                    blocks = 10 * blocks + (int) (blockseq[l] - 48);
                }
                type = 1;
                Flist[i].blocks = blocks;
                Flist[i].list = (struct block*) malloc(sizeof (struct block)*(blocks));
                biter = 0;
                //printf("blocks %d \n",blocks);
            } else if (type == 1) // read the fragment id, changed to allow dynamic length feb202011
            {
                Flist[i].id = (char*) malloc(t + 1);
                for (l = 0; l < t; l++) Flist[i].id[l] = blockseq[l];
                Flist[i].id[l] = '\0';
                //Flist[i].id = (char*)malloc(1); Flist[i].id[0] = '0'; this doesnt reduce the memory requirement
                type = 2;
            } else if (type == 2 && biter < blocks) {
                offset = 0;
                for (l = 0; l < t; l++) {
                    offset = 10 * offset + (int) (blockseq[l] - 48);
                }
                type = 3;
                Flist[i].list[biter].offset = offset - 1;
                //printf("block %d %d ",biter,offset-1);
            } else if (type == 2 && biter == blocks) {
                offset = 0;
                Flist[i].calls = 0;
                for (l = 0; l < blocks; l++) {
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].pv[t1] = pow(0.1, (float) (blockseq[offset + t1] - QVoffset) / 10);
                    //for (t1=0;t1<Flist[i].list[l].len;t1++) printf("qv %f %d ",pow(0.1,(float)(blockseq[offset+t1]-QVoffset)/10),blockseq[offset+t1]-33);
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].qv[t1] = blockseq[offset + t1];
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].p1[t1] = log10(1.0 - Flist[i].list[l].pv[t1]); // added 03/03/15
                    offset += Flist[i].list[l].len;
                    Flist[i].calls += Flist[i].list[l].len;
                }
            } else if (type == 3) {
                Flist[i].list[biter].hap = (char*) malloc(t + 1);
                Flist[i].list[biter].qv = (char*) malloc(t + 1);
                Flist[i].list[biter].len = t;
                Flist[i].list[biter].pv = (float*) malloc(sizeof (float)*Flist[i].list[biter].len);
                Flist[i].list[biter].p1 = (float*) malloc(sizeof (float)*Flist[i].list[biter].len);

                for (l = 0; l < t; l++) Flist[i].list[biter].hap[l] = blockseq[l];
                //for (l=0;l<t;l++) Flist[i].list[biter].pv[l] = 0.005; for (l=0;l<t;l++) Flist[i].list[biter].qv[l] = 'A';

                //for (l=0;l<t;l+=2) Flist[i].list[biter].qv[l/2] = blockseq[l+1];
                //for (l=0;l<t;l+=2) Flist[i].list[biter].pv[l/2] = pow(0.1,(float)(blockseq[l+1]-QVoffset)/10); 
                //for (t1=0;t1<Flist[i].list[biter].len;t1++) Flist[i].list[biter].pv[t1] = 0.01;

                //Flist[i].list[biter].post = (float*)malloc(sizeof(float)*Flist[i].list[biter].len); // how many times it matches 
                //for (t1=0;t1<Flist[i].list[biter].len;t1++) Flist[i].list[biter].post[t1] =0;
                type = 2;
                biter++;
            }
            t = 0;
        }
    }
    fclose(ff);
    qsort(Flist, fragments, sizeof (struct fragment), fragment_compare);
    //for (i=0;i<fragments;i++)  fprintf(stdout,"fragment %d blocks %d offset %d\n",i,Flist[i].blocks,Flist[i].list[0].offset);
    /****************************** READ FRAGMENT QUALITY FILE*************************************************/
    //	ff = fopen(qualityfile,"r"); if (ff == NULL || QV == -1) fprintf(stderr,"couldn't open fragment QV file \n");
    return fragments;
}


// counts # of variants in VCF file, allows for arbitrary long lines

int count_variants_vcf(char* vcffile) {
    FILE* fp = fopen(vcffile, "r");
    if (fp == NULL) {
        fprintf(stderr, "could not open file %s \n", vcffile);
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
    fprintf(stderr, "read %d variants from %s file \n", variants, vcffile);
    return variants;
}

// 1000 genmes file have large deletions, need to allow longer VCF line
// this function is distinct from the function used to read VCF file for parsing haplotype informative reads, why ???
// read variants from VCF file, this code doesn't check the vcf file and assumes that column 10 contains the genotypes for the individual we are interested in phasing

int read_vcffile(char* vcffile, struct SNPfrags* snpfrag, int snps) {
    char buffer[100000];
    char temp[1000];
    int i = 0, j = 0, s = 0, e = 0, var = 0;
    FILE* fp = fopen(vcffile, "r");
    while (fgets(buffer, 100000, fp)) {
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
        //variant->format = (char*)malloc(e-s+1); for (j=s;j<e;j++) variant->format[j-s] = buffer[j]; variant->format[j-s] = '\0';

        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n') i++;
        e = i;
        snpfrag[var].genotypes = (char*) malloc(e - s + 1);
        for (j = s; j < e; j++) snpfrag[var].genotypes[j - s] = buffer[j];
        snpfrag[var].genotypes[j - s] = '\0';
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
        fprintf(stderr, "couldn't open variant file %s\n", variantfile);
        exit(0);
    }
    while (fgets(buffer, MAXBUF, ff) != NULL) snps++;
    fclose(ff);
    return snps;
}

int read_variantfile(char* variantfile, struct SNPfrags* snpfrag, int snps) {
    // more ../../../Haplotype-assembly/Dec8hairs/NA18507.variants.sorted.chr6 | awk 'BEGIN {i =1; } {split($6,g,"/"); if (g[1] != g[2]) {print $1"_"i,$2,$3,$4,$5,$6,$7; i++;} }' > NA18507.chr6.hetvariants.inputforhapcut
    // read initial variant file in format:  vartype varid chrom position allele1 allele2 genotype score 
    // the number of lines in this file should be === snps unless we want to allow homozygous variants in this file 
    // the samfile parser can also read this same file for extracting the fragment matrix 
    int i = 0;
    char varid[256];
    char chromosome[256];
    char allele0[256];
    char allele1[256];
    char genotype[64];
    int score;
    // need to allocate memory for variables id, chromosome, allele0, allele1 
    FILE* sf = fopen(variantfile, "r");
    for (i = 0; i < snps; i++) {
        fscanf(sf, "%s %s %d %s %s %s %d \n", varid, chromosome, &snpfrag[i].position, allele0, allele1, genotype, &score);
        snpfrag[i].id = (char*) malloc(strlen(varid));
        strcpy(snpfrag[i].id, varid);
        snpfrag[i].chromosome = (char*) malloc(strlen(chromosome));
        strcpy(snpfrag[i].chromosome, chromosome);
        snpfrag[i].allele0 = (char*) malloc(strlen(allele0));
        strcpy(snpfrag[i].allele0, allele0);
        snpfrag[i].allele1 = (char*) malloc(strlen(allele1));
        strcpy(snpfrag[i].allele1, allele1);
        snpfrag[i].genotypes = (char*) malloc(strlen(genotype));
        strcpy(snpfrag[i].genotypes, genotype);
    }
    //	for (i=0;i<snps;i++) fscanf(sf,"%s %s %d %s %s %s %d \n",snpfrag[i].id,snpfrag[i].chromosome,&snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1,genotype,&score);
    fprintf(stderr, "read variants from variantfile %s \n", variantfile);
    fclose(sf); // bug fixed feb 1 2012 
    return 1;
}

int read_haplotypefile(char* hapfile, struct SNPfrags* snpfrag, int snps, char* HAP1, char* initHAP, int* bn) {
    /****************************** READ HAPLOTYPE SOLUTION*************************************************/
    int i = 0, j = 0;
    char id[100];
    char c1, c2;
    int offset, len, phased, blocks;
    FILE* sf = fopen(hapfile, "r");
    struct BLOCK* blist;
    if (sf == NULL) fprintf(stderr, "couldn't open initial haplotype file file \n");
    else {
        j = 0;
        while (1) {
            fscanf(sf, "%s ", id);
            if (strcmp(id, "BLOCK:") != 0) break;
            fscanf(sf, "%s %d %s %d %s %d \n", id, &offset, id, &len, id, &phased);
            j++;
            for (i = 0; i < len; i++) fscanf(sf, "%s %c %c \n", id, &c1, &c2);
            fscanf(sf, "%s \n", id);
        }
        fclose(sf);

        blocks = j;
        blist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*blocks);
        sf = fopen(hapfile, "r");
        j = 0;
        while (1) {
            fscanf(sf, "%s ", id);
            if (strcmp(id, "BLOCK:") != 0) break; //fprintf(stdout,"%s %d\n",id,j-1);
            fscanf(sf, "%s %d %s %d %s %d \n", id, &offset, id, &len, id, &phased);
            blist[j].offset = offset - 1;
            blist[j].length = len;
            blist[j].phased = phased;
            //if (pflag) fprintf(stdout,"BLOCK--- %9d len %5d phased %5d \n",offset,len,phased); 
            j++;
            for (i = 0; i < len; i++) {
                fscanf(sf, "%s %c %c \n", id, &c1, &c2);
                if (c1 != '-') {
                    HAP1[offset + i - 1] = c1;
                    bn[offset + i - 1] = offset;
                    initHAP[offset + i - 1] = c1;
                    snpfrag[offset + i - 1].blockno = j;
                }
                strcpy(snpfrag[offset + i - 1].id, id); // IMPORTANT copy SNP id from haplotype solution to SNP ID 
                // offset is the id of each block since it is supposed to be unique  
            }
            fscanf(sf, "%s \n", id);
        }
        fclose(sf);
    }
    return 1;
    /***************************** READ HAPLOTYPE SOLUTION*************************************************/
}

