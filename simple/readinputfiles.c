#include "readinputfiles.h"
extern int HOMOZYGOUS_PRIOR;

int fragment_compare(const void *a, const void *b) {
    const struct fragment *ia = (const struct fragment*) a;
    const struct fragment *ib = (const struct fragment*) b;
    if (ia->list[0].offset == ib->list[0].offset) {
        return ia->list[ia->blocks - 1].offset + ia->list[ia->blocks - 1].len - ib->list[ib->blocks - 1].offset - ib->list[ib->blocks - 1].len;
        //return ia->blocks - ib->blocks;
    } else return ia->list[0].offset - ib->list[0].offset;
}


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
    char GQ[5];
    char* gen;
    int i = 0, j = 0, k=0, len=0, s = 0, e = 0, var = 0, GQ_ix, format_ix;
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
        len = j-s;
        gen = snpfrag[var].genotypes; //  for convenience
        if (GQ_ix != -1) {
            // get to the index where GQ is located.
            k = 0; // where we are in gen
            for (j=0; j<GQ_ix; j++){
                while(gen[k] != ':')
                    k++;
                k++; // step past the ':'
            }
            // reached GQ field. read it in.
            
            j=0;
            while (k<len && gen[k] != ':') {
                GQ[j] = gen[k];
                k++;
                j++;
            }
            GQ[j] = '\0';

            snpfrag[var].homozygous_prior = atof(GQ) / -10; // log prior probability of homozygousity
        } else{
            snpfrag[var].homozygous_prior = HOMOZYGOUS_PRIOR;
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
        fprintf(stderr, "couldn't open variant file %s\n", variantfile);
        exit(0);
    }
    while (fgets(buffer, MAXBUF, ff) != NULL) snps++;
    fclose(ff);
    return snps;
}
