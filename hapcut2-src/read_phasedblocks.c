#include "common.h"


/****************************** READ HAPLOTYPE SOLUTION BLOCK by BLOCK*************************************************/
int read_haplotypefile(char* hapfile, struct SNPfrags* snpfrag, int snps, char* HAP1, int* bn) {
    int i = 0, j = 0;
    char id[100];
    char c1, c2;
    int offset, len, phased, blocks;
    FILE* sf = fopen(hapfile, "r");
    struct BLOCK* blist;
    if (sf == NULL) fprintf_time(stderr, "couldn't open initial haplotype file file \n");

    else {
        // obtain the number of blocks 
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

        // now read the actual blocks and phase for each variant
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

