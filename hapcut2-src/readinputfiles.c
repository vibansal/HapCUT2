#include "readinputfiles.h"
#include "common.h"
#include <assert.h>

extern float HOMOZYGOUS_PRIOR;
extern int NEW_FRAGFILE_FORMAT;

int get_num_fragments(char* fragmentfile)
{
    char buffer[MAXBUF];
     FILE* ff = fopen(fragmentfile, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open fragment file %s\n", fragmentfile);
        exit(0);
    }
    int fragments = 0;
    while (fgets(buffer, MAXBUF, ff) != NULL){
        if (!((buffer[0] == '0')&&(buffer[1] == ' ')))
            fragments++;
    }
    fclose(ff);
    return fragments;
}

/*
int fragment_compare(const void *a, const void *b) { // used for sorting fragment list by position
    const struct fragment *ia = (const struct fragment*) a;
    const struct fragment *ib = (const struct fragment*) b;
    if (ia->list[0].offset == ib->list[0].offset) {
        return ia->list[ia->blocks - 1].offset + ia->list[ia->blocks - 1].len - ib->list[ib->blocks - 1].offset - ib->list[ib->blocks - 1].len;
        //return ia->blocks - ib->blocks;
    } else return ia->list[0].offset - ib->list[0].offset;
}*/


int read_fragment_matrix(char* fragmentfile, struct fragment* Flist, int fragments,int OFFSET_flist) {
    int i = 0, j = 0, k = 0, t = 0, t1 = 0, done = 0;
    int blocks = 0, type = 0, l = 0, biter = 0, offset = 0,dtype=0,isize = 0;
    char buffer[MAXBUF];
    char blockseq[5000000];
    for (i=0;i<MAXBUF;i++) buffer[i] = 0;
    for (i=0;i<5000000;i++) blockseq[i] = 0;
    char ch;
    int num_fields;
    int expected_num_fields;

    FILE* ff = fopen(fragmentfile, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open fragment file \n");
        return -1;
    }

    for (i = OFFSET_flist; i < OFFSET_flist+fragments; i++) {
        //		fprintf(stdout,"%s \n",buffer);
        done = 0;
        while (!done){
            j = 0;
            ch = fgetc(ff);
            num_fields = 1;
            while (ch != '\n') {
                buffer[j] = ch;
                j++;
                ch = fgetc(ff);
                if (ch == ' '){
                    num_fields++;
                }
            }
            buffer[j] = '\0';

            // if there are 0 blocks then ignore this fragment completely
            done = !((buffer[0] == '0') && (buffer[1] == ' '));
        }

        Flist[i].data_type = 0; // default data type
        Flist[i].htrans_prob = -80;

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

                if (NEW_FRAGFILE_FORMAT){
                    expected_num_fields = (6 + 2*blocks);
                }else{
                    expected_num_fields = (3 + 2*blocks);
                }

                if (num_fields < expected_num_fields){
                    fprintf_time(stderr, "ERROR: Invalid fragment file, too few fields at line %d.\n",i);
                    if (NEW_FRAGFILE_FORMAT)
                        fprintf_time(stderr, "If this is Hi-C data, are you using the new format for Hi-C data with extractHAIRS --HiC 1 option?\n");
                    exit(1);
                }else if (num_fields > expected_num_fields){
                    fprintf_time(stderr, "ERROR: Invalid fragment file, too many fields at line %d.\n",i);
                    if (!NEW_FRAGFILE_FORMAT)
                        fprintf_time(stderr, "If the file is the new HiC-related format (from extractHAIRS --hic option), be sure to use --hic, --hic_htrans_file, or --nf HapCUT2 options.\n");
                    exit(1);
                }
            } else if (type == 1) // read the fragment id, changed to allow dynamic length feb202011
            {
                Flist[i].id = (char*) malloc(t + 1);
                for (l = 0; l < t; l++) Flist[i].id[l] = blockseq[l];
                Flist[i].id[l] = '\0';
                //Flist[i].id = (char*)malloc(1); Flist[i].id[0] = '0'; this doesnt reduce the memory requirement
                if (NEW_FRAGFILE_FORMAT)
                    type = 4; // read in data type next
                else
                    type = 2; // old format, skip right to reading alleles
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

                type = 2;
                biter++;
            } else if (type == 4){
                // read in data type (0=normal, 1=HiC)
                dtype = 0;
                for (l = 0; l < t; l++) {
                    dtype = 10 * dtype + (int) (blockseq[l] - 48);
                }
                Flist[i].data_type = dtype;

                if (dtype == 2){
                    fprintf(stderr,"ERROR: Unlinked 10X fragments given as input. Unlinked 10X Fragments must be processed with LinkFragments.py script.\n\nExample:\n extractHAIRS --bam <BAM> --vcf <VCF> --out unlinked_fragments.frag --10X 1;\npython LinkFragments.py --fragments unlinked_fragments.frag --bam <BAM> --vcf <VCF> --out linked_fragments.frag\n");
                    exit(1);
                }

                type = 5;
            } else if (type == 5){
                // read in the position of mate 2
                if (blockseq[0] == '-'){ // negative index means no mate 2
                    Flist[i].mate2_ix = -1;
                }else{
                    offset = 0;
                    for (l = 0; l < t; l++) {
                        offset = 10 * offset + (int) (blockseq[l] - 48);
                    }
                    Flist[i].mate2_ix = offset - 1;
                }

                type = 6;
            } else if (type == 6){

                // read in the absolute insert size
                if (blockseq[0] == '-'){ // negative index means no mate 2
                    Flist[i].isize = -1;
                }else{
                    isize = 0;
                    for (l = 0; l < t; l++) {
                        isize = 10 * isize + (int) (blockseq[l] - 48);
                    }
                    Flist[i].isize = isize;
                }

                type = 2; // go back to reading in fragment info
            }
            t = 0;
        }
    }
    fclose(ff);
    qsort(Flist, fragments, sizeof (struct fragment), fragment_compare);

    return fragments;
}
