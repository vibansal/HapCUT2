#include "hic.h"
#include <assert.h>

int count_htrans_bins(char* htrans_file) {
    int bins = 0;
    char buffer[MAXBUF];
    FILE* ff = fopen(htrans_file, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open htrans file %s\n", htrans_file);
        exit(0);
    }
    while (fgets(buffer, MAXBUF, ff) != NULL) bins++;
    fclose(ff);
    return bins;
}

int read_htrans_file(char* htrans_file, float* htrans_probs, int num_bins) {
    int i = 0, j = 0;

    char buffer[MAXBUF];
    char ch;

    FILE* ff = fopen(htrans_file, "r");
    if (ff == NULL) {
        fprintf_time(stderr, "couldn't open htrans file \n");
        return -1;
    }

    for (i = 0; i < num_bins; i++) {

        // read bin size into bins
        j = 0;
        ch = fgetc(ff);
        while (ch != '\t') {
            buffer[j] = ch;
            j++;
            ch = fgetc(ff);
        }
        buffer[j] = '\0';
        //bins[i] = atoi(buffer);

        // read htrans probability into htrans_probs
        j = 0;
        while (ch != '\n') {
            buffer[j] = ch;
            j++;
            ch = fgetc(ff);
        }
        buffer[j] = '\0';
        htrans_probs[i] = atof(buffer);

    }

    return 0;
}


void init_HiC(struct fragment* Flist,int fragments,char* htrans_data_file)
{
        int MAXIS = -1;
        int i =0;
        // determine the probability of an h-trans interaction for read
        for (i=0; i<fragments;i++){

            Flist[i].htrans_prob = -80;

            if (Flist[i].isize > MAXIS)
                MAXIS = Flist[i].isize;
        }
        HTRANS_MAXBINS = MAXIS/HTRANS_BINSIZE + 1;

    // read in file with estimated probabilities of Hi-C h-trans interactions with distance
    if (strcmp(htrans_data_file, "None") != 0){
        int num_bins        = count_htrans_bins(htrans_data_file);
        float* htrans_probs = (float*) malloc(sizeof(float) * num_bins);
        read_htrans_file(htrans_data_file, htrans_probs, num_bins);
        for (i=0; i<fragments;i++)      Flist[i].htrans_prob = log10(htrans_probs[Flist[i].isize / HTRANS_BINSIZE]);
        free(htrans_probs);
        }
}

