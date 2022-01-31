#include "hic.h"
#include <assert.h>
float HIC_EM_THRESHOLD = 0.99; 
/*
functions that are specific to processing of Hi-C reads for phasing: estimating h-trans error rates and I/O
HTRANS_BINSIZE = 5000;
HTRANS_MAXBINS = 10000; // this value will be overwritten at startup
HTRANS_READ_LOWBOUND = 500;
HTRANS_MAX_WINDOW = 4000000; // maximum window size for h-trans estimation
*/

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

// estimate probabilities of h-trans to feed back into HapCUT algorithm, HIC data 
int estimate_htrans_probs(struct fragment* Flist, int fragments, char* HAP, struct SNPfrags* snpfrag,char* htrans_OUTFILE){

    float* MLE_sum   = calloc(HTRANS_MAXBINS,sizeof(float));
    float* MLE_count = calloc(HTRANS_MAXBINS,sizeof(float));
    float* p_htrans      = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_sum   = calloc(HTRANS_MAXBINS,sizeof(float));
    float* adj_MLE_count = calloc(HTRANS_MAXBINS,sizeof(float));

    int i=0,j=0,k=0,f=0,bin;
    int i_minus = 0, i_plus = 0;
    int e_window_size = HTRANS_BINSIZE; //track the effective window size

    // consistent counts, inconsistent counts
    int i1=-1, i2=-1;
    char a1='-', a2='-', h1='-', h2='-';
    float q1=0,q2=0;
    // count the number of consistent vs inconsistent for a given insert size bin
    for (f = 0; f < fragments; f++){

        // consider mate pairs only
        if (Flist[f].mate2_ix == -1 || Flist[f].isize == -1){
            continue;
        }

        // insert size bin
        bin = Flist[f].isize / HTRANS_BINSIZE;

        if (bin < 0){
            fprintf_time(stderr,"ERROR: bin less than 0");
            exit(1);
        }
        if (bin >= HTRANS_MAXBINS){
            fprintf_time(stderr,"ERROR: bin greater than HTRANS_MAXBINS");
            exit(1);
        }

        // keep things very simple by only sampling 1-snp mates
        if (Flist[f].calls < 2){
            continue;
        }

        i1 = Flist[f].list[0].offset;
        a1 = Flist[f].list[0].hap[0];
        q1 = Flist[f].list[0].pv[0];

        i2 = Flist[f].mate2_ix;
        for (j=0; j<Flist[f].blocks; j++){
            for (k=0; k<Flist[f].list[j].len; k++){
                if (Flist[f].list[j].offset+k == Flist[f].mate2_ix){
                    a2 = Flist[f].list[j].hap[k];
                    q2 = Flist[f].list[j].pv[k];
                    break;
                }
            }
        }

        h1 = HAP[i1];
        h2 = HAP[i2];

        if (h1 == '-' || h2 == '-'
         || snpfrag[i1].post_hap < log10(HIC_EM_THRESHOLD)|| snpfrag[i2].post_hap < log10(HIC_EM_THRESHOLD)
         || (snpfrag[i1].bcomp != snpfrag[i2].bcomp)
         || (snpfrag[i1].bcomp == -1)
         || (snpfrag[i2].bcomp == -1)
         || (snpfrag[i1].frags == 1)
         || (snpfrag[i2].frags == 1)){
            continue;
        }

        MLE_count[bin]++;

        assert(i1 < Flist[f].mate2_ix);
        assert(i2 == Flist[f].mate2_ix);
        assert (h1 == '1' || h1 == '0');
        assert (h2 == '1' || h2 == '0');
        assert (a1 == '1' || a1 == '0');
        assert (a2 == '1' || a2 == '0');

        if ((a1 == a2) == (h1 == h2)){ // phase match
            MLE_sum[bin] += ((1-q1)*q2 + (1-q2)*q1);
        }else{                         // phase mismatch
            MLE_sum[bin] += ((1-q1)*(1-q2) + q1*q2);
        }
    }

    // using neighboring bins to calculate mean htrans for each bin if the number of reads in bin is less < 500 (READ_LOWBOUND)
    for (i = 0; i < HTRANS_MAXBINS; i++){
        adj_MLE_count[i] = MLE_count[i];
        adj_MLE_sum[i] = MLE_sum[i];
        i_minus = i;
        i_plus = i;
        e_window_size = HTRANS_BINSIZE; //track the effective window size
        for (j = 0; j< 100000; j++){
            if (adj_MLE_count[i] >= HTRANS_READ_LOWBOUND) break; // 500 observations only ??
            i_minus--;
            i_plus++;
            if (i_minus >= 0){
                adj_MLE_count[i] += MLE_count[i_minus];
                adj_MLE_sum[i] += MLE_sum[i_minus];
                e_window_size += HTRANS_BINSIZE;
            }
            if(i_plus < HTRANS_MAXBINS){
                adj_MLE_count[i] += MLE_count[i_plus];
                adj_MLE_sum[i] += MLE_sum[i_plus];
                e_window_size += HTRANS_BINSIZE;
            }
            if (e_window_size >= HTRANS_MAX_WINDOW) break; // cut off window expansion if it's larger than some amount
        }
	//fprintf(stdout,"i %d %d-%d count %f %f \n",i,i_minus*HTRANS_BINSIZE,i_plus*HTRANS_BINSIZE,adj_MLE_count[i],adj_MLE_sum[i]);
    }

    // compute the MLE for each bin
    for (i = 0; i < HTRANS_MAXBINS; i++){
        if (adj_MLE_count[i] == 0 || adj_MLE_sum[i] == 0)
            p_htrans[i] = -80;
        else
            p_htrans[i] = log10(adj_MLE_sum[i] / adj_MLE_count[i]);
    }

    // assign the probabilities to fragments based in insert size
    for (f = 0; f < fragments; f++){
        if (Flist[f].mate2_ix == -1){
            Flist[f].htrans_prob = -80;
            continue;
        }
        Flist[f].htrans_prob = p_htrans[Flist[f].isize / HTRANS_BINSIZE];
    }

    if (strcmp(htrans_OUTFILE,"None") != 0){
        FILE* fp;
        fp = fopen(htrans_OUTFILE, "w");
        // compute the MLE for each bin
        for (i = 0; i < HTRANS_MAXBINS; i++){
            fprintf(fp,"%d-%d\t%f\n",i*HTRANS_BINSIZE,(i+1)*HTRANS_BINSIZE,
                    pow(10,p_htrans[i]));
        }
        fclose(fp);
    }

    free(MLE_sum);
    free(MLE_count);
    free(adj_MLE_sum);
    free(adj_MLE_count);
    free(p_htrans);
    return 0;
}
