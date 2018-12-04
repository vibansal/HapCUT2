
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

// 11/28/2018 
// should include indels for realignment even if we don't phase them... 
// translated from https://github.com/pjedge/longshot/blob/master/src/realignment.rs
// this will replace nw.c code, current nw.c code works in probability space using best alignment... 

//define states of HMM, emission probability vector and transition probability matrix 
// states = match, insertion, deletion 
// simple dynamic programming over 2-D matrix between v and w (read and haplotype) 




int ALLOW_END_GAPS = 0; // for anchored alignment...
int BAND_WIDTH = 20; // default 

typedef struct // probabilities are in logspace 
{
        int states; // match, insertion,deletion 
        double** TRS; // transition probabilities between states 
        double match; double mismatch; double insertion; double deletion; // emission probs
} Align_Params;

extern Align_Params AP;

#define logsum(a, b) (((a) > (b)) ? ((a) + log(1.0 + exp(b-a))) : ((b) + log(1.0 + exp( (a) - (b)))))

double logsum1(double a, double b)
{
        if (a > b) return a + log(1 + exp(b-a));
        else return b + log(1+exp(a-b));
}

// do we need anchors, just require the read-seq to be aligned end-2-end while allowing end-gaps for reference haplotype 
double sum_all_alignments(char* v, char* w,Align_Params* params, int min_band_width)
{
        int lv = strlen(v),lw = strlen(w);
        int len_diff = lv-lw; if (len_diff < 0) len_diff *= -1; // take absolute value of difference
        int band_width = min_band_width + len_diff;
	//fprintf(stdout,"diff %d bw %d \n",len_diff,band_width);

        double* lower_prev = calloc(sizeof(double),lw+1);
        double* middle_prev = calloc(sizeof(double),lw+1);
        double* upper_prev = calloc(sizeof(double),lw+1);
        double* lower_curr = calloc(sizeof(double),lw+1);
        double* middle_curr = calloc(sizeof(double),lw+1);
        double* upper_curr = calloc(sizeof(double),lw+1);

        double log0 = -100000;
	int i=0,j=0;
        for (i=0;i<lw+1;i++)
        {
                lower_prev[i]= log0;            middle_prev[i]= log0; upper_prev[i]= log0;
                lower_curr[i]= log0;            middle_curr[i]= log0; upper_curr[i]= log0;
        }
        middle_prev[0] = 0.0;

        if (ALLOW_END_GAPS ==1)
        {
                upper_prev[1] = params->TRS[0][2];
                for (j=2;j<lw+1;j++) upper_prev[j] = upper_prev[j-1] + params->TRS[2][2];
        }

        for (i=1;i<lv;i++)
        {
                int band_middle = (lw*i)/lv;
                int band_start = 1;
                if (band_middle - band_width/2 >= 1) band_start = band_middle-band_width/2;
                int band_end = lw;
                if (band_middle + band_width/2 <= lw) band_end = band_middle + band_width/2;

                if (ALLOW_END_GAPS ==1)
                {
                        if (band_start ==1)
                        {
                                middle_curr[0] = log0;
                                if (i==1) lower_curr[0] = params->TRS[0][1];
                                else lower_curr[0] = lower_prev[0] + params->TRS[1][1];
                        }
                }

                for (j=band_start;j<band_end+1;j++)
                {
                        int lower_continue = lower_prev[j] + params->TRS[1][1]; // insertion to insertion 
                        int lower_from_middle = middle_prev[j] + params->TRS[0][1]; // match to insertion 
                        lower_curr[j] = params->insertion + logsum(lower_continue,lower_from_middle);

                        int upper_continue = upper_curr[j-1] + params->TRS[2][2];
                        int upper_from_middle = middle_curr[j-1] + params->TRS[0][2];
                        upper_curr[j] = params->deletion + logsum(upper_continue,upper_from_middle);

                        int middle_from_lower = lower_prev[j-1] + params->TRS[1][0];
                        int middle_continue = middle_prev[j-1] + params->TRS[0][0];
                        int middle_from_upper = upper_prev[j-1] + params->TRS[2][0];
                        double s = logsum(middle_from_lower,middle_continue);
                        double match_emission = params->match;
                        if (v[i-1] != w[j-1]) match_emission = params->mismatch;
                        middle_curr[j] = match_emission + logsum(s,middle_from_upper);
                }
                for (j=band_start;j<band_end+1;j++)
                {
                        upper_prev[j] = upper_curr[j]; middle_prev[j] = middle_curr[j]; lower_prev[j] = lower_curr[j];
                }
                upper_curr[band_start]= log0; middle_curr[band_start] = log0; lower_curr[band_start] = log0;
        }

	double ll = middle_prev[lw];
        free(lower_prev); 
	free(middle_prev); 
	free(upper_prev);
        free(lower_curr); 
	free(middle_curr); 
	free(upper_curr);
        return ll/log(10);
}


void test_realignment()
{
	Align_Params AP; 
	AP.states = 3; 
	AP.TRS = calloc(sizeof(double*),AP.states); // match =0, ins =1, del = 2
	int i=0;
	for (i=0;i<AP.states;i++) AP.TRS[i] = calloc(sizeof(double),AP.states);
	AP.match = log(0.979); AP.mismatch = log(0.007); AP.deletion = log(1); AP.insertion =log(1); 
	AP.TRS[0][0] = log(0.879); AP.TRS[0][1] = log(0.076); AP.TRS[0][2] = log(0.045); 
	AP.TRS[1][0] = log(0.865); AP.TRS[1][1] = log(0.135); 
	AP.TRS[2][0] = log(0.730); AP.TRS[2][2] = log(0.27); 

	char* ref = malloc(1024); char* alt = malloc(1024); char* read = malloc(1024);
	strcpy(ref,"GCTGGTGTAATGCAATG"); 	strcpy(alt,"GCTGGTGCAATGCAATG");  	strcpy(read,"GCTGGTTAATGCAATG");

	// real example where sum-all-paths makes a difference
	strcpy(read,"GGGCAGGCCCGCTGAG"); 	strcpy(ref,"GGGCAGCCCCTGAG");  	strcpy(alt,"GGGCAGCCGCTGAG");

	double score0 =  sum_all_alignments(ref,read,&AP,BAND_WIDTH);
	double score1 =  sum_all_alignments(alt,read,&AP,BAND_WIDTH);
	double scoreS = logsum(score0,score1);
	double phred;
	if (score0 > score1)  phred = -10.0*(score1-scoreS)/log(10); 
	else phred = -10.0*(score0-scoreS)/log(10);
	fprintf(stdout,"match %f %f %f \n",AP.TRS[0][0],AP.TRS[0][1],AP.TRS[0][2]);
	fprintf(stdout,"ins %f %f del %f %f\n",AP.TRS[1][0],AP.TRS[1][1],AP.TRS[2][0],AP.TRS[2][2]);
	fprintf(stdout,"code for pair HMM realignment of two sequences %f %f %f\n",log10(exp(score0)),log10(exp(score1)),phred);
}

/*
int main(int argc, char** argv) 
{
	test_realignment();
	return 0;
}
*/

