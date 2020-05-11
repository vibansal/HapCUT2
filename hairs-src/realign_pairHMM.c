
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

// can do without log-space since sum will not be very small...

// 11/28/2018 
// translated from https://github.com/pjedge/longshot/blob/master/src/realignment.rs

//define states of HMM, emission probability vector and transition probability matrix 
// states = match, insertion, deletion 
// simple dynamic programming over 2-D matrix between v and w (read and haplotype) 

int ALLOW_END_GAPS = 0; // for anchored alignment...
int BAND_WIDTH = 20; // default 

//#define logsum(a, b) (((a) > (b)) ? ((a) + log(1.0 + exp(b-a))) : ((b) + log(1.0 + exp( (a) - (b)))))
#define logsum(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))

double logsum1(double a, double b)
{
	if (a > b) return a + log10(1.0 + pow(10,b-a));
	else return b + log10(1.0+pow(10,a-b));
}

// do we need anchors, just require the read-seq to be aligned end-2-end while allowing end-gaps for reference haplotype 
double sum_all_alignments_logspace(char* v, char* w,Align_Params* params, int min_band_width)
{
	int lv = strlen(v),lw = strlen(w);
	int len_diff = lv-lw; if (len_diff < 0) len_diff *= -1; // take absolute value of difference
	int band_width = min_band_width + len_diff;
	int b1,b2;
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

	//if (ALLOW_END_GAPS ==1)
	{
		upper_prev[1] = params->lTRS[0][2];
		for (j=2;j<lw+1;j++) 
		{
			upper_prev[j] = upper_prev[j-1] + params->lTRS[2][2];
			middle_prev[j] = log0;
		}
	}

	for (i=1;i<lv;i++)
	{
		int band_middle = (lw*i)/lv;
		int band_start = 1;
		if (band_middle - band_width/2 >= 1) band_start = band_middle-band_width/2;
		int band_end = lw;
		if (band_middle + band_width/2 <= lw) band_end = band_middle + band_width/2;

		//if (ALLOW_END_GAPS ==1)
		{
			if (band_start ==1)
			{
				middle_curr[0] = log0;
				if (i==1) lower_curr[0] = params->lTRS[0][1];
				else lower_curr[0] = lower_prev[0] + params->lTRS[1][1];
			}
		}

		for (j=band_start;j<band_end+1;j++)
		{
			double lower_continue = lower_prev[j] + params->lTRS[1][1]; // insertion to insertion 
			double lower_from_middle = middle_prev[j] + params->lTRS[0][1]; // match to insertion 
			lower_curr[j] = params->linsertion + logsum(lower_continue,lower_from_middle);

			double upper_continue = upper_curr[j-1] + params->lTRS[2][2];
			double upper_from_middle = middle_curr[j-1] + params->lTRS[0][2];
			upper_curr[j] = params->ldeletion + logsum(upper_continue,upper_from_middle);

			double middle_from_lower = lower_prev[j-1] + params->lTRS[1][0];
			double middle_continue = middle_prev[j-1] + params->lTRS[0][0];
			double middle_from_upper = upper_prev[j-1] + params->lTRS[2][0];

			double s = logsum(middle_from_lower,middle_continue);
			double match_emission = params->lmatch;
			if (v[i-1] != w[j-1]) match_emission = params->lmismatch;
			b1 = BTI[(int)v[i-1]]-1; 
                        b2 = BTI[(int)w[j-1]]-1; 
			if (b1 >= 0 && b2 < 4 && b2 >= 0 && b2 < 4) 
			{
				//fprintf(stdout,"match %d %d %f \n",b1,b2,params->MEM[b1][b2]);
				match_emission = params->lMEM[b1][b2]; 
			}
			middle_curr[j] = match_emission + logsum(s,middle_from_upper);
		}
		for (j=band_start;j<band_end+1;j++)
		{
			upper_prev[j] = upper_curr[j]; middle_prev[j] = middle_curr[j]; lower_prev[j] = lower_curr[j];
		}
		upper_curr[band_start]= log0; middle_curr[band_start] = log0; lower_curr[band_start] = log0;
	}

	double ll = middle_prev[lw];
	free(lower_prev); 	free(middle_prev);  	free(upper_prev);
	free(lower_curr);  	free(middle_curr);  	free(upper_curr);
	//fprintf(stderr,"v %s w %s %f \n",v,w,ll/log(10));
	return ll;
}

// direct multiplication instead of logspace, double limited to 
double sum_all_alignments_fast(char* v, char* w,Align_Params* params, int min_band_width)
{
	int lv = strlen(v),lw = strlen(w);
	int len_diff = lv-lw; if (len_diff < 0) len_diff *= -1; // take absolute value of difference
	int band_width = min_band_width + len_diff;
	int b1,b2;
	//fprintf(stdout,"diff %d bw %d \n",len_diff,band_width);

	double* lower_prev = calloc(sizeof(double),lw+1);
	double* middle_prev = calloc(sizeof(double),lw+1);
	double* upper_prev = calloc(sizeof(double),lw+1);
	double* lower_curr = calloc(sizeof(double),lw+1);
	double* middle_curr = calloc(sizeof(double),lw+1);
	double* upper_curr = calloc(sizeof(double),lw+1);

	double log0 = 0.0;
	int i=0,j=0;
	for (i=0;i<lw+1;i++)
	{
		lower_prev[i]= log0;            middle_prev[i]= log0; upper_prev[i]= log0;
		lower_curr[i]= log0;            middle_curr[i]= log0; upper_curr[i]= log0;
	}
	middle_prev[0] = 1.0;

	//if (ALLOW_END_GAPS ==1)
	{
		upper_prev[1] = params->TRS[0][2];
		for (j=2;j<lw+1;j++) 
		{
			upper_prev[j] = upper_prev[j-1] * params->TRS[2][2];
			middle_prev[j] =log0;
		}
	}

	for (i=1;i<lv;i++) // main loop
	{
		int band_middle = (lw*i)/lv;
		int band_start = 1;
		if (band_middle - band_width/2 >= 1) band_start = band_middle-band_width/2;
		int band_end = lw;
		if (band_middle + band_width/2 <= lw) band_end = band_middle + band_width/2;

		//if (ALLOW_END_GAPS ==1)
		{
			if (band_start ==1)
			{
				middle_curr[0] = log0; upper_curr[0] = log0; // different between two functions ?? check 
				if (i==1) lower_curr[0] = params->TRS[0][1];
				else lower_curr[0] = lower_prev[0] * params->TRS[1][1];
			}
		}

		for (j=band_start;j<band_end+1;j++)
		{
			double lower_continue = lower_prev[j] * params->TRS[1][1]; // insertion to insertion 
			double lower_from_middle = middle_prev[j] * params->TRS[0][1]; // match to insertion 
			lower_curr[j] = params->insertion * (lower_continue+lower_from_middle);

			double upper_continue = upper_curr[j-1] * params->TRS[2][2];
			double upper_from_middle = middle_curr[j-1] * params->TRS[0][2];
			upper_curr[j] = params->deletion * (upper_continue+upper_from_middle);

			double middle_from_lower = lower_prev[j-1] * params->TRS[1][0];
			double middle_continue = middle_prev[j-1]* params->TRS[0][0];
			double middle_from_upper = upper_prev[j-1] * params->TRS[2][0];

			double s = (middle_from_lower+middle_continue);
			double match_emission = params->match;
			if (v[i-1] != w[j-1]) match_emission = params->mismatch;
			b1 = BTI[(int)v[i-1]]-1; b2 = BTI[(int)w[j-1]]-1; 
			if (b1 >= 0 && b2 < 4 && b2 >= 0 && b2 < 4) 
			{
				//fprintf(stdout,"match %d %d %f \n",b1,b2,params->MEM[b1][b2]);
				match_emission = params->MEM[b1][b2]; 
			}
			middle_curr[j] = match_emission * (s+middle_from_upper);
		}
		for (j=band_start;j<band_end+1;j++)
		{
			upper_prev[j] = upper_curr[j]; middle_prev[j] = middle_curr[j]; lower_prev[j] = lower_curr[j];
		}
		if (band_start >=2)
		{
			upper_prev[band_start-2] = -1000000;
	                middle_prev[band_start-2] = -1000000;
                        lower_prev[band_start-2] = -1000000; // NaN
		}

		upper_curr[band_start]= log0; middle_curr[band_start] = log0; lower_curr[band_start] = log0;
	}
	double ll;
	if (middle_prev[lw] ==0) // if underflow has occurred, middle_prev[lw] will be == 0
	{
		ll = sum_all_alignments_logspace(v,w,params,min_band_width);
	}
	else ll = log10(middle_prev[lw]);
	free(lower_prev); 	free(middle_prev);  	free(upper_prev);
	free(lower_curr);  	free(middle_curr);  	free(upper_curr);
	//fprintf(stderr,"v %s w %s %f \n",v,w,ll/log(10));
	return ll;
}

