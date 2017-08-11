#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <stddef.h>
#include<float.h>
#include<ctype.h>

/*
// the scoring scheme is: match=2, mismatch=-5, gapOpen=-2 and gapExt=-1.  BWA-MEM
#define MATCH 2     //-0.0044
#define MISMATCH -5  //-2
#define GAP_OPEN -2 //-3
#define GAP_EXTEND -1 //-1
*/

extern int VERBOSE;

extern float MATCH;
extern float MISMATCH;
extern float GAP_OPEN;
extern float GAP_EXTEND;
extern float INSERTION_OPEN;
extern float INSERTION_EXTEND;
extern float DELETION_OPEN;
extern float DELETION_EXTEND;

////////////////////////////////////////////////////////////////////////////////////////////////////

//dynamic programming: needleman-wunsch for global alignment of sequences, with affine gap penalty
//str1 is the reference haplotype, str2 is the read
double nw(char* str1, char* str2,int VERBOSE)
{
	if (VERBOSE) printf("applying alignment:\n");

	// gap characters are appended so that score can be initialized to infinity
	int i,j;
	char* h = calloc(strlen(str1)+2,sizeof(char)); /* prepend gap chars*/
	sprintf(h, "%c%s", '_', str1);
	int len_h = strlen(h); //define to save on fn calls
	if (VERBOSE) printf("h: %s\n", h);
	char* v = calloc(strlen(str2)+2,sizeof(char));
	sprintf(v, "%c%s", '_', str2);
	int len_v = strlen(v);
	if (VERBOSE) printf("v: %s\n", v);

	typedef struct {
		double v; //score value
		int x;
		int y;
		char w; //where did we come from?
	} xy;

	xy** m = calloc(len_v, sizeof(xy*)); //backtracking
	xy** ix = calloc(len_v, sizeof(xy*));
	xy** iy = calloc(len_v, sizeof(xy*));

	for(i = 0; i < len_v; i++) { //calloc
		m[i] = calloc(len_h, sizeof(xy));
		ix[i] = calloc(len_h, sizeof(xy));
		iy[i] = calloc(len_h, sizeof(xy));
	}

	int full =0; // to decide if alignment is end-to-end or not..
	if(full == 0)
	{
		for(i = 0; i < len_v; i++) for(j = 0; j < len_h; j++) m[i][j].v = ix[i][j].v = iy[i][j].v = -INFINITY;
		m[0][0].v = 0;
		for(j = 0; j < len_h; j++) iy[0][j].v = DELETION_OPEN + DELETION_EXTEND*j;
		for(i = 0; i < len_v; i++) ix[i][0].v = INSERTION_OPEN + INSERTION_EXTEND*i;
	}
	else
	{
	//	for(i = 0; i < len_v; i++) m[i][0].v = ix[i][0].v = iy[i][0].v = -INFINITY;
	}

	int margin = 1000; //don't implement banding for full alignment

	//dynamic programming loop
	for(j = 1; j < len_h; j++) {
		for(i = 1; i < len_v; i++) {
			//implement a margin....
			if(j <= i+margin && j >= i-margin) {
				double max_mij = -INFINITY;
				double s = (toupper(h[j]) == toupper(v[i])) ? MATCH : MISMATCH; //do bases match??

				/*m*/
				if(m[i-1][j-1].v + s >= max_mij) {
					max_mij = m[i-1][j-1].v + s;
					m[i][j].w = 'm';
				}
				if(ix[i-1][j-1].v + s >= max_mij) {
					max_mij = ix[i-1][j-1].v +s;
					m[i][j].w = 'x';
				}
				if(iy[i-1][j-1].v + s >= max_mij) {
					max_mij = iy[i-1][j-1].v +s;
					m[i][j].w = 'y';
				}
				m[i][j].v = max_mij;

				/*I_x*/
				double max_ix = -INFINITY;
				if(m[i-1][j].v + INSERTION_OPEN >= max_ix) {
					max_ix = m[i-1][j].v + INSERTION_OPEN + INSERTION_EXTEND;
					ix[i][j].w = 'm';
				}
				if(ix[i-1][j].v + INSERTION_EXTEND >= max_ix) {
					max_ix = ix[i-1][j].v + INSERTION_EXTEND;
					ix[i][j].w = 'x';
				}
				ix[i][j].v = max_ix;

				/*I_y */
				double max_iy = -INFINITY;
				if(m[i][j-1].v + DELETION_OPEN >= max_iy) {
					max_iy = m[i][j-1].v + DELETION_OPEN + DELETION_EXTEND;
					iy[i][j].w = 'm';
				}
				if(iy[i][j-1].v + DELETION_EXTEND >= max_iy) {
					max_iy = iy[i][j-1].v + DELETION_EXTEND;
					iy[i][j].w = 'y';
				}
				iy[i][j].v = max_iy;

			}

		}
	}

	/*alteration: now find the largest value in the last row (could be M, I_x or I_y) */
	//and include the margin explicitly in this: the result should be the same
	i = len_v-1;
	xy max_mat; max_mat.v = -INFINITY; max_mat.x = i;
	//char w;

	//    for(j = 0; j < len_h; j++) {
	//we want to be within the margin and also within the matrix
	for(j = (i-margin > 0 ? i-margin : 0); j < (i+margin < len_h ? i+margin : len_h); j++) { /* NEW! */
		//            printf("consider (%d,%d)\n", i, j);
		if(m[i][j].v >= max_mat.v) {
			max_mat.y = j; max_mat.v = m[i][j].v; max_mat.w = m[i][j].w;
			//w = 'm';
		}
		if(ix[i][j].v >= max_mat.v) {
			max_mat.y = j; max_mat.v = ix[i][j].v; max_mat.w = ix[i][j].w;
			//w = 'x';
		}
		if(iy[i][j].v >= max_mat.v) {
			max_mat.y = j; max_mat.v = iy[i][j].v; max_mat.w = iy[i][j].w;
			//w = 'y';
		}
	}

	//al->score = max_mat.v;
	if (VERBOSE) fprintf(stdout,"alignment score %f \n",max_mat.v);

	//dealloc
	for(i = 0; i < len_v; i++) {
		free(m[i]); free(ix[i]); free(iy[i]);
	}
	free(m); free(ix); free(iy); free(h), free(v);
	return max_mat.v; // the alignment score
	//no return
}
