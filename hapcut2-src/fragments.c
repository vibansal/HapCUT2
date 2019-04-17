#include "fragments.h"
#include "common.h"
#include "frag_likelihood.c"

void free_fragment(struct fragment* FRAG)
{
    int j=0;
    for (j = 0; j < FRAG->blocks; j++){
    free(FRAG->list[j].hap);
            free(FRAG->list[j].qv);
            free(FRAG->list[j].pv);
            free(FRAG->list[j].p1);
    }
    free(FRAG->list); free(FRAG->id);
}

void free_fragmentlist(struct fragment* Flist, int fragments)
{
    int i=0,j=0;
    for (i=0;i<fragments;i++) free_fragment(&Flist[i]);
    free(Flist);
}

// simply print in the same format as the input file
void print_fragment(struct fragment* FRAG,FILE* OUTFILE)
{
   int j=0,k=0;
   fprintf(OUTFILE,"%d %s ",FRAG->blocks,FRAG->id);
   for (j = 0; j < FRAG->blocks; j++) {
	fprintf(OUTFILE,"%d ",FRAG->list[j].offset); 
	for (k = 0; k < FRAG->list[j].len; k++) fprintf(OUTFILE,"%c",FRAG->list[j].hap[k]); 
	fprintf(OUTFILE," ");
   } 
   for (j = 0; j < FRAG->blocks; j++) { 
	for (k = 0; k < FRAG->list[j].len; k++) fprintf(OUTFILE,"%c",FRAG->list[j].qv[k]); 
   }
   fprintf(OUTFILE,"\n");
}

// current output format is different, no blocks but list of variants, bad results for Hi-C (understandable)
int filter_fragments(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,struct fragment* nFlist)
{
   int j=0,k=0,snp_ix=0,i=0;
   int total=0,het=0;
   int n=0;
   int q =0;

   struct fragment* FRAG;
   for (i=0;i<fragments;i++)
   {
	FRAG = &Flist[i]; het=0; total =0;
	   for (j = 0; j < FRAG->blocks; j++) 
	   {
		for (k = 0; k < FRAG->list[j].len; k++) 
		{
			snp_ix = FRAG->list[j].offset + k; // index of current position
			if (snpfrag[snp_ix].ignore == '0') het++;
			total++;
		}
	   }
	   if (het < 2) continue; 
           //fprintf(stdout,"stats %d %d \n",total,het); 

	   nFlist[n].id = calloc(sizeof(char),strlen(FRAG->id)+1); strcpy(nFlist[n].id,FRAG->id);
   	   nFlist[n].list = calloc(sizeof(struct block),het); 
   	   nFlist[n].blocks = het; nFlist[n].data_type = FRAG->data_type;
	   q=0;
	   
	   for (j = 0; j < FRAG->blocks; j++)
	   {
		for (k = 0; k < FRAG->list[j].len; k++)
		{
			snp_ix = FRAG->list[j].offset + k; // index of current position
			if (snpfrag[snp_ix].ignore == '0') 
			{
				nFlist[n].list[q].offset = snp_ix; nFlist[n].list[q].len = 1;
				nFlist[n].list[q].hap = calloc(sizeof(char),2); nFlist[n].list[q].hap[0] = FRAG->list[j].hap[k];
				nFlist[n].list[q].qv = calloc(sizeof(char),2); nFlist[n].list[q].qv[0] = FRAG->list[j].qv[k];
				nFlist[n].list[q].pv = calloc(sizeof(float),1); nFlist[n].list[q].pv[0] = FRAG->list[j].pv[k];
				nFlist[n].list[q].p1 = calloc(sizeof(float),1); nFlist[n].list[q].p1[0] = FRAG->list[j].p1[k];
				q++;
			}
		}
	   }
           nFlist[n].calls = q;
           //fprintf(stdout,"NEW "); print_fragment(&nFlist[n],stdout);
	   n++;
   }
   return n;
}
   


/*
   fprintf(OUTFILE,"%s %d %c %d\n",FRAG->id,FRAG->PS,FRAG->HP,FRAG->PQ);
   fprintf(OUTFILE,"%d %c %d %f %f ",FRAG->PS,FRAG->HP,FRAG->PQ,p0,p1);
   fprintf(OUTFILE,"%d %s ",FRAG->blocks,FRAG->id); 
   */

// output the block-ID (first SNP in block), haplotype-assignment and probability of assignment  | only for long-reads or linked reads
int fragment_assignment(struct fragment* FRAG, struct SNPfrags* snpfrag,char* h)
{
    int f=0,j = 0, k = 0,alleles=0,offset=0,component;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    char tag;

    alleles=0; p0=0; p1 =0,offset=-1; component = -1;

    for (j = 0; j < FRAG->blocks; j++) {
	for (k = 0; k < FRAG->list[j].len; k++) {
	    if (h[FRAG->list[j].offset + k] == '-' || (int) FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
	    if (offset < 0) offset = FRAG->list[j].offset; 
	    prob = (QVoffset - (int) FRAG->list[j].qv[k]); prob /= 10;
	    prob2 = FRAG->list[j].p1[k];
	    alleles++;

	    if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) {
		p0 += prob2;
		p1 += prob;
	    } else {
		p0 += prob;
		p1 += prob2;
	    }
	}
   }
   if (p0 > p1) { tag ='0'; prob = pow(10,p0-addlogs(p0,p1)); FRAG->PQ = (int)(10*(addlogs(p0,p1)-p1)); }
   else { tag = '1'; prob = pow(10,p1-addlogs(p0,p1)); FRAG->PQ = (int)(10*(addlogs(p0,p1)-p0)); } 
   if (prob >=0.9 && offset >=0) { 
	   component = snpfrag[offset].component; // unphased SNPs should be ignored
	   // print PS,hap(0|1),probability followed by original fragment copy
	   if (component < 0) { fprintf(stderr,"error \n"); } 
	   else
	   {
	     FRAG->PS = snpfrag[component].position; FRAG->HP = tag; 
	   }
	   return 1;
   }
   return 0;
}

// for each fragment, output the block-ID (first SNP in block), haplotype-assignment and probability of assignment  | only for long-reads or linked reads
void fragment_assignments(struct fragment* Flist,int fragments, struct SNPfrags* snpfrag,char* h,char* outfile)
{
    int f=0,j = 0, k = 0,alleles=0,offset=0,component;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    char tag;
    int valid = 0;
    FILE* OUTFILE = fopen(outfile,"w");

    for (f=0;f<fragments;f++)
    {
	valid = fragment_assignment(&Flist[f], snpfrag,h);
        if (valid ==1) fprintf(OUTFILE,"%s %d %c %d\n",Flist[f].id,Flist[f].PS,Flist[f].HP,Flist[f].PQ);
    }
    fclose(OUTFILE);
}

