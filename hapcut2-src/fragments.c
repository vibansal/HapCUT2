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
    int i=0;
    for (i=0;i<fragments;i++) free_fragment(&Flist[i]);
    free(Flist);
}

int fragment_compare(const void *a, const void *b) { // used for sorting fragment list by position
    const struct fragment *ia = (const struct fragment*) a;
    const struct fragment *ib = (const struct fragment*) b;
    if (ia->list[0].offset == ib->list[0].offset) {
        return ia->list[ia->blocks - 1].offset + ia->list[ia->blocks - 1].len - ib->list[ib->blocks - 1].offset - ib->list[ib->blocks - 1].len;
        //return ia->blocks - ib->blocks;
    } else return ia->list[0].offset - ib->list[0].offset;
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

// filter the full fragment list to only keep fragments covering heterozygous variants and having at least 2 alleles covered
// the format of the new fragment list is not technically correct, each block is a singleton
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
	    		if ( (int)FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
			if (snpfrag[snp_ix].phase == '1') het++;
			total++;
		}
	   }
	   if (het < 2) continue; 
           //fprintf(stdout,"stats %d %d \n",total,het); 
	   // copy the ID, insert size, data_type and mate2_ix from original fragment
	   nFlist[n].id = calloc(sizeof(char),strlen(FRAG->id)+1); strcpy(nFlist[n].id,FRAG->id);
	   nFlist[n].isize = FRAG->isize; 
	   nFlist[n].mate2_ix = FRAG->mate2_ix;
   	   nFlist[n].blocks = het; nFlist[n].data_type = FRAG->data_type;
   	   nFlist[n].list = calloc(sizeof(struct block),nFlist[n].blocks); 
	   q=0;
	   
	   for (j = 0; j < FRAG->blocks; j++)
	   {
		for (k = 0; k < FRAG->list[j].len; k++)
		{
			snp_ix = FRAG->list[j].offset + k; // index of current position
	    		if ( (int)FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
			if (snpfrag[snp_ix].phase == '1') 
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
   
// output the block-ID (first SNP in block), haplotype-assignment and probability of assignment  | only for long-reads or linked reads
int assign_fragment_2hap(struct fragment* FRAG, struct SNPfrags* snpfrag,char* h)
{
    int j = 0, k = 0,alleles=0,snpid=0,component;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    char tag;

    alleles=0; p0=0; p1 =0,snpid=-1; component = -1;

    for (j = 0; j < FRAG->blocks; j++) {
	for (k = 0; k < FRAG->list[j].len; k++) {
	    if (h[FRAG->list[j].offset + k] == '-' || (int) FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
	    if (snpid < 0) snpid = FRAG->list[j].offset+k;  // assign the variant-id to the fragment  
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
   if (prob >=0.95 && snpid >=0) { 
	   component = snpfrag[snpid].component; 
	   if (component >= 0) // only for phased variants 
	   {
	     FRAG->PS = snpfrag[component].position; FRAG->HP = tag; 
	     return 1;
	   }
   }
   return 0;
}

// for each fragment, output the block-ID (first SNP in block), haplotype-assignment and probability of assignment  | only for long-reads or linked reads
void fragment_assignments(struct fragment* Flist,int fragments, struct SNPfrags* snpfrag,char* h,char* outfile)
{
    int f=0;
    int valid = 0;
    FILE* OUTFILE = fopen(outfile,"w");

    for (f=0;f<fragments;f++)
    {
	valid = assign_fragment_2hap(&Flist[f], snpfrag,h);
        if (valid ==1) fprintf(OUTFILE,"%s %d %c %d\n",Flist[f].id,Flist[f].PS,Flist[f].HP,Flist[f].PQ);
    }
    fclose(OUTFILE);
}

