#include "common.h"

// for each fragment, output the block-ID (first SNP in block), haplotype-assignment and probability of assignment  | only for long-reads or linked reads
void fragment_assignments(struct fragment* Flist,int fragments, struct SNPfrags* snpfrag,char* h,char* outfile)
{
    int f=0,j = 0, k = 0,alleles=0,offset=0,component;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    char tag;
    FILE* OUTFILE = fopen(outfile,"w");

    for (f=0;f<fragments;f++)
    {
	    alleles=0; p0=0; p1 =0,offset=-1; component = -1;

	    for (j = 0; j < Flist[f].blocks; j++) {
		for (k = 0; k < Flist[f].list[j].len; k++) {
		    if (h[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
		    if (offset < 0) offset = Flist[f].list[j].offset; 
		    prob = (QVoffset - (int) Flist[f].list[j].qv[k]); prob /= 10;
		    prob2 = Flist[f].list[j].p1[k];
		    alleles++;

		    if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
			p0 += prob2;
			p1 += prob;
		    } else {
			p0 += prob;
			p1 += prob2;
		    }
		}
	   }
	   if (p0 > p1) { tag ='0'; prob = pow(10,p0-addlogs(p0,p1)); Flist[f].PQ = (int)(10*(addlogs(p0,p1)-p1)); }
	   else { tag = '1'; prob = pow(10,p1-addlogs(p0,p1)); Flist[f].PQ = (int)(10*(addlogs(p0,p1)-p0)); } 
	   if (prob >=0.9 && offset >=0) { 
	   component = snpfrag[offset].component; // unphased SNPs should be ignored
	   // print PS,hap(0|1),probability followed by original fragment copy
	   if (component < 0) { continue; fprintf(stderr,"error \n"); } 
	   Flist[f].PS = snpfrag[component].position; Flist[f].HP = tag; 
	   fprintf(OUTFILE,"%s %d %c %d\n",Flist[f].id,Flist[f].PS,Flist[f].HP,Flist[f].PQ);
	   /*
	   fprintf(OUTFILE,"%d %c %d %f %f ",Flist[f].PS,Flist[f].HP,Flist[f].PQ,p0,p1);
	   fprintf(OUTFILE,"%d %s ",Flist[f].blocks,Flist[f].id); 
	   for (j = 0; j < Flist[f].blocks; j++) {
		fprintf(OUTFILE,"%d ",Flist[f].list[j].offset); 
		for (k = 0; k < Flist[f].list[j].len; k++) fprintf(OUTFILE,"%c",Flist[f].list[j].hap[k]); 
		fprintf(OUTFILE," ");
	   } 
	   for (j = 0; j < Flist[f].blocks; j++) { 
		for (k = 0; k < Flist[f].list[j].len; k++) fprintf(OUTFILE,"%c",Flist[f].list[j].qv[k]); 
	   }
	   fprintf(OUTFILE,"\n");
           */
	   }
    }
    fclose(OUTFILE);
}

