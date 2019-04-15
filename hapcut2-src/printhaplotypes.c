#include "common.h"
extern int DISCRETE_PRUNING;
extern int ERROR_ANALYSIS_MODE;
extern int SPLIT_BLOCKS;
extern int SKIP_PRUNE;
extern float THRESHOLD;
extern int HIC;

// THIS FUNCTION PRINTS THE CURRENT HAPLOTYPE ASSEMBLY in a new file block by block
int print_hapfile(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fname, int score, char* outfile) {
    // print a new file containing one block phasing and the corresponding fragments
    int i = 0, t = 0, k = 0, span = 0;
    char c=0, c1=0, c2=0;
    //char fn[200]; sprintf(fn,"%s-%d.phase",fname,score);
    FILE* fp;
    int edgelist[4096];
    fp = fopen(outfile, "w");

    for (i = 0; i < blocks; i++) {
        span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
        fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
        fprintf(fp, "SPAN: %d fragments %d\n", span, clist[i].frags);
        for (k = 0; k < clist[i].phased; k++) { // only variants linked in component are phased 

            t = clist[i].slist[k];
            if (h1[t] =='0')
                c= '1';
            else if (h1[t] =='1')
                c = '0';
            else c = h1[t];
            // print this line to keep consistency with old format
            // if SNP was pruned then print '-'s
            if ((!ERROR_ANALYSIS_MODE)&&(!SKIP_PRUNE)
                &&((snpfrag[t].post_hap < log10(THRESHOLD) && !DISCRETE_PRUNING)
                ||(snpfrag[t].pruned_discrete_heuristic && DISCRETE_PRUNING))){
                fprintf(fp, "%d\t-\t-\t", t + 1);
            }else if (snpfrag[t].genotypes[0] == '0' && snpfrag[t].genotypes[2] == '0'){
                fprintf(fp, "%d\t0\t0\t", t + 1);   // homozygous 00
            }else if (snpfrag[t].genotypes[0] == '1' && snpfrag[t].genotypes[2] == '1'){
                fprintf(fp, "%d\t1\t1\t", t + 1);   // homozygous 11
            }else if (snpfrag[t].genotypes[0] == '2' || snpfrag[t].genotypes[2] == '2') {

                if (h1[t] == '0') {
                    c1 = snpfrag[t].genotypes[0];
                    c2 = snpfrag[t].genotypes[2];
                } else if (h1[t] == '1') {
                    c2 = snpfrag[t].genotypes[0];
                    c1 = snpfrag[t].genotypes[2];
                }
                fprintf(fp, "%d\t%c\t%c\t", t + 1, c1, c2); // two alleles that are phased in VCF like format
            } else {
                fprintf(fp, "%d\t%c\t%c\t", t + 1, h1[t], c);
            }

            // generate the string for the switch confidence
            char switch_conf[100];
            if ((SPLIT_BLOCKS || ERROR_ANALYSIS_MODE)&&(!HIC)){
                float switch_conf_fl = -10.0 * subtractlogs(0,snpfrag[t].post_notsw);
                if (switch_conf_fl > 100.0){
                    switch_conf_fl = 100.0;
                }
                if (!(switch_conf_fl >= 0.0 && switch_conf_fl <= 100.0)){
                    fprintf(stderr, "Invalid switch confidence score\n");
                    exit(1);
                }
                sprintf(switch_conf,"%0.2f",switch_conf_fl);
            }else{
                strcpy(switch_conf,".");
            }

            // generate the string for the snp confidence
            char snp_conf[100];
            char discrete_conf[100];
            if (SKIP_PRUNE){
                strcpy(snp_conf, ".");
                strcpy(discrete_conf, ".");
            }else{
                float snp_conf_fl = -10.0 * subtractlogs(0,snpfrag[t].post_hap);
                if (snp_conf_fl > 100.0){
                    snp_conf_fl = 100.0;
                }
                if (!(snp_conf_fl >= 0.0 && snp_conf_fl <= 100.0)){
                    fprintf(stderr, "Invalid SNV confidence score\n");
                    exit(1);
                }
                sprintf(snp_conf,"%0.2f",snp_conf_fl);
                sprintf(discrete_conf,"%d",snpfrag[t].pruned_discrete_heuristic);
            }

            fprintf(fp, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d", snpfrag[t].chromosome, snpfrag[t].position, snpfrag[t].allele0, snpfrag[t].allele1, snpfrag[t].genotypes, discrete_conf, switch_conf, snp_conf,snpfrag[t].frags);

	   //if (HiC && snpfrag[t].frags > 0 && snpfrag[t].tedges < 20 )
	 //  {
	//	int t1=0,f1=0,l1=0,e1=0;
		/*
		for (f1=0;f1< snpfrag[t].frags;f1++); 
		{
			for (t1=0;t1<Flist[f1].blocks;t1++) 
			{
			    for (l1=0;l1<Flist[f1].len[t1];l1++) edgelist[e1++] = Flist[f1].offset + l1; 
			}
		}*/
	  // }
	   fprintf(fp,"\n");  

        }
        if (i < blocks - 1) fprintf(fp, "******** \n");
    }
    fclose(fp);
    return 0;
}

