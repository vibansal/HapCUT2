#include "hapcontig.h"
#include "common.h"
extern int DISCRETE_PRUNING;
extern int ERROR_ANALYSIS_MODE;
extern int SPLIT_BLOCKS;
extern int SKIP_PRUNE;
extern float THRESHOLD;
extern int HIC;
extern int VERBOSE;

// populate the connected component data structure only for non-trivial connected components, at least 2 variants
void generate_contigs(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int components, struct BLOCK* clist) {
    // bcomp maps the SNP to the component number in clist since components << snps
    // note that the previous invariant about the first SNP (ordered by position) being the root of the component is no longer true !!
    // should we still require it ??
    int i = 0, component = 0;
    for (i = 0; i < snps; i++) snpfrag[i].bcomp = -1;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1) // root node of component
        {
            snpfrag[i].bcomp = component;
            clist[component].slist = calloc(sizeof (int), snpfrag[i].csize);
            clist[component].phased = 0;
            component++;
        }
    }
    //fprintf_time(stderr,"non-trivial components in graph %d \n",components);
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component < 0) continue; // to allow for initialization to -1 in other code feb 15 2013
        if (snpfrag[i].csize <= 1 && snpfrag[i].component == i) continue; // ignore singletons that are not connected to other variants

        if (snpfrag[i].component != i) snpfrag[i].bcomp = snpfrag[snpfrag[i].component].bcomp;
        if (snpfrag[i].bcomp < 0) continue;
        component = snpfrag[i].bcomp;
        if (clist[component].phased == 0) clist[component].offset = i;
        clist[component].slist[clist[component].phased] = i;
        clist[component].phased++;
        clist[component].lastvar = i;
    }
    for (i = 0; i < components; i++) clist[i].length = clist[i].lastvar - clist[i].offset + 1;
    for (i = 0; i < components; i++) clist[i].frags = 0;
    for (i = 0; i < fragments; i++) {
        if (snpfrag[Flist[i].list[0].offset].bcomp < 0)continue; // ignore fragments that cover singleton vertices
        clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++;
    }
    for (i = 0; i < components; i++) clist[i].flist = calloc(clist[i].frags,sizeof (int));
    for (i = 0; i < components; i++) clist[i].frags = 0;
    for (i = 0; i < fragments; i++) {
        if (snpfrag[Flist[i].list[0].offset].bcomp < 0)continue;
        clist[snpfrag[Flist[i].list[0].offset].bcomp].flist[clist[snpfrag[Flist[i].list[0].offset].bcomp].frags] = i;
        clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++;
    }
    for (i = 0; i < components; i++) {
        if (VERBOSE) fprintf_time(stdout, "comp %d first %d last %d phased %d fragments %d \n", i, clist[i].offset, clist[i].lastvar, clist[i].phased, clist[i].frags);
    }
}

// THIS FUNCTION PRINTS THE CURRENT HAPLOTYPE ASSEMBLY in a new file block by block
int print_contigs(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* outfile) {
    // print a new file containing one block phasing and the corresponding fragments
    int i = 0, t = 0, k = 0, span = 0;
    char c=0, c1=0, c2=0;
    FILE* fp;
    fp = fopen(outfile, "w");
    float logthresh=log10(THRESHOLD);

    int flip_alleles=0;
   // switch_haplotype_order(clist,blocks,h1); // make sure that haplotype in column '2' has '0' allele for first variant in each block

    for (i = 0; i < blocks; i++) {
        span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
        fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
        fprintf(fp, "SPAN: %d fragments %d\n", span, clist[i].frags);
     
	if (h1[clist[i].slist[0]] == '0') flip_alleles = 0;
        else flip_alleles = 1;
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
                &&((snpfrag[t].post_hap < logthresh && !DISCRETE_PRUNING)
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
		if (flip_alleles ==0)  fprintf(fp, "%d\t%c\t%c\t", t + 1, c1, c2); // two alleles that are phased in VCF like format
		else  fprintf(fp, "%d\t%c\t%c\t", t + 1, c2, c1); // two alleles that are phased in VCF like format
            } else {
                if (flip_alleles ==0) fprintf(fp, "%d\t%c\t%c\t", t + 1, h1[t], c);
                else fprintf(fp, "%d\t%c\t%c\t", t + 1,c, h1[t]);
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

            // generate the string for the snp confidence, single mismatch
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
	   fprintf(fp,"\n");  

        }
        if (i < blocks - 1) fprintf(fp, "******** \n");
    }
    fclose(fp);
    return 0;
}

int cmpint (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

float calculate_N50(struct BLOCK* clist, int blocks, struct SNPfrags* snpfrag, char* HAP1) 
{
    int i=0;
    int span=0,totalspan=0;
    int* lengths = calloc(sizeof(int),blocks);

    for (i = 0; i < blocks; i++) {
        span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
        lengths[i] = span; totalspan += span;
        //fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
        //for (k = 0; k < clist[i].phased; k++)  t = clist[i].slist[k];
    }
    qsort(lengths,blocks,sizeof(int),cmpint);
    //int N50block = 0;
    int sumlengths=0;
    float N50length=0;
    for (i=0;i<blocks;i++)
    {
	span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
	sumlengths += span;
	if (sumlengths*2 > totalspan) 
	{
		//N50block = i; 
		N50length = span;
		break;
	}
    }   
    N50length /=1000;
   // fprintf_time(stderr,"N50 haplotype length is %0.2f kilobases \n",N50length);
    free(lengths);
    return N50length;
}

