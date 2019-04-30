#include "common.h"
#include <zlib.h>

extern int DISCRETE_PRUNING;
extern int ERROR_ANALYSIS_MODE;
extern int SPLIT_BLOCKS;
extern int SKIP_PRUNE;
extern float THRESHOLD;
extern int HIC;

void print_hapcut_options() {
    fprintf(stdout, "\nHapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies\n\n");
    fprintf(stdout, "USAGE : ./HAPCUT2 --fragments fragment_file --VCF variantcalls.vcf --output haplotype_output_file\n\n");
    fprintf(stderr, "Basic Options:\n");
    fprintf(stderr, "--fragments, --f <FILENAME>:        file with haplotype-informative reads generated using the extracthairs program\n");
    fprintf(stderr, "--VCF <FILENAME>:                   variant file in VCF format (use EXACT SAME file that was used for the extracthairs program)\n");
    fprintf(stderr, "--output, --o <OUTFILE> :          file to which phased haplotype segments/blocks will be output \n");
    fprintf(stderr, "--outvcf <0/1> :                    output phased variants to VCF file (<OUTFILE>.phased.vcf), default: 0 \n");
    fprintf(stderr, "--converge, --c <int>:              cut off iterations (global or maxcut) after this many iterations with no improvement. default: 5\n");
    fprintf(stderr, "--verbose, --v <0/1>:               verbose mode: print extra information to stdout and stderr. default: 0\n");

    fprintf(stderr, "\nRead Technology Options:\n");
    fprintf(stderr, "--hic <0/1> :                       increases accuracy on Hi-C data; models h-trans errors directly from the data. default: 0\n");
    fprintf(stderr, "--hic_htrans_file, --hf <FILENAME>  optional tab-delimited input file where second column specifies h-trans error probabilities for insert size bins 0-50Kb, 50Kb-100Kb, etc.\n");
    fprintf(stderr, "--qv_offset, --qo <33/48/64> :      quality value offset for base quality scores, default: 33 (use same value as for extracthairs)\n");
    fprintf(stderr, "--long_reads, --lr <0/1> :          reduces memory when phasing long read data with many SNPs per read. default: automatic.\n");

    fprintf(stderr, "\nHaplotype Post-Processing Options:\n");
    fprintf(stderr, "--threshold, --t <float>:           PHRED SCALED threshold for pruning low-confidence SNPs (range 0-100, larger values prune more.). default: 6.98\n");
    fprintf(stderr, "--skip_prune, --sp <0/1>:           skip default likelihood pruning step (prune SNPs after the fact using column 11 of the output). default: 0\n");
    //fprintf(stderr, "--split_blocks, --sb <0/1>:         split blocks using simple likelihood score to reduce switch errors. default: 0\n");
    //fprintf(stderr, "--split_threshold, --st <float>:    PHRED SCALED threshold for splitting blocks (range 0-100, larger values split more). default: 0\n");
    fprintf(stderr, "--call_homozygous, --ch <0/1>:      call positions as homozygous if they appear to be false heterozygotes. default: 0\n");
    fprintf(stderr, "--discrete_pruning, --dp <0/1>:     use discrete heuristic to prune SNPs. default: 0\n");
    fprintf(stderr, "--error_analysis_mode, --ea <0/1>:  compute switch confidence scores and print to haplotype file but don't split blocks or prune. default: 0\n");

    fprintf(stderr, "\nAdvanced Options:\n");
    fprintf(stderr, "--new_format, --nf <0/1>:           use new Hi-C fragment matrix file format (but don't do h-trans error modeling). default: 0\n");
    fprintf(stderr, "--max_iter, --mi <int> :            maximum number of global iterations. Preferable to tweak --converge option instead. default: 10000\n");
    fprintf(stderr, "--maxcut_iter, --mc <int> :         maximum number of max-likelihood-cut iterations. Preferable to tweak --converge option instead. default: 10000\n");
    fprintf(stderr, "--htrans_read_lowbound, --hrl <int> with --hic on, h-trans probability estimation will require this many matepairs per window. default: 500\n");
    fprintf(stderr, "--htrans_max_window, --hmw <int>    with --hic on, the insert-size window for h-trans probability estimation will not expand larger than this many basepairs. default: 4000000\n");

    fprintf(stderr, "\n\nHi-C-specific Notes:\n");
    fprintf(stderr, "  (1) When running extractHAIRS, must use --hic 1 option to create a fragment matrix in the new Hi-C format.\n");
    fprintf(stderr, "  (2) When running HapCUT2, use --hic 1 if h-trans probabilities are unknown. Use --hic_htrans_file if they are known\n");
    fprintf(stderr, "  (3) Using --hic_htrans_file is faster than --hic and may yield better results at low read coverage (<30x).\n");
    fprintf(stderr, "\n");

}

// THIS FUNCTION PRINTS THE CURRENT HAPLOTYPE ASSEMBLY in a new file block by block
int print_hapfile(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fname, int score, char* outfile) {
    // print a new file containing one block phasing and the corresponding fragments
    int i = 0, t = 0, k = 0, span = 0;
    char c=0, c1=0, c2=0;
    //char fn[200]; sprintf(fn,"%s-%d.phase",fname,score);
    FILE* fp;
    //int edgelist[4096];
    fp = fopen(outfile, "w");

    for (i = 0; i < blocks; i++) {
        span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
        fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
        fprintf(fp, "SPAN: %d fragments %d\n", span, clist[i].frags);
        for (k = 0; k < clist[i].phased; k++) {

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

void print_VCF_header(FILE* outfile)
{
        //fprintf(outfile,"##source=HapCUT2 phased haplotype blocks\n");
        //fprintf(outfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        fprintf(outfile,"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ID of Phase Set for Variant\">\n");
        fprintf(outfile,"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype\">\n");
        //fprintf(outfile,"##FORMAT=<ID=JQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability of a phasing switch error in gap prior to this variant\">\n");
        //fprintf(outfile,"##FORMAT=<ID=BX,Number=.,Type=String,Description=\"Barcodes for this variant\">\n");
        fprintf(outfile,"##FORMAT=<ID=PD,Number=1,Type=Integer,Description=\" phased Read Depth\">\n");
        //fprintf(outfile,"##commandline=\"....\"";

}

char* concatStrings(char** var_list,int n,char sep)
{
	// concatenate string GT + ':' + DP + ':' from an array
	int i=0,l=0,j=0,k=0;
	for (i=0;i<n;i++) l += strlen(var_list[i]) +1;
	char* concat = calloc(l+1,sizeof(char));
	for (i=0;i<n;i++)
	{
		for (j=0;j<strlen(var_list[i]);j++) concat[k++] = var_list[i][j];
		concat[k++] = sep;
	}
	concat[k] = '\0';
	return concat;
}

// split a string using a single separator, '\n' and '\0' are also delimitors at end of string
int splitString(char* input,char sep,char** var_list)
{
	int n=0,i=0,s=0,j=0;
	while (1)
	{
		if (input[i] == sep || input[i] == '\n' || input[i] == '\0')
		{
			if (i-s > 0) {
				var_list[n] = malloc(i-s+1);
				for (j = s; j < i; j++) var_list[n][j - s] = input[j];

        var_list[n][j-s] = '\0';
				n +=1;
			}
			s = i+1; // start of next string
		}
		if (input[i] =='\n' || input[i] == '\0') break;
		i++;
	}
	return n;
}


// print VCF header copied from original VCF file
// PS tag for each variant based on position of first SNP in component
// 1/2 genotypes or 0/2 genotypes need to be handled properly
int output_vcf(char* vcffile, struct SNPfrags* snpfrag, int snps,char* H1,struct fragment* Flist, int fragments,char* outvcf,int BARCODES, char* ONLY_CHROM) {
    char* buffer= malloc(500000);
    FILE* out = fopen(outvcf,"w");
    gzFile fp = gzopen(vcffile, "r");
    int var=0;

    char h1='0',h2='0';
    int phased=0,component=0,PS=0,PQ=100,PD=0;
    char PGT[4];

    char** var_list = calloc(1024,sizeof(char*)); char** info_list = calloc(1024,sizeof(char*)); char** geno_list = calloc(1024,sizeof(char*));
    char* newC9 = malloc(4096); char* newC10 = malloc(4096); // column 9 and 10 of new VCF file
    char* barcodes_0 = malloc(1024*64); char* barcodes_1 = malloc(1024*64); char* barcodes_2 = malloc(1024*64);
    char* temp = malloc(10024); char* seek;
    int i=0,vn=0,in=0,gn=0,j=0,k=0;
    float snp_conf_fl;
    int b0=0,b1=0;

    int lines=0,formatline=0,infolines=0;
    while (gzgets(fp, buffer, 500000)) {
	lines +=1;
        if (buffer[0] == '#') {

        if (strcmp(ONLY_CHROM,"None")!=0) {
          continue;
        }

		if (lines ==1 && strncmp(buffer,"##fileformat",12) == 0) formatline =1;
		if (lines >= 1 && formatline ==0)  // no header present
		{
       			fprintf(out,"##fileformat=VCFv4.0\n");	formatline = 1;
		}
		if (strncmp(buffer,"##INFO=<ID",10) ==0)
		{
			if (infolines ==0) fprintf(out,"##INFO=<ID=hapcut2,Number=1,Type=Integer,Description=\"phased by HapCUT2 or not\">\n");
			infolines++;
		}

		if (strncmp(buffer,"#CHROM",6) ==0) print_VCF_header(out);
		if (strstr(buffer,"FORMAT=<ID=PS") == NULL && strstr(buffer,"FORMAT=<ID=PD") == NULL && strstr(buffer,"FORMAT=<ID=PQ") == NULL) fprintf(out,"%s",buffer);
		continue;
	}


  if (strcmp(ONLY_CHROM,"None")!=0 && strcmp(ONLY_CHROM,snpfrag[var].chromosome)!=0) {
      var++;
      continue;
  }

	vn= splitString(buffer,'\t',var_list); // split the line
	in = splitString(var_list[8],':',info_list); // split INFO column (9)
	gn = splitString(var_list[9],':',geno_list); // split GENOTYPE column (10)

	//thirdallele =0;
	//if (geno_list[0][0] == '2' || geno_list[0][2] =='2') thirdallele = 1; // need to fix phased genotype, hapcut only uses two alleles labeled 0 & 1
	if (H1[var] == '0') { h1 = '0'; h2 = '1'; }
	else if (H1[var] == '1') { h1 = '1'; h2 = '0'; }
	else {h1 = '.'; h2 = '.'; } // unphased by hapcut2


	component = snpfrag[var].component;
	if (H1[var] == '-' || component < 0 || snpfrag[component].csize < 2) phased=0;
	else phased = 1;

	if (phased ==0)
	{
		sprintf(PGT,"%c/%c",geno_list[0][0],geno_list[0][2]);
	}
	else
	{
		PS = snpfrag[component].position;
		snp_conf_fl = -10.0 * subtractlogs(0,snpfrag[var].post_hap);
		if (snp_conf_fl > 100.0)     snp_conf_fl = 100.0;
		PQ = (int)(snp_conf_fl+0.5); // round to nearest int
		PD = snpfrag[var].frags; // should we ignore singleton fragments
		sprintf(PGT,"%c|%c",h1,h2);

		if ( (!ERROR_ANALYSIS_MODE) && (!SKIP_PRUNE) && (snpfrag[var].post_hap < log10(THRESHOLD) && !DISCRETE_PRUNING ))
		{
			h1='.'; h2 = '.'; sprintf(PGT,"./."); phased =2;
		}
	}

	strcpy(newC9,"GT"); strcpy(newC10,PGT);
	for (i=0;i<gn;i++)
	{
		// preserve original genotype as 'OGT' field
		if ( ((strcmp(info_list[i],"AD") ==0 || strcmp(info_list[i],"DP") ==0 || strcmp(info_list[i],"GQ") ==0) ) && phased >=1  )
		{
			strcat(newC9,":"); strcat(newC9,info_list[i]); strcat(newC10,":"); strcat(newC10,geno_list[i]);
		}
		if ( strcmp(info_list[i],"PQ") !=0 && strcmp(info_list[i],"PS") !=0 && strcmp(info_list[i],"GT") !=0 && strcmp(info_list[i],"PD") !=0 && phased ==0  )
		{
			strcat(newC9,":"); strcat(newC9,info_list[i]); strcat(newC10,":"); strcat(newC10,geno_list[i]);
		}
	}
	if (phased ==1)
	{
		//fprintf(out,"var %d PGLL %f %f %f %f %f\n",var,snpfrag[var].PGLL[0],snpfrag[var].PGLL[1],snpfrag[var].PGLL[2],snpfrag[var].PGLL[3],snpfrag[var].PGLL[4]);
		b0=0;b1=0;
		if (BARCODES ==1) {
			for (i=0;i<snpfrag[var].frags;i++)
			{
				if (Flist[snpfrag[var].flist[i]].HP == '0')
				{
					if (b0 ==0) strcpy(barcodes_0,"");
					else strcat(barcodes_0,";");
					j=0; k=0; seek = Flist[snpfrag[var].flist[i]].id;
					while (seek[j] != '\0')
					{
						if (seek[j++] == ':') k +=1;
						if (k ==2) break;
					}
					strcat(barcodes_0,Flist[snpfrag[var].flist[i]].id+j);
					b0 +=1;
				}
				if (Flist[snpfrag[var].flist[i]].HP == '1')
				{
					if (b1 ==0) strcpy(barcodes_1,"");
					else strcat(barcodes_1,";");
					j=0; k=0; seek = Flist[snpfrag[var].flist[i]].id;
					while (seek[j] != '\0')
					{
						if (seek[j++] == ':') k +=1;
						if (k ==2) break;
					}
					strcat(barcodes_1,Flist[snpfrag[var].flist[i]].id+j);
					b1 +=1;
				}
			}
		}

		// print PS,BQ,PQ,PD tags
		strcat(newC9,":PS"); strcat(newC10,":"); sprintf(temp,"%d",PS); strcat(newC10,temp);
		strcat(newC9,":PQ"); strcat(newC10,":"); sprintf(temp,"%d",PQ); strcat(newC10,temp);
		strcat(newC9,":PD"); strcat(newC10,":"); sprintf(temp,"%d",PD); strcat(newC10,temp);
		if(BARCODES ==1) {  strcat(newC9,":BX"); strcat(newC10,":"); strcat(newC10,barcodes_0); strcat(newC10,","); strcat(newC10,barcodes_1); }
		// MEC score for each variant, phased genotype likelihoods
	}

	// HAPCUT = 0/1/2 (pruned, not phased, not considered...)
	for (i=0;i<7;i++) fprintf(out,"%s\t",var_list[i]);
	if (strcmp(var_list[7],".") ==0) sprintf(temp,"hapcut2=%d",phased);
	else sprintf(temp,"%s;hapcut2=%d",var_list[7],phased);
	fprintf(out,"%s\t%s\t%s\n",temp,newC9,newC10);

	//for (i=0;i<in;i++) fprintf(stdout,"%s ",info_list[i]); fprintf(stdout,"\n");
	//for (i=0;i<gn;i++) fprintf(stdout,"%s(%s) ",geno_list[i],info_list[i]); fprintf(stdout,"\n\n");
	//for (i=0;i<vn;i++) fprintf(stdout,"%s\n",var_list[i]); fprintf(stdout,"## var %d\n\n",var);

        var++;
	for (i=0;i<vn;i++) free(var_list[i]);
	for (i=0;i<gn;i++) free(geno_list[i]);
	for (i=0;i<in;i++) free(info_list[i]);
    }
    gzclose(fp); fclose(out); free(buffer); free(barcodes_0); free(barcodes_1); free(barcodes_2);
    return 1;
}
