
void print_VCF_header(FILE* outfile)
{
        //fprintf(outfile,"##source=HapCUT2 phased haplotype blocks\n");
        //fprintf(outfile,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        fprintf(outfile,"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ID of Phase Set for Variant\">\n");
        fprintf(outfile,"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype\">\n");
        //fprintf(outfile,"##FORMAT=<ID=JQ,Number=1,Type=Integer,Description=\"Phred QV indicating probability of a phasing switch error in gap prior to this variant\">\n");
        //fprintf(outfile,"##FORMAT=<ID=BX,Number=.,Type=String,Description=\"Barcodes for this variant\">\n");
        fprintf(outfile,"##FORMAT=<ID=PD,Number=1,Type=Integer,Description=\"phased Read Depth\">\n");
        //fprintf(outfile,"##commandline=\"....\"";

}

// print VCF header copied from original VCF file 
// PS tag for each variant based on position of first SNP in component 
// 1/2 genotypes or 0/2 genotypes need to be handled properly 
int output_vcf(char* vcffile, struct SNPfrags* snpfrag, int snps,char* H1,struct fragment* Flist, int fragments,char* outvcf,int BARCODES) {
    char* buffer= malloc(1000000);
    FILE* out = fopen(outvcf,"w"); 
    FILE* fp = fopen(vcffile, "r");
    int var=0;
	
    char h1='0',h2='0';
    int phased=0,component=0,PS=0,PQ=100,PD=0; 
    char PGT[4];

    char** var_list = calloc(1024,sizeof(char*)); 
    char** info_list = calloc(1024,sizeof(char*)); 
    char** geno_list = calloc(1024,sizeof(char*));
    char* newC9 = malloc(4096); char* newC10 = malloc(4096); // column 9 and 10 of new VCF file 
    char* temp = malloc(10024); 
    //char* seek;
    int i=0,vn=0,in=0,gn=0;
    float snp_conf_fl;	
    //int b0=0,b=0;
 
    int lines=0,formatline=0,infolines=0;
    while (fgets(buffer, 1000000, fp)) { // read input VCF
	lines +=1; 
        if (buffer[0] == '#') {
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
		if (strncmp(buffer,"#CHROM",6) ==0) print_VCF_header(out); // print additional lines in the header of VCF
		if (strstr(buffer,"FORMAT=<ID=PS") == NULL && strstr(buffer,"FORMAT=<ID=PD") == NULL && strstr(buffer,"FORMAT=<ID=PQ") == NULL) fprintf(out,"%s",buffer);
		continue;
	}

	vn= splitString(buffer,'\t',var_list); // split the line 
	in = splitString(var_list[8],':',info_list); // split INFO column (9)
	gn = splitString(var_list[9],':',geno_list); // split GENOTYPE column (10), assumption is single sample VCF
	
	if (geno_list[0][0] == '2' || geno_list[0][2] =='2')  // tri-allelic genotype
	{
		if (H1[var] == '0') { h1 = geno_list[0][0]; h2 = geno_list[0][2]; } 
		else if (H1[var] == '1') { h1 = geno_list[0][2]; h2 = geno_list[0][0]; } 
	}
	else // standard bi-allelic genotype
	{
		if (H1[var] == '0') { h1 = '0'; h2 = '1'; } 	
		else if (H1[var] == '1') { h1 = '1'; h2 = '0'; } 
		else {h1 = '.'; h2 = '.'; } // unphased by hapcut2	
	}
            
	component = snpfrag[var].component; 
	if (H1[var] == '-' || component < 0 || snpfrag[component].csize < 2) phased=0; 
	else phased = 1; 

	if (phased ==0) sprintf(PGT,"%c/%c",geno_list[0][0],geno_list[0][2]);
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
		//b0=0;b1=0;
		// print PS,BQ,PQ,PD tags
		strcat(newC9,":PQ"); strcat(newC10,":"); sprintf(temp,"%d",PQ); strcat(newC10,temp);
		strcat(newC9,":PD"); strcat(newC10,":"); sprintf(temp,"%d",PD); strcat(newC10,temp);
		strcat(newC9,":PS"); strcat(newC10,":"); sprintf(temp,"%d",PS); strcat(newC10,temp);
		//if(BARCODES ==1) {  strcat(newC9,":BX"); strcat(newC10,":"); strcat(newC10,barcodes_0); strcat(newC10,","); strcat(newC10,barcodes_1); } 
		// MEC score for each variant, phased genotype likelihoods 
	}
	else // unphased variant, add PS tag to default genotype and set to missing
	{
		strcat(newC9,":PS"); strcat(newC10,":."); 
	}
	
	for (i=0;i<8;i++) fprintf(out,"%s\t",var_list[i]);  // first 8 columns as is
         
	// HAPCUT = 0/1/2 (pruned, not phased, not considered...)
	//if (strcmp(var_list[7],".") ==0) sprintf(temp,"hapcut2=%d",phased); 
	//else sprintf(temp,"%s;hapcut2=%d",var_list[7],phased); 
	fprintf(out,"%s\t%s\n",newC9,newC10);
	
	//for (i=0;i<in;i++) fprintf(stdout,"%s ",info_list[i]); fprintf(stdout,"\n");
	//for (i=0;i<gn;i++) fprintf(stdout,"%s(%s) ",geno_list[i],info_list[i]); fprintf(stdout,"\n\n");
	//for (i=0;i<vn;i++) fprintf(stdout,"%s\n",var_list[i]); fprintf(stdout,"## var %d\n\n",var); 
        var++;
	for (i=0;i<vn;i++) free(var_list[i]);
	for (i=0;i<gn;i++) free(geno_list[i]);
	for (i=0;i<in;i++) free(info_list[i]);
    }
    fclose(fp); 
    fclose(out); 
    free(temp); free(newC9); free(newC10);
    free(buffer); 
    return 1;
}


// CODE to print barcodes for linked-reads, currently not used
/*
    char* barcodes_0 = malloc(1024*64); char* barcodes_1 = malloc(1024*64); char* barcodes_2 = malloc(1024*64);
    free(barcodes_0); free(barcodes_1); free(barcodes_2);
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

*/


	/* 
	htsFile *vcfpointer = NULL;
        bcf_hdr_t *vcf_header = NULL;
        bcf1_t *vcf_record = bcf_init();
        vcfpointer = bcf_open(vcffile, "r");
        if(vcfpointer == NULL) { fprintf(stderr,"vcf file not available \n"); return -1; 
        }
        vcf_header = bcf_hdr_read(vcfpointer);
        if(vcf_header == NULL) {fprintf(stderr,"vcf header not available \n"); return -1;
        }
	int* gq=NULL; int ngq_arr=0;

        while(bcf_read(vcfpointer, vcf_header, vcf_record) == 0) {
	var++;

	int ngq = bcf_get_format_int32(vcf_header, vcf_record, "GQ", &gq, &ngq_arr);

	// does not work if header does not have chromosome lengths [W::vcf_parse] contig '20' is not defined in the header
	fprintf(stdout," %s ",bcf_hdr_id2name(vcf_header, vcf_record->rid));
	fprintf(stdout,"%d %d\n",vcf_record->pos,*gq);

        // if variant is phased, need to update genotype and other fields
        // if vairant is unphased, need to set genotype to '0/1' 
        // what about homozygous genotypes 

        }

        bcf_hdr_destroy(vcf_header);
        bcf_destroy(vcf_record); 
        bcf_close(vcfpointer);

	*/
