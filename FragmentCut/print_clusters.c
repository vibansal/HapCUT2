
// functions used for parsing long-fragment (dilution-pool) sequence data 

int cmp_fosmid(const void *a,const void *b)
{
	const struct FOSMID* fa = (const struct FOSMID*)a;  const struct FOSMID* fb = (const struct FOSMID*)b; 
	return fa->firstblock- fb->firstblock; 
}

// find mate in readlist and also identify PCR duplicates partially
void find_matepair(struct alignedread** readlist,int s,int e,int max_insert_size)
{
        int i=0,j=0;
        int prevstart = 0, prevIS = 0, duplicates = 0;
        for (i=s;i<e;i++)
        {
                if (readlist[i]->IS <0) continue;
                else if (readlist[i]->IS >= 0) // insert size is positive 
                {
                        if (readlist[i]->position == prevstart  && readlist[i]->IS == prevIS)
                        {
                                readlist[i]->flag |= 1024; //
                                duplicates++;
                        }
                        prevstart = readlist[i]->position; prevIS = readlist[i]->IS;

                }
                readlist[i]->mateindex = -1;
                if (readlist[i]->IS > 0 && readlist[i]->IS < max_insert_size)
                {
                        // search forward in list for mate
                        for (j = i+1; j < e && readlist[j]->position <= readlist[i]->mateposition;j++)
                        {
                                if (readlist[j]->position == readlist[i]->mateposition && readlist[j]->mateposition == readlist[i]->position && strcmp(readlist[i]->readid, readlist[j]->readid) ==0)
                                { readlist[i]->mateindex =  j; readlist[j]->mateindex = i; break; }
                        }

                }
        }
}


// compact print
int print_fragment_vars1(FRAGMENT* fragment,VARIANT* varlist)
{
	int i=0,j=0,prev = -1,hets=0,uniquevar=1,reads=0;
	int ct[4] = {0,0,0,0}; // counts and quality for each allele 
	for (i=0;i<fragment->variants+1;i++) // some variants can be duplicate 
	{
		if (i ==0) fprintf(stdout,"offset:%d | ",fragment->alist[i].varid);
		if (i == fragment->variants || (prev != fragment->alist[i].varid && prev >=0)) 
		{
			fprintf(stdout,"%d ",varlist[prev].position);
			if (ct[0] > 0) { for (j=0;j<ct[0];j++) fprintf(stdout,""); fprintf(stdout,"0:c=%d:%0.1f ",ct[0],(float)ct[1]/ct[0]); } 
			if (ct[2] > 0) { for (j=0;j<ct[2];j++) fprintf(stdout,""); fprintf(stdout,"1:c=%d:%0.1f ",ct[2],(float)ct[3]/ct[2]); }
			if (ct[0] > 0 && ct[2] > 0) fprintf(stdout,"HET | "); else fprintf(stdout,"| "); 
			if (ct[0] > 0 && ct[2] > 0)
			{
				if ((float)ct[0] < 0.95*(float)(ct[2]+ct[0]) && (float)ct[0] > 0.05*(float)(ct[0]+ct[2]) ) hets++;
				else if (ct[0] >= 2 && ct[2] >= 2) hets++;
			} 
			ct[0] = 0; ct[1] = 0; ct[2] = 0; ct[3] = 0; uniquevar +=1; 
		
		} 
		if (i==fragment->variants) break;
		if (fragment->alist[i].allele == '0' && fragment->alist[i].qv-QVoffset >= MINQ) { ct[0] +=1; ct[1] += fragment->alist[i].qv-QVoffset; } 
		if (fragment->alist[i].allele == '1' && fragment->alist[i].qv-QVoffset >= MINQ) { ct[2] +=1; ct[3] += fragment->alist[i].qv-QVoffset; } 

		prev = fragment->alist[i].varid; 
	}
	fprintf(stdout," | hets: %d/%d\n",hets,uniquevar);
	return hets;
}

// not used currently
void print_fragment_vars(FRAGMENT* fragment,VARIANT* varlist)
{
	fprintf(stdout,"---------variants------------\n");
	int i=0,j=0,prev = -1,hets=0,prevhet = -1,uniquevar=1;
	for (i=0;i<fragment->variants;i++) // some variants can be duplicate 
	{
		j = fragment->alist[i].varid;
		if (prev != j) 
		{
			if (prev >= 0) fprintf(stdout,"\n");
			fprintf(stdout,"%d:%d %s/%s %c:%d | ",j,varlist[j].position,varlist[j].allele1,varlist[j].allele2,fragment->alist[i].allele,fragment->alist[i].qv-QVoffset);
			uniquevar++;
		}
		else // same variant multiple reads
		{
			if (fragment->alist[i].allele != fragment->alist[i-1].allele) 
			{
				fprintf(stdout,"%c:%d:HET | ",fragment->alist[i].allele,fragment->alist[i].qv-QVoffset); 
				if (prevhet != j) hets +=1; prevhet = j; 
			}
			else fprintf(stdout,"%c:%d | ",fragment->alist[i].allele,fragment->alist[i].qv-QVoffset);
		}
		prev = j; 
	}
	fprintf(stdout,"\n----- hetvariants:%d/%d------ \n",hets,uniquevar);
}


// compare haplotype fragment to phased VCF file (partial or complete)
int compare_fragment(FRAGMENT* fragment,VARIANT* varlist,FILE* outfile)
{
	int i=0; char flag = 'x';
	fragment->blocks =1; int matches = 0; int mismatches = 0; int vars=0; int mm=0;
	int m= 0; int sw = 0, bf = 0; 

        for (i=0;i<fragment->variants;i++)
	{
		if (varlist[fragment->alist[i].varid].genotype[1] == '/') continue; 
		if (fragment->alist[i].allele == varlist[fragment->alist[i].varid].genotype[0]) matches++; else mismatches++;
		vars++;

		if (fragment->alist[i].allele != varlist[fragment->alist[i].varid].genotype[0] && m ==0) m = -1; // init
		else if (fragment->alist[i].allele == varlist[fragment->alist[i].varid].genotype[0] && m ==0) m = 1;  // init
		else if (fragment->alist[i].allele == varlist[fragment->alist[i].varid].genotype[0] && m == -1) { m = 1; sw++; } 
		else if (fragment->alist[i].allele != varlist[fragment->alist[i].varid].genotype[0] && m == 1) { m = -1; sw++; } 
	}

	if (matches ==0 && mismatches ==0) return 0; // no phasing information to compare
	else if (matches ==0) fprintf(outfile,"PERFECT-MATCH ");
	else if (mismatches ==0) fprintf(outfile,"PERFECT-MISMATCH ");
	else if (matches <=1) fprintf(outfile,"ERROR_SINGLE_MATCH "); 
	else if (mismatches <=1) fprintf(outfile,"ERROR_SINGLE_MISMATCH "); 
	else if (sw ==1) fprintf(outfile,"ERROR_SINGLE_SWITCH");
	else if (matches > 1 && mismatches > 1) fprintf(outfile,"ERROR_SWITCH "); 
	else fprintf(outfile,"MULT-ERRORS "); 

	if (matches < mismatches) mm = matches; else mm = -1*mismatches;
	fprintf(outfile,"COMPARE_PHASE: %d/%d switches:%d %s..%s | ",mm,vars,sw,varlist[fragment->alist[0].varid].genotype,varlist[fragment->alist[fragment->variants-1].varid].genotype);
        for (i=0;i<fragment->variants;i++) 
	{
		flag = '0';
		if (varlist[fragment->alist[i].varid].genotype[1] == '/') continue; 
		if (fragment->alist[i].allele != varlist[fragment->alist[i].varid].genotype[0] && matches >= mismatches) flag = '-'; 
		if (fragment->alist[i].allele == varlist[fragment->alist[i].varid].genotype[0] && matches < mismatches) flag = '+'; 
		if (flag !='0') fprintf(outfile,"%c:",flag); 
		fprintf(outfile,"%c:%c%c%c:%d ",fragment->alist[i].allele,varlist[fragment->alist[i].varid].genotype[0],varlist[fragment->alist[i].varid].genotype[1],varlist[fragment->alist[i].varid].genotype[2],(int)fragment->alist[i].qv-33);
	}
	fprintf(outfile,"\n");
	return 1;
}

// segfault due to size limit 4096 on length of variant list 08/24/16 
int generate_single_fragment(struct alignedread** readlist,FRAGMENT* flist,VARIANT* varlist,struct FOSMID* fosmid,char* tag)
{
	int s = fosmid->firstread; int e = fosmid->lastread; int length = fosmid->ws; int reads_window = fosmid->reads_window; 
	int start = fosmid->start; int end = fosmid->end;	double read_density = fosmid->mean*1000/BLOCK_SIZE; 

	int j=0,i=0,k=0;
	int unique_variants =1; int hets=0; int counts[4]; int qv=0; int prev = readlist[s]->position;
	FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*4096);
	if (BED_FILE ==1) reads_window = 0; 
	for (k=s;k<e && fragment.variants < 4024;k++)
	{
		i = k; 
		if (readlist[i]->mateindex >= 0 && readlist[i]->IS < 0) continue;  // if mateindex exists and insert size is negative, we will consider it below 

		if (BARCODE ==1 && (readlist[i]->barcode == NULL || strcmp(readlist[i]->barcode,readlist[s]->barcode) != 0)) continue; // ignore this read 
		if (BED_FILE ==1) reads_window++;
		//fprintf(stdout,"%d ",readlist[i]->position-prev); prev = readlist[i]->position;   
		if (readlist[i]->findex >= 0) // read overlaps at least one variant
		{
			for (j=0;j<flist[readlist[i]->findex].variants;j++) 
			{		
				fragment.alist[fragment.variants].varid = flist[readlist[i]->findex].alist[j].varid; 
				fragment.alist[fragment.variants].allele = flist[readlist[i]->findex].alist[j].allele;
				fragment.alist[fragment.variants].qv = flist[readlist[i]->findex].alist[j].qv;
				fragment.variants++;
			}
		}
		// consider the mate of the read but this requires the reads to be linked to each other pairs...
		i = readlist[i]->mateindex; 
		if (i >= 0 && readlist[i]->findex >=0)
		{
			for (j=0;j<flist[readlist[i]->findex].variants;j++) 
			{			
				fragment.alist[fragment.variants].varid = flist[readlist[i]->findex].alist[j].varid; 
				fragment.alist[fragment.variants].allele = flist[readlist[i]->findex].alist[j].allele;
				fragment.alist[fragment.variants].qv = flist[readlist[i]->findex].alist[j].qv;
				fragment.variants++;
			}
		}
		if (fragment.variants >= 4096) { fprintf(stderr,"need to increase variant list size \n"); break; } 
	}
	//fprintf(stdout," mapping positions for interval %s %d-%d len %d \n",readlist[s]->barcode,start,end,end-start);
	if (BED_FILE ==1) read_density = (double)reads_window*1000/length;

	// even if there are two variants covered, it could be a single variant covered twice... so sort list and check
	if (fragment.variants > 2) qsort(fragment.alist,fragment.variants,sizeof(allele),compare_alleles); 
	unique_variants = 0; if (fragment.variants > 0) unique_variants = 1; 
	for (i=1;i<fragment.variants;i++)
	{
		j = fragment.alist[i].varid;
		if (i > 0 && fragment.alist[i].varid != fragment.alist[i-1].varid) unique_variants++;
		if (i > 0 && j == fragment.alist[i-1].varid && fragment.alist[i].allele != fragment.alist[i-1].allele) hets++;
	}
	fosmid->hets = hets; fosmid->unique_variants = unique_variants; 

	// changed to allow printing fragments with single variants 08/25/16 
	if (hets >=(float)unique_variants*0.2  || unique_variants < 1)  // fragment only has single variant or has 2 or more heterzygous variants 
	{
		// such fragments need to be printed in output for debugging and to tag chimeric ones... 
		//fprintf(fragment_file,"#0 -:%d_%d_%d_%0.2f_%d",readlist[s]->position,readlist[e-1]->position,length,read_density,reads_window);
		if ((float)hets/unique_variants >= 0.05 && hets >=2) fprintf(stdout,"HET-VARS ");
		fprintf(stdout,"fragment uvars:%d hets:%d \n",unique_variants,hets);

		//fprintf(fragment_file,"0 -:%d_%d_%d_%0.2f_%d",start,end,length,read_density,reads_window);
		//if (unique_variants < 2) fprintf(fragment_file," vars:%d no-info\n",unique_variants);
		//else fprintf(fragment_file," vars:%d hets:%d\n",unique_variants,hets);

		if (unique_variants >= 2) print_fragment_vars1(&fragment,varlist); 
		free(fragment.alist);
		return 0; 
	}

	FRAGMENT fp; fp.variants = 0; fp.alist = (allele*)malloc(sizeof(allele)*(unique_variants+1)); 
	
	counts[0] = counts[1] = counts[2] = counts[3] = 0;
	counts[(int)fragment.alist[0].allele-48]++; counts[(int)fragment.alist[0].allele-48+2] += (int)fragment.alist[0].qv-QVoffset; 

	j=0;
	for (i=1;i<=fragment.variants;i++)
	{
		if (i==fragment.variants || fragment.alist[i].varid != fragment.alist[i-1].varid) 
		{
			// print consensus base
			if (counts[0] > counts[1] && counts[1] <= 1) 
			{
				fp.alist[j].varid = fragment.alist[i-1].varid; fp.alist[j].allele = '0';
				qv = (QVoffset + counts[2]-counts[3]);
				if (counts[2]-counts[3] >= 60) qv = 60 + QVoffset;
				fp.alist[j].qv = (char)(qv);
				if (qv -QVoffset >= MINQ) j++;
			}
			else if (counts[1] > counts[0] && counts[0] <= 1)
			{
				fp.alist[j].varid = fragment.alist[i-1].varid; fp.alist[j].allele = '1';
				qv = (QVoffset + counts[3]-counts[2]);
				if (counts[3]-counts[2] >= 60) qv = 60 + QVoffset;
				fp.alist[j].qv = (char)(qv);
				if (qv -QVoffset >= MINQ) j++;
			}
			counts[0] = counts[1] = counts[2] = counts[3] = 0;
		}
		if (i < fragment.variants) 
		{
			counts[(int)fragment.alist[i].allele-48]++; 
			counts[(int)fragment.alist[i].allele-48+2] += (int)fragment.alist[i].qv-QVoffset; 
		}
	}
	/*
	*/
	fprintf(stdout,"fragment %d %d \n",unique_variants,j);
	fp.id = (char*)malloc(1024); 
	if (BARCODE ==1) sprintf(fp.id,"%s:%d_%d_%d_%0.2f_%d_%s",varlist[fp.alist[0].varid].chrom,start,end,length,read_density,reads_window,readlist[s]->barcode);
	else sprintf(fp.id,"%s:%d_%d_%d_%0.2f_%d_%0.2f",varlist[fp.alist[0].varid].chrom,start,end,length,read_density,reads_window,fosmid->delta);
	
	int adjusted_hets = print_fragment_vars1(&fragment,varlist); 
	fp.variants = j; 
	if (j >=2) 
	{	
		if (adjusted_hets > 0) fprintf(stdout,"FRAGMENT%s:%d:%d ",tag,adjusted_hets,fp.variants); else fprintf(stdout,"FRAGMENT%s ",tag); 
		print_fragment(&fp,varlist,stdout); 
		//fprintf(stderr,"fragfile %s \n",fragment_file);
		//if (fragment_file != stdout) print_fragment(&fp,varlist,fragment_file); // print fragment to output file, it has no other information, can be used as direct hapcut input
		if (COMPARE_PHASE ==1) compare_fragment(&fp,varlist,stdout);
	}
	else if (j ==1) // fragment with single variant, print it so that it can potentially be merged with adjacent fragments
	{
		fprintf(stdout,"FRAG_SINGLE "); print_fragment(&fp,varlist,stdout); 
		if (COMPARE_PHASE ==1) compare_fragment(&fp,varlist,stdout);
	}
	free(fp.alist); free(fp.id);

	free(fragment.alist); return 1;
}


int print_read(struct alignedread** readlist, int i,int prevpos,FRAGMENT* flist,VARIANT* varlist)
{
	int j=0,v=0;
	fprintf(stdout,"dist %5d %5d %5d ",readlist[i]->position-prevpos,readlist[i]->cluster,readlist[i]->blockid); 
	fprintf(stdout,"%s %s %d-%d %d %d ",readlist[i]->readid,readlist[i]->chrom,readlist[i]->position,readlist[i]->position+readlist[i]->IS,readlist[i]->IS,readlist[i]->mquality);
	for (j=0;j<readlist[i]->cigs;j++) fprintf(stdout,"%d%c",readlist[i]->cigarlist[j]>>4,INT_CIGAROP[readlist[i]->cigarlist[j]&0xf]); 
	if (readlist[i]->findex >= 0) 
	{
		//fprintf(stdout," vars %d ",flist[readlist[i]->findex].variants);
		for (j=0;j<flist[readlist[i]->findex].variants;j++) 
		{			
			v =flist[readlist[i]->findex].alist[j].varid; 
			fprintf(stdout," %d:%c:%d:%s/%s:%c%c%c ",v+1,flist[readlist[i]->findex].alist[j].allele,varlist[v].position,varlist[v].RA,varlist[v].AA,varlist[v].genotype[0],varlist[v].genotype[1],varlist[v].genotype[2]);
		}
	}

	if (readlist[i]->mateindex < 0 || readlist[readlist[i]->mateindex]->findex <0) { fprintf(stdout," \n"); return 1; } 

	i = readlist[i]->mateindex; 
	//fprintf(stdout," || mate:%d:%d:",readlist[i]->position,readlist[i]->mquality);
	//for (j=0;j<readlist[i]->cigs;j++) fprintf(stdout,"%d%c",readlist[i]->cigarlist[j]>>4,INT_CIGAROP[readlist[i]->cigarlist[j]&0xf]); 
	if (readlist[i]->findex >= 0) 
	{
		//fprintf(stdout," vars %d ",flist[readlist[i]->findex].variants);
		for (j=0;j<flist[readlist[i]->findex].variants;j++) 
		{			
			v =flist[readlist[i]->findex].alist[j].varid; 
			fprintf(stdout," %d:%c:%d:%s/%s:%c%c%c ",v+1,flist[readlist[i]->findex].alist[j].allele,varlist[v].position,varlist[v].RA,varlist[v].AA,varlist[v].genotype[0],varlist[v].genotype[1],varlist[v].genotype[2]);
			fprintf(stdout," %d:%c:%d ",v+1,flist[readlist[i]->findex].alist[j].allele,varlist[v].position);
		}
	}
	fprintf(stdout," \n");
	return 1;
}

// print reads within a range (s-e) 
void print_reads_window(struct alignedread** readlist, int s,int e,FRAGMENT* flist,VARIANT* varlist,int if_variant_read)
{
	fprintf(stdout,"cluster chr%s:%d-%d length %d reads %d density %0.4f\n",readlist[s]->chrom,readlist[s]->position,readlist[e-1]->position,readlist[e-1]->position-readlist[s]->position,e-s,(float)(e-s)/(readlist[e-1]->position-readlist[s]->position));
	int i=0,prevpos = -1;
	for (i=s;i<e;i++)
	{
		if (readlist[i]->IS < 0 ||  ((readlist[i]->flag & 1024) ==1024)) continue; 
		if (readlist[i]->IS > 3000) 
		{
			//fprintf(stdout,"P-DEL ");  print_read(readlist,i,prevpos,flist,varlist);
		}
		//else if (readlist[i]->position - prevpos > 1000 && prevpos > 0) print_read(readlist,i,prevpos,flist,varlist);
		//else if (if_variant_read ==0) print_read(readlist,i,prevpos,flist,varlist);
		//else if (readlist[i]->findex >= 0 || (readlist[i]->mateindex > 0 && readlist[readlist[i]->mateindex].findex >=0)) print_read(readlist,i,prevpos,flist,varlist);
		prevpos = readlist[i]->position;
	}
}
