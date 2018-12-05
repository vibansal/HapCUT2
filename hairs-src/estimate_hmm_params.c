
//char INT_CIGAROP[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};

int BTI[] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 1, 0, 2,  0, 0, 0, 3,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  4, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  // A | C | G | T
        0, 1, 0, 2,  0, 0, 0, 3,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  4, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  // a | c | g | t 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

typedef struct // probabilities are in logspace 
{
        int states; // match, insertion,deletion 
        double** TRS; // transition probabilities between states 
        double** MEM; // probability of A->C, A->G (match state)
        double match; double mismatch; double insertion; double deletion; // emission probs
} Align_Params;

extern Align_Params *AP;

// allocate memory and initializes
Align_Params* init_params()
{
        int i=0,j=0;
        Align_Params *AP = malloc(sizeof(Align_Params));
        (*AP).states = 3;
        (*AP).TRS = calloc((*AP).states,sizeof(double*)); // match =0, ins =1, del = 2
        for (i=0;i<AP->states;i++) AP->TRS[i] = calloc(sizeof(double),AP->states);
        //  initialize alignment parameters data structure 
        AP->match = log(0.979); AP->mismatch = log(0.007);
        AP->deletion = log(1); AP->insertion =log(1);
        AP->TRS[0][0] = log(0.879); AP->TRS[0][1] = log(0.076); AP->TRS[0][2] = log(0.045);
        AP->TRS[1][0] = log(0.865); AP->TRS[1][1] = log(0.135);
        AP->TRS[2][0] = log(0.730); AP->TRS[2][2] = log(0.27);
        AP->MEM = calloc(4,sizeof(double*)); for (i=0;i<4;i++) AP->MEM[i] = calloc(4,sizeof(double));
        for (i=0;i<4;i++)
        {
                for (j=0;j<4;j++)
                {
                        if (i==j) AP->MEM[i][j] = AP->match; else AP->MEM[i][j] = AP->mismatch;
                }
        }
        return AP;
}

// for every HP of length 'k', tabulate # of times read-correctly and indel error 

int estimate_counts_read(struct alignedread* read,REFLIST* reflist,int* emission_counts,int* trans_counts,int* indel_lengths)
{
	//static int emission_counts[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // AA AC AG AT | CA,CC,CG,CT |... 16 pairs im match state 
	//static int trans_counts[9] = {0,0,0,0,0,0,0,0,0}; // M->M, M->I,M-D | I->M, I->D, I->I 
	// 0 = M | I =1, D = 2 from cigar encoding 
	// read->strand is also needed for emission counts 
	int b1=0,b2=0;

	int current = reflist->current; 
	if (current < 0 || reflist->used[current] ==0) return -1;
	int f=0;
        int i=0,t=0, l1=0,l2=0; 
	int l=0; int m=0; int op;
	int prevop = -1;
        for (i=0;i<read->cigs;i++)
        {
                op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4;
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
                {
			if (prevop == 1) { trans_counts[1*3+0] +=1; indel_lengths[20] +=1; }
			else if (prevop == 2) { trans_counts[2*3+0] +=1; indel_lengths[0] +=1; }
			trans_counts[0*3+0] += l-1; // l-1 self loops 

                        for (t=0;t<l;t++)
                        {
				//fprintf(stdout,"current %d len %d %d %d %c\n",current,reflist->lengths[current],read->position,l2,reflist->sequences[current][0]);
				b1 = BTI[reflist->sequences[current][read->position+l2+t-1]]; 
				b2 = BTI[read->sequence[l1+t]]; 
				if (b1 > 0 && b2 > 0 && b1 <= 4 && b2 <= 4) 
				{
					if ((read->flag & 16) ==0) emission_counts[(b1-1)*4 + (b2-1)] +=1; 
					else emission_counts[(4-b1)*4 + (4-b2)] +=1; 
				}
                        }
                }
		if (op == BAM_CDEL) 
		{
			/*
			for (t=-5;t<0;t++) fprintf(stderr,"%c",reflist->sequences[current][read->position+l2+t-1]);
			fprintf(stderr,"|"); for (t=0;t<l;t++) fprintf(stderr,"%c",reflist->sequences[current][read->position+l2+t-1]);
			fprintf(stderr,"|"); for (t=l;t<l+10;t++) fprintf(stderr,"%c",reflist->sequences[current][read->position+l2+t-1]);
			fprintf(stderr," DEL %d %d\n",l,l2+read->position-1);
			*/
			
			if (prevop == BAM_CMATCH || prevop == BAM_CEQUAL || prevop == BAM_CDIFF) trans_counts[0*3 + 2] +=1; 
			if (l < 10) trans_counts[2*3+2] +=l-1;
			if (l < 20) indel_lengths[l] +=1;
		}
		if (op == BAM_CINS)
		{
			/*
			for (t=-5;t<0;t++) fprintf(stderr,"%c",reflist->sequences[current][read->position+l2+t-1]);
			fprintf(stderr,"|"); for (t=0;t<l;t++) fprintf(stderr,"%c",read->sequence[l1+t]);
			fprintf(stderr,"|"); for (t=l;t<l+10;t++) fprintf(stderr,"%c",reflist->sequences[current][read->position+l2+t-1]);
			fprintf(stderr," INS %d %d\n",l,read->position+l2-1);
			*/
			if (prevop == BAM_CMATCH || prevop == BAM_CEQUAL || prevop == BAM_CDIFF) trans_counts[0*3 + 1] +=1; 
			if (l < 10) trans_counts[1*3+1] +=l-1;
			if (l < 20) indel_lengths[l+20] +=1;
		}
                if (op == BAM_CMATCH || op >= 7) { l1 +=l; l2 +=l; }
                else if (op == BAM_CDEL || op == BAM_CREF_SKIP) l2 +=l;
                else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)  l1 += l;
		prevop = op;
        }
	return 1;
}

void print_error_params(int* emission_counts,int* trans_counts,int* indel_lengths,Align_Params* AP)
{
	char ITB[] = {'A','C','G','T'};
	char state[] = {'M','I','D'};
	int pr =0,b1=0,b2=0;
	for (b1=0;b1<4;b1++)
	{
		int total =0.01;
		for (b2=0;b2<4;b2++) total += emission_counts[b1*4+b2]; 
		for (b2=0;b2<4;b2++) 
		{
			fprintf(stderr,"%c -> %c %0.4f | ",ITB[b1],ITB[b2],(float)emission_counts[b1*4+b2]/total);
			AP->MEM[b1][b2] = log((float)emission_counts[b1*4+b2]/total);
		}
		fprintf(stderr,"\n");
	}
	for (b1=0;b1<3;b1++)
	{
		int total =0.01;
		for (b2=0;b2<3;b2++) total += trans_counts[b1*3+b2]; 
		for (b2=0;b2<3;b2++) 
		{
			fprintf(stderr,"%c -> %c %0.4f | ",state[b1],state[b2],(float)trans_counts[b1*3+b2]/total);
			AP->TRS[b1][b2] = log((float)trans_counts[b1*3+b2]/total);
		}
		fprintf(stderr,"\n");
	}
	//for (pr=0;pr<16;pr++) fprintf(stdout,"%c:%c %d \n",ITB[pr/4],ITB[pr%4],emission_counts[pr]);
	fprintf(stderr,"emission for match state \nreads transition counts \n");
	//for (pr=0;pr<9;pr++) fprintf(stdout,"%c:%c %d \n",state[pr/3],state[pr%3],trans_counts[pr]);
	//for (pr=0;pr<10;pr++) fprintf(stderr,"%d del: %d ins: %d \n",pr,indel_lengths[pr],indel_lengths[pr+20]);

}


// reflist only has sequences for contigs in VCF, so this can cause segfault
int realignment_params(char* bamfile,REFLIST* reflist,char* regions,Align_Params* AP)
{
	struct alignedread* read = (struct alignedread*) malloc(sizeof (struct alignedread));
	int i=0;
        int emission_counts[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // AA AC AG AT | CA,CC,CG,CT |... 16 pairs im match state 
        int trans_counts[] = {0,0,0,0,0,0,0,0,0}; // M->M, M->I,M-D | I->M, I->D, I->I 
        int *indel_lengths = calloc(40,sizeof(int));
        for (i=0;i<40;i++) indel_lengths[i]=0;
        // 0 = M | I =1, D = 2 from cigar encoding 
	char* newregion = calloc(4096,sizeof(char));

	samfile_t *fp;
	bam1_t *b;
	bam_iter_t iter;
        bam_index_t *idx; // bam file index
	int ret; int ref=-1,beg=0,end=0;

	if ((fp = samopen(bamfile, "rb", 0)) == 0) 
	{
        	fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
	        return -1;
    	}

        if (regions != NULL) 
	{
		strcpy(newregion,regions); 
		for (i=0;i<strlen(newregion);i++) 
		{
			if (newregion[i] == ':') break;
		}
		newregion[i+1] = '\0';
		//fprintf(stderr,"newrion %s %s \n",regions,newregion);
		if ((idx = bam_index_load(bamfile)) ==0) { fprintf(stderr,"unable to load bam index for file %s\n",bamfile); return -1; }
		bam_parse_region(fp->header,newregion,&ref,&beg,&end);
	        if (ref < 0) { fprintf(stderr,"invalid region for bam file %s \n",regions); return -1; }
		iter = bam_iter_query(idx,ref,beg,end);
	}

	b = bam_init1();  // samread(fp,b)

	int prevtid =-1;
	int reads=0,ureads=0,useful=0;

	reflist->current = -1;

	while (reads < 10000)	
	{
		if (ref < 0) ret = samread(fp,b);  // read full bam file 
	        else ret = bam_iter_read(fp->x.bam,iter,b);  // specific region 
		if (ret < 0) break;
	        fetch_func(b, fp, read);
                if ((read->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | 2048)) || read->mquality < MIN_MQ) 
		{
	            free_readmemory(read);  	    continue;
        	}

		if (read->tid != prevtid) 
		{
			reflist->current = -1;
			for (i = 0; i < reflist->ns; i++) 
			{
				if (strcmp(reflist->names[i], read->chrom) == 0) 
				{
				    reflist->current = i;
				    break;
				}
			}
            	}
		//if (reflist->current < 0) continue;
		reads +=1;
		//print_read_debug(read);		
	        useful = estimate_counts_read(read,reflist,emission_counts,trans_counts,indel_lengths);
		if (useful > 0) ureads +=1;
		prevtid = read->tid;
	}
	fprintf(stderr,"using %d reads to estimate realignment parameters for HMM \n",reads);
	if (ureads > 1000) print_error_params(emission_counts,trans_counts,indel_lengths,AP);
	
	bam_destroy1(b); samclose(fp); free(newregion);

}
