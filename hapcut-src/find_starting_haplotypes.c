
///////////////////////////////////////////////////// CODE TO GENERATE INITIAL HAPLOTYPE SOLUTION USING A CLUSTERING METHOD//////////////////////////////


int sample_block(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int start, int end,char* h1,int snps)
{
	int j=0,t=0,t1=0,t2=0,choice=0; char c;
	double parr0=0,parr1=0, sum=0,rd=0,b;
	if (snpfrag[snpfrag[start].component].best_mec <= 0) return -1; if (snpfrag[start].frags < 1) return -1;
	for (j=snpfrag[start].ff;j<fragments;j++)
	{
		if (Flist[j].list[0].offset > end) break;
		if (Flist[j].component != snpfrag[start].component) continue;
		for (t1=0;t1<Flist[j].blocks;t1++)
		{
			for (t2=0;t2<Flist[j].list[t1].len;t2++)
			{
				t = Flist[j].list[t1].offset + t2;
				if (t >= start && t <= end)
				{
					b = log10(1-Flist[j].list[t1].pv[t2]) - log10(Flist[j].list[t1].pv[t2]);
					//if (PMODEL ==0)   b = log10(1-snpfrag[t].pv) - log10(snpfrag[t].pv);
					c = Flist[j].clust;  if (Flist[j].list[t1].hap[t2] != h1[t]) c=  (char)(49 - (int)(c-48));
					if (c =='0') parr0 += b; else parr0 -= b;
				}
			}
		}
	}
	if (parr0 >0 ) { parr1 -= parr0; parr0 = 0; } rd = drand48(); sum = pow(10,parr0) + pow(10,parr1);
	if (rd*sum >= pow(10,parr0))
	{
		for (j=start;j<=end;j++) { if (h1[j] == '-') continue; if (h1[j] =='0') h1[j] = '1'; else h1[j] = '0';  }  choice = 1;
	}
	return choice;
}

void frag_cluster_initialize(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,char* h1,int snps,struct BLOCK* clist,int comps)
{
	int i=0,j=0,k=0,a,b,t1=0,t2=0,match=0, mismatch=0,iter=0,assigned=0,last;
	char* thap = (char*)malloc(snps); 
	for (i=0;i<snps;i++) thap[i] = '-';
	for (i=0;i<snps;i++) h1[i] = '-';
	// cluster fragments within each connected component using pairwise overlap and initialize the starting haplotype from this solution 

	for (k=0;k<comps;k++)
	{
		for (i=0;i<clist[k].frags;i++) Flist[clist[k].flist[i]].clust = 'x'; // initialized 
		Flist[clist[k].flist[0]].clust = '0'; assigned=10; iter =0;
		while (assigned >0 || iter < 5)
		{
			assigned =0; iter++;
			for (i=0;i<clist[k].frags-1;i++)
			{
				t1 = clist[k].flist[i]; 	//if (Flist[t1].l == 'x') continue;
				for (j=i+1;j<clist[k].frags;j++)
				{
					last = Flist[t1].list[Flist[t1].blocks-1].offset+Flist[t1].list[Flist[t1].blocks-1].len; 
					t2 = clist[k].flist[j]; match =0; mismatch =0; 
					if (Flist[t2].list[0].offset >= last) break; 
					if (Flist[t1].clust != 'x' && Flist[t2].clust != 'x') continue;
					if (Flist[t1].clust == 'x' && Flist[t2].clust == 'x') continue;

					for (a=0;a<Flist[t1].blocks;a++)
					{
						for (b=0;b<Flist[t1].list[a].len;b++) thap[Flist[t1].list[a].offset+b] = Flist[t1].list[a].hap[b]; 
					}
					for (a=0;a<Flist[t2].blocks;a++)
					{
						for (b=0;b<Flist[t2].list[a].len;b++) { if (thap[Flist[t2].list[a].offset+b] == Flist[t2].list[a].hap[b]) thap[Flist[t2].list[a].offset+b] = 'M'; else thap[Flist[t2].list[a].offset+b] = 'm';}
					}
					for (a=0;a<Flist[t1].blocks;a++)
					{
						for (b=0;b<Flist[t1].list[a].len;b++) { if (thap[Flist[t1].list[a].offset+b] == 'M') match++; else if (thap[Flist[t1].list[a].offset+b] == 'm') mismatch++; }
					}
					if (match + mismatch ==0) continue;
					if (match != mismatch  && (Flist[t1].clust =='x' || Flist[t2].clust =='x')) assigned++;
					if (match - mismatch >= 1 && Flist[t1].clust != 'x') Flist[t2].clust = Flist[t1].clust;
					else if (mismatch - match >=1 && Flist[t1].clust =='0') Flist[t2].clust = '1';
					else if (mismatch -match >=1 && Flist[t1].clust =='1') Flist[t2].clust = '0';
					else if (match - mismatch >= 1 && Flist[t2].clust != 'x') Flist[t1].clust = Flist[t2].clust;
					else if (mismatch - match >=1 && Flist[t2].clust =='0') Flist[t1].clust = '1';
					else if (mismatch -match >=1 && Flist[t2].clust =='1') Flist[t1].clust = '0';
				}
			} //fprintf(stdout,"cluster %d frags %d assigned %d \n",k,clist[k].frags,assigned);

		}
	}
	for (i=0;i<fragments;i++) 
	{
		//														fprintf(stdout,"%d %d %c ",i,Flist[i].blocks,Flist[i].clust);for (t=0;t<Flist[i].blocks;t++) fprintf(stdout,"| %d %s ",Flist[i].list[t].offset,Flist[i].list[t].hap); fprintf(stdout,"\n");
		for (t1=0;t1<Flist[i].blocks;t1++)	
		{
			for (t2=0;t2<Flist[i].list[t1].len;t2++)
			{
				if (Flist[i].clust == '0' && h1[Flist[i].list[t1].offset + t2] == '-') h1[Flist[i].list[t1].offset + t2] = Flist[i].list[t1].hap[t2]; 
				if (Flist[i].clust == '1' && h1[Flist[i].list[t1].offset + t2] == '-' && h1[Flist[i].list[t1].offset + t2] =='1') h1[Flist[i].list[t1].offset + t2] ='0';
				if (Flist[i].clust == '1' && h1[Flist[i].list[t1].offset + t2] == '-' && h1[Flist[i].list[t1].offset + t2] =='0') h1[Flist[i].list[t1].offset + t2] ='1';
				// Flist[i].list[t1].pv[t2] = 0.05; bug here do not change pv HERE !!!!
			}
		}
	}
	for (i=0;i<fragments;i++) { if (Flist[i].clust == 'x') Flist[i].clust = (char)(48+(int)(drand48()*2-0.00001)); }
	//for (i=0;i<snps;i++) { if (snpfrag[i].frags < 1 || h1[i] != '-') continue; h1[i] = '0';   sample_block(Flist,fragments,snpfrag,i,i,h1,snps); }
	for (i=0;i<snps;i++) { if (snpfrag[i].frags < 1 ) continue; h1[i] = '0';   sample_block(Flist,fragments,snpfrag,i,i,h1,snps); }


}

int edge_compare(const void *a,const void *b)
{
	const struct edge *ia = (const struct edge*)a;
	const struct edge *ib = (const struct edge*)b;
	return ia->snp - ib->snp;
}



////// reading fragment and variant file functions //////////////////////////////

int fragment_compare(const void *a,const void *b)
{
	const struct fragment *ia = (const struct fragment*)a;
	const struct fragment *ib = (const struct fragment*)b;
	if (ia->list[0].offset == ib->list[0].offset)
	{
		return ia->list[ia->blocks-1].offset + ia->list[ia->blocks-1].len - ib->list[ib->blocks-1].offset - ib->list[ib->blocks-1].len;
		//return ia->blocks - ib->blocks;
	}
	else return ia->list[0].offset - ib->list[0].offset;
}
