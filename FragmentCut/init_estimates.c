int cmpfunc (const void * a, const void * b) {   return ( *(int*)a - *(int*)b );}

// function to estimate # of fragments, intra-fragment gap and background read density 
int estimate_lambda_b(struct BLOCK_READ* RL,int blocks,double estimates[])
{
	int i=0,j=0,k=0,prev=-1,peaks=0;
	for (j=0;j<blocks;j++) 
	{
		if (RL[j].valid == '0') continue; 
		if (RL[j].reads > 0) k++; 
	}
	estimates[0] = 3000; estimates[1] = 75000;
	int* peakdist = calloc(k+1,sizeof(int)); k =0; 

	for (j=0;j<blocks;j++) 
	{
		if (RL[j].valid == '0') continue; 
		if (RL[j].reads >0) 
		{ 
			if (prev >= 0) peakdist[k++] = j-prev; 
			prev = j; 
		} 
	}	
	peaks =k;
	qsort(peakdist,k,sizeof(int),cmpfunc); 

	// exponential distribution has tail prob = e^(-lambda x), use values of 75k and 150k to calculate ratio 
	double tail_counts[5] = {0,0,0,0,0}; double tail_points[5] = {10000,20000,30000,40000,50000}; // 50K,75K, 100K, 150K
	for (j=peaks-1;j>=0;j--) 
	{
		if (peakdist[j]*BLOCK_SIZE > tail_points[0]) tail_counts[0]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[1]) tail_counts[1]+=1;
		if (peakdist[j]*BLOCK_SIZE > tail_points[2]) tail_counts[2]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[3]) tail_counts[3]+=1;
		if (peakdist[j]*BLOCK_SIZE >tail_points[4]) tail_counts[4]+=1;
	}
	double exp_mean1  = (tail_points[2]-tail_points[0])/log(tail_counts[0]/tail_counts[2]); 
	double exp_mean2  = (tail_points[3]-tail_points[1])/log(tail_counts[1]/tail_counts[3]); 
	double exp_mean3  = (tail_points[4]-tail_points[2])/log(tail_counts[2]/tail_counts[4]); 
	
	double intragap_mean = peakdist[(peaks-(int)tail_counts[1])/2]*BLOCK_SIZE/log(2); // ln(2)/lambda 
	double intragap_mean1 = peakdist[(peaks)/2]*BLOCK_SIZE/log(2); // ln(2)/lambda 

	for (i=0;i<4;i++) fprintf(stderr,"%0.0f:%0.0f ",tail_points[i],tail_counts[i]);
	fprintf(stderr,"tail counts | mean of intra-gap %0.2f %0.2f %d mean_exp %f %f %f\n",intragap_mean,intragap_mean1,peaks,exp_mean1,exp_mean2,exp_mean3);
	for (i=0;i<4;i++) fprintf(stdout,"%0.0f:%0.0f ",tail_points[i],tail_counts[i]);
	fprintf(stdout,"tail counts mean1 %0.3f %d mean_exp %f %f %f\n",intragap_mean,peaks,exp_mean1,exp_mean2,exp_mean3);
	fprintf(stderr,"mean_exp %f %f %f\n",(float)BLOCK_SIZE/exp_mean1,(float)BLOCK_SIZE/exp_mean2,(float)BLOCK_SIZE/exp_mean3);
	free(peakdist);

	estimates[0] = intragap_mean; estimates[1]= exp_mean1;  
	//if (tail_counts[3] >= 50) estimates[1] = exp_mean2; else estimates[1] = exp_mean1;
}


// print intra-fragment gaps and coverage per block in compact format, useful for quick visualization 
void print_block(struct BLOCK_READ* RL,int start,int end,int flag)
{
	int j=0;int empty =0,empty1=0; float gaps=0; 
	int prev = 0, peaks = 0, sign = 0; // number of changes in coverage between adjacent bins, only positive, if two changes very close, disregard...
	if (flag ==0) 
	{
		for (j=start;j<=end;j++) 
		{ 
			if (RL[j].valid == '1' && RL[j].counts ==0) fprintf(stdout,"|."); 
			else if (RL[j].valid == '1') fprintf(stdout,"|%d",RL[j].counts); 
			else fprintf(stdout,"*"); 
		} fprintf(stdout," counts\n");
	}
	else
	{
		double seglen[3]={0,1.0e-6,0}, gaplen[3] = {0,1.0e-6,0}; int s=0;
		fprintf(stdout,"\ncounts-gap    ");
		for (j=start;j<=end;j++) 
		{ 
			if (RL[j].valid == '0') empty1++; 
			else 
			{
				// mean gap and segment lengths	
				if ((RL[j].counts==0 || j == end) && s > 0) { seglen[0] +=s; seglen[2] += s*s; seglen[1] += 1; s=0; } 
				if (RL[j].counts > 0) 
				{
					if (empty > 0) { gaplen[0] += empty; gaplen[2] += empty*empty; gaplen[1] +=1; } 
					s+=1;
				}

				if (RL[j].counts ==0) 
				{ 
					if (empty ==0) gaps += RL[j].mappability; empty++; 
				} 
				else
				{
					if (empty1+ empty < 2) fprintf(stdout,"|%d|",RL[j].counts);
					else if (empty1 >= 10) fprintf(stdout,"----%d+%d----|%d|",empty*BLOCK_SIZE,empty1*BLOCK_SIZE,RL[j].counts);  
					else fprintf(stdout,"----%d----|%d|",empty*BLOCK_SIZE,RL[j].counts);  
					empty =0; empty1=0; 
				}
				//if (RL[j].peak > 0) { fprintf(stdout,"+"); peaks++; } 
				if (RL[j].counts > prev && sign ==0) { sign = 1; } 
				else if (RL[j].counts < prev && sign ==1) { sign = 0; fprintf(stdout,"+"); peaks++; } 
				//if (RL[j].counts > prev && sign ==0) { fprintf(stdout,"+"); peaks++; sign = 1; } 
				//else if (RL[j].counts < prev) sign = 0; 
				prev = RL[j].counts; 
			}
		} 
		double std1 = sqrt(seglen[2]/seglen[1] - seglen[0]*seglen[0]/(seglen[1]*seglen[1]));
		double std2 = sqrt(gaplen[2]/gaplen[1] - gaplen[0]*gaplen[0]/(gaplen[1]*gaplen[1]));
		if (empty > 0) fprintf(stdout,"----%d----",empty*BLOCK_SIZE);  
		if (flag >= 2) fprintf(stdout," changes %d gaps %0.2f %d ld %0.2f seg:%0.3f:%d:%0.3f gap:%0.3f:%d:%0.2f \n\n",peaks,gaps,(end-start-empty1)*BLOCK_SIZE,(float)((end-start-empty1)*BLOCK_SIZE)/(peaks+1.0e-4),seglen[0]/seglen[1],(int)seglen[1],std1,BLOCK_SIZE*gaplen[0]/gaplen[1],(int)gaplen[1],std2);
		else fprintf(stdout,"inter-frag %0.2f %d %0.2f\n\n",gaps,(end-start-empty1)*BLOCK_SIZE,(float)((end-start-empty1)*BLOCK_SIZE)/(gaps+1.0e-4));
		
	}
}
