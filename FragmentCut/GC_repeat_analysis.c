
// # of clusters = number of reads for which left distance is >= big_threshold and right distance <= small_threshold 

// calculate distance between adjacent 'GOOD' reads (first in pair, high mapping quality, not PCR duplicate) 
void distance_nextread_dist(struct alignedread** readlist,int s,int e)
{
        int i=s,j=0;
        while (readlist[i]->IS < 0 ||  ((readlist[i]->flag & 1024) ==1024) || readlist[i]->mquality < MIN_MQ) i++; // find first good read

	while (i < e)
	{
                j = i+1; while (j < e && (readlist[j]->IS < 0 || (readlist[j]->flag & 1024) ==1024 || readlist[j]->mquality < MIN_MQ)) j++;  // find next good read in list
                if (j >= e) break; 
		//fprintf(stdout,"dnread i:%d j:%d %0.4f\n",i,j,log10(readlist[j]->position-readlist[i]->position+1));
		i=j; 
        }
}


// functions related to GC content and repeat masking 

// this is done before finding the fragments 
int calculate_block_stats(struct BLOCK_READ* RL, int blocks, REFLIST* reflist,REFLIST* genomemask)
{
	int i=0,j=0,k=0,c1,c2,c; int X0=0,X1=0,GC=0,valid=0; int subsize = BLOCK_SIZE/BS1;

	float counts[GCBINS]; for (i=0;i<GCBINS;i++) counts[i]=0; // readcounts for GC % blocks (0-5,5-10,10-15....95-100)
	int GCbins[GCBINS]; for (i=0;i<GCBINS;i++) GCbins[i] = 0; 
	float totalcounts =0; int totalbins =0;

	// use 0/1 counts for sub-bins to smooth count function | calculate mappability of block and GC content 
	for (i=0;i<blocks;i++) 
	{
		if (genomemask->current >= 0) // genome mask file exists
		{
			X0 =0; X1 = 0; k=0;
			for (j=0;j<BLOCK_SIZE;j++)
			{
				c = genomemask->sequences[genomemask->current][i*BLOCK_SIZE+j]- 63; 
				c1 = c>>3; c2 = c & 7; 
				c1 = c1? 1<<(c1-1) : 0; c2 = c2? 1<<(c2-1) : 0;  
				if (c1 == 1 && c2 < 2) X0++; 
				if (c1 > 1 || c2 > 1 ) X1++; 
				if ((j+1)%subsize==0)  // remove entire sub-block if any base in that sub-block has low mappability... too stringent 
				{
					if (BS1 >=5) 
					{
						if (X1 > 0 && RL[i].subcounts[k] > 0) RL[i].subcounts[k] = 0;
						if (SINGLE_READS ==1 && X1 > 0 && RL[i].subcounts[k+BS1] > 0) RL[i].subcounts[k+BS1] = 0;
						if (X1 > 0) RL[i].adjusted_size -= subsize; 
					}
					k++; X1 = 0; 
				}
			}
			RL[i].mappability = (float)X0/BLOCK_SIZE; 
			if (RL[i].adjusted_size >= BLOCK_SIZE) RL[i].adjusted_size = (int)(RL[i].mappability*(float)BLOCK_SIZE); 
		}
		// smoothed read counts	
		RL[i].reads = 0; 
		for (k=0;k<BS1;k++) 
		{
			if (RL[i].subcounts[k] >= 1) RL[i].reads +=1;
			if (SINGLE_READS==1 && RL[i].subcounts[k+BS1] >= 1) RL[i].reads +=1;
		}

		if (reflist->ns > 0) 
		{
			GC = 0; RL[i].GC = 0;
			for (j=0;j<BLOCK_SIZE;j++)
			{
				if (reflist->sequences[reflist->current][i*BLOCK_SIZE+j] == 'G' || reflist->sequences[reflist->current][i*BLOCK_SIZE+j] =='C') GC +=1;
			}
			RL[i].GC = (float)GC/BLOCK_SIZE; GC = (int)(RL[i].GC*GCBINS);
			if (RL[i].mappability > 0.99) { counts[GC] += RL[i].reads; GCbins[GC] +=1; } 
			//fprintf(stderr,"%d %d mappability %f GC %d %d %d \n",i*BLOCK_SIZE,i*BLOCK_SIZE+j,RL[i].mappability,GC,RL[i].reads,RL[i].adjusted_size);	
		}

		if (RL[i].mappability < MIN_MAPPABILITY) RL[i].valid = '0';
		if (RL[i].adjusted_size < subsize && BS1 >= 5) RL[i].valid = '0'; 
		if (RL[i].GC <= MIN_GC || RL[i].GC >= MAX_GC) RL[i].valid = '0'; 
		else if (reflist->ns > 0) 
		{
			totalcounts += RL[i].reads; totalbins++; 
		}

		if (RL[i].valid == '1') valid++;

		//fprintf(stdout," block %d reads %d %d first %d last %d \n",i+j,RL[i+j].reads,empty,RL[i+j].firstread,RL[i+j].lastread);
	}
	double maxGC = 0; double meanRD = totalcounts/(totalbins+1);
	for (i=1;i<19;i++) 
	{
		fprintf(stdout,"GCanalysis before DP %2d-%2d reads %6.0f bins %6d average %0.3f GM %0.3f\n",i*5,i*5+5,counts[i],GCbins[i],(float)counts[i]/(GCbins[i]+1),meanRD); 
		counts[i] /= (GCbins[i]+1); if ( counts[i] > maxGC) maxGC = counts[i]; 
	}

	for (i=0;i<blocks;i++) 
	{
		GC = (int)(RL[i].GC*GCBINS);
		if (GC_CORRECTION ==1) RL[i].GC_corr = meanRD/counts[GC]; // for low GC content bins, this value is > 1 and vice versa 
		else RL[i].GC_corr = 1.0;
		RL[i].GC_corr_log = log(RL[i].GC_corr); 
		// correction for GC coverage...
		//RL[i].GC_corr = 1; RL[i].reads *= meanRD/counts[(int)(RL[i].GC*BLOCK_SIZE/5)]; // correction for GC content...
		//fprintf(stderr,"%d %d mappability %f GC %0.3f %d %d %0.3f\n",i*BLOCK_SIZE,i*BLOCK_SIZE+j,RL[i].mappability,RL[i].GC,RL[i].reads,RL[i].adjusted_size,RL[i].GC_corr);	
	}

	return valid;
}



// calculate background read density using blocks of size 1KB, fraction of blocks that are empty (0 reads) = e^-l according to poisson distribution 
// this is used for likelihood calculations
double calculate_background_RD(struct BLOCK_READ* RL,int blocks,int validblocks)
{
	int i=0,j=0; 
	double empty_blocks=0,nonempty_blocks=0, W = 0; int non_empty = 0;
	// use bigger block size of 1000bp to calculate background density, if block of 1KB is empty -> due to background....
	// formula = lambda = -1*log(empty_blocks/total_blocks - heavy blocks)
	// some blocks could be empty due to mappability, exclude them or at least weight them appropriately
	for (i=0;i<blocks && i+BF < blocks;i+=BF) 
	{
		non_empty = 0; for (j=0;j<BF;j++) non_empty += RL[i+j].reads*RL[i+j].GC_corr; 
		W = 0; for (j=0;j<BF;j++) 
		{ 
			if (RL[i+j].valid == '1') W += RL[i+j].mappability;
		}
		if (non_empty ==0) empty_blocks += W/BF; 
		else if (non_empty < 3*BF) nonempty_blocks += W/BF; // less than 3 reads per 100 bp BLOCK
	}
	double read_density = log((nonempty_blocks+empty_blocks)) - log(empty_blocks);  read_density /= BF;
	//read_density = 0.007; // average number of reads per 100 bp window 
	fprintf(stderr,"1KB-blocks %0.1f %0.1f %0.8f valid blocks %d:%f\n",nonempty_blocks+empty_blocks,empty_blocks,read_density,validblocks,(float)validblocks/blocks);
	//////// block based analysis gives us probability that a block has 1 or more reads -> background poisson distribution for reads... ////////

	return read_density;

}

// GC read depth statistics across all bins, ignores low mappability bins (< 0.99) for calculation 
int calculate_GC_stats(struct BLOCK_READ* RL, int blocks,struct FOSMID* fosmidlist,int fosmids)
{
	int i=0,j=0;
	float counts[GCBINS]; for (i=0;i<GCBINS;i++) counts[i]=0.01; // readcounts for GC % blocks (0-5,5-10,10-15....95-100)
	float GCbins[GCBINS]; for (i=0;i<GCBINS;i++) GCbins[i] = 0; 
	float totalcounts =0; float totalbins =0; float mean; int GC; 
	float reads=0,WS = 0;

	for (i=0;i<fosmids;i++)
	{
		mean = fosmidlist[i].mean;
		for (j=fosmidlist[i].firstblock;j<=fosmidlist[i].lastblock;j++)
		{
			if (RL[j].mappability < 0.99) continue; 
			reads += RL[j].reads; WS += BLOCK_SIZE; 
		}
		mean = reads/WS; 
		for (j=fosmidlist[i].firstblock;j<=fosmidlist[i].lastblock;j++)
		{
			if (RL[j].mappability < 0.99) continue; 
			GC = (int)(RL[j].GC*GCBINS);
			counts[GC] += RL[j].reads; GCbins[GC] += mean;  
			totalcounts += RL[j].reads; totalbins += mean; 
		}
	}
	double maxGC = 0; double meanRD = totalcounts/totalbins;
	for (i=1;i<19;i++) 
	{
		fprintf(stdout,"GC %2d-%2d reads %6.0f bins %6f average %0.3f GM %0.3f\n",i*5,i*5+5,counts[i],GCbins[i],(float)counts[i]/(GCbins[i]+1),meanRD); 
		counts[i] /= (GCbins[i]+1); 
	}


}
