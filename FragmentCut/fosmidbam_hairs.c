#include<stdint.h>
#include "fragments.h"
#include "print_clusters.c"
#include "GC_repeat_analysis.c"
#include "analyze_clusters.c"
//#include "sissor-code.c"
#include "init_estimates.c"

// use "tail" of inter-fragment start distances to estimate lambda for # of fragments per pool.... rather than # of fragments estimated... 
// calculate mean/median/variance of fragment lengths, read-density-per-fragment mean, lambda_b, fragment_penalty 
int find_best_segmentation(struct BLOCK_READ* RL,int blocks,int* BL,int validblocks,struct FOSMID* fosmidlist,int refflag,double* bgdensity,double length_pdf[],double block_penalty,double* mean_fl)
{
	fprintf(stderr,"DP iteration:\t");
	int i=0,j=0,k=0,c=0,t=0; int fosmids=0; int flag =0; int block0_start = -1;  int first=0,last=0,l=0;
	double ws=0,bg_reads=0,bg_bins=0; // read counts outside boundaries of fragments
        int genome_covered = 0; // genome spanned by fragments in basepairs 
	double mean_fraglen=0;
	double bg_ll=0,seg_ll=0;
	double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);
	double read_density = *bgdensity; // passed to function, updated after function completes
	double mean[3];
	double p1,delta;
	int delta_counts[50]; for (i=0;i<50;i++) delta_counts[i] = 0; 
	
	for (i=0;i<300;i++) length_pdf[i] = 0; 

	i=validblocks-1; while (i >= 0) // start from last block in the list 
	{
		j = RL[BL[i]].previous[0];    // both i and j are indexes in smaller list BL[i] and BL[j] are indexes into full list 
		if (j < i-1) // (j-i) is a haplotype fragment
		{
			//fprintf(stderr,"long block %d %d \n",BL[i],BL[j]); // this corresponds to the blocks not included in segmentation...
			// need to add code to add windows ignored due to GC content but adjacent to valid islands 
			first = RL[BL[j+1]].firstread; k=1; while(first < 0) { first = RL[BL[j+1]+k].firstread; k++; } 
			last = RL[BL[i]].lastread; k= 1;  while(last < 0) { last = RL[BL[i]-k].lastread; k++; } 

			// last block 'i' is part of segment [] rather than [)	
			fosmidlist[fosmids].firstblock = BL[j+1]; fosmidlist[fosmids].lastblock = BL[i]; 
			fosmidlist[fosmids].fv = j+1; fosmidlist[fosmids].lv = i; 
			fosmidlist[fosmids].start = BL[j+1]*BLOCK_SIZE; fosmidlist[fosmids].end = (BL[i]+1)*BLOCK_SIZE; 
			fosmidlist[fosmids].firstread = first; fosmidlist[fosmids].lastread = last;
			fosmidlist[fosmids].mean  = RL[BL[i]].mean; fosmidlist[fosmids].reads_window = RL[BL[i]].reads_window; 
			fosmidlist[fosmids].ws =  (BL[i]-BL[j+1]+1)*BLOCK_SIZE; fosmidlist[fosmids].score_window = RL[BL[i]].score_window; 


			mean_fraglen += fosmidlist[fosmids].ws; 
		
			seg_ll =  fosmidlist[fosmids].score_window;
			bg_ll =0; for (t=fosmidlist[fosmids].fv;t<=fosmidlist[fosmids].lv;t++) bg_ll += RL[BL[t]].bgscore; 
			delta = seg_ll+block_penalty-bg_ll; p1 = exp(-1.0*delta)/(exp(-1.0*delta)+1.0); 
			fosmidlist[fosmids].delta = delta;
			if (delta >= 49) delta_counts[49]++; else if (delta < 0) fprintf(stderr,"BUG\n"); else delta_counts[(int)delta]++;
			if (delta <6) 
			{
				//fprintf(stderr,"BUG delta %f %f %f i:%d\n",delta,seg_ll,bg_ll,fosmids);
				//for (t=fosmidlist[fosmids].fv;t<=fosmidlist[fosmids].lv;t++) { bg_reads += (double)RL[BL[t]].reads*p1; bg_bins += p1; } //RL[BL[t]].adjusted_size*p1/BLOCK_SIZE;  } 
			}
			if (p1 < 0.01) // only add solid fragments  
			{
				if (fosmidlist[fosmids].ws < 200000) { length_pdf[(int)(fosmidlist[fosmids].ws/1000)+1] +=1; length_pdf[0] +=1; }
				fosmids++; 
				genome_covered += BL[i]-BL[j+1]+1; // under fragments
			}	
			flag =0; block0_start = -1;
		}
		else // empty block background
		{
                        if(RL[BL[i]].adjusted_size >= BLOCK_SIZE-1) { bg_reads += RL[BL[i]].reads; bg_bins += 1; } 
			//RL[BL[i]].adjusted_size/BLOCK_SIZE; // NEED to adjust for mappability and GC, adjusting for GC reduces bg_density... 
			// perhaps GC correction not applicable since it is based on fragment reads
			if (block0_start < 0) block0_start = BL[j]; // empty block... 
			flag++; 
		}
		i = j;
	}
        *bgdensity = bg_reads/bg_bins;  mean_fraglen /= fosmids; *mean_fl = mean_fraglen;

        fprintf(stderr,"fragments %d bgr:%0.0f/%0.0f rd(perblock):%0.5f %0.5f genome-cov %d/%0.4f\t",fosmids,bg_reads,bg_bins,*bgdensity,*bgdensity/BLOCK_SIZE,genome_covered,(double)genome_covered/blocks);
	fprintf(stderr,"mean_frag_len %0.2f \n\n",mean_fraglen); 
	for (i=0;i<50;i++) fprintf(stderr,"%d:%d ",i+1,delta_counts[i]); fprintf(stderr," delta-counts\n");
	for (i=1;i<100;i++) fprintf(stderr,"%d:%0.0f| ",i,length_pdf[i]); fprintf(stderr,"\n");
	return fosmids;
}

int print_fosmid_list(struct FOSMID* fosmidlist,int fosmids,struct alignedread** readlist,FRAGMENT* flist,VARIANT* varlist,struct BLOCK_READ* RL,int blocks,int* BL,double frag_penalty,double read_density,double length_pdf[],double mean_fl)
{
	int i=0,k=0,l=0,j=0; double ws=0;
	int fragments_in_window= 0; int total_length =0; int distance_next=0;
	int clusters = 0; int curr_cluster = 0,prev=0;
	int cstart=-1,cend =-1; 
	double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE);
	double bg_ll,seg_ll,mean[3]; int vblocks =0;
	double delta=0;

        struct FOSMID* newfosmidlist = calloc(2*fosmids,sizeof(struct FOSMID)); int nfosmids =0; // new fosmid list 
	int nf=0; int middle_frags=0;
	for (nf=0;nf<2*fosmids;nf++) newfosmidlist[nf].middle = '0'; 


	for (i=0;i<fosmids;i++) 
	{
		if (cstart <0) { cstart = fosmidlist[i].firstblock; cend = fosmidlist[i].lastblock; } 

		fosmidlist[i].cluster = curr_cluster; 
		k = fosmidlist[i].lastblock; 	ws = (k+1)*BLOCK_SIZE-fosmidlist[i].firstblock*BLOCK_SIZE;
		if (FOSMIDS != 2)
		{
			seg_ll =  fosmidlist[i].score_window;
			bg_ll =0; for (j=fosmidlist[i].fv;j<=fosmidlist[i].lv;j++) bg_ll += RL[BL[j]].bgscore; 
			delta = seg_ll+frag_penalty-bg_ll;
			if (delta < 5) fprintf(stdout,"small-:"); fprintf(stdout,"comparison seg:%0.2f bg:%0.2f delta %0.2f length %0.2f\n",seg_ll,bg_ll,delta,ws);
		}

		vblocks=0; for (j=fosmidlist[i].firstblock;j<fosmidlist[i].lastblock;j++) {  if (RL[j].valid == '1') vblocks++; } 
		print_block(RL,fosmidlist[i].firstblock,fosmidlist[i].lastblock,2);

		fprintf(stdout,"%d ==block %d...%d length %0.0f window:%0.2f 1 | ",fosmidlist[i].cluster,fosmidlist[i].firstblock*BLOCK_SIZE,k*BLOCK_SIZE,ws,fosmidlist[i].score_window);
		fprintf(stdout,"%d....%d reads %0.0f density(reads/bp) %0.3f mean(/kb) %0.3f val-blocks %0.3f delta %0.3f\n",fosmidlist[i].firstblock,k,fosmidlist[i].reads_window,(float)fosmidlist[i].ws/fosmidlist[i].reads_window,fosmidlist[i].mean*1000/BLOCK_SIZE,(float)vblocks/(fosmidlist[i].lastblock-fosmidlist[i].firstblock),delta);

		print_reads_window(readlist,fosmidlist[i].firstread,fosmidlist[i].lastread,flist,varlist,1);	
		// function to print fragment to file 
		generate_single_fragment(readlist,flist,varlist,&fosmidlist[i],"");

		if (i < fosmids-1) distance_next = (fosmidlist[i+1].firstblock-fosmidlist[i].lastblock)*BLOCK_SIZE; else distance_next = 0;
		fosmidlist[i].distance_next = distance_next; 
		fragments_in_window++; total_length += ws;

		if (i < fosmids-1) print_block(RL,fosmidlist[i].lastblock+1,fosmidlist[i+1].firstblock-1,1);

		if (distance_next < MIN_SEPARATION) 
		{
			fprintf(stdout,"....MERGE.... %d\n\n",distance_next); total_length += distance_next; cend = fosmidlist[i+1].lastblock;  
		}
		else // analyze cluster and print output
		{
			nf = nfosmids; 
			if ((fragments_in_window >=2 && fragments_in_window <= 4) || (fragments_in_window ==1 && total_length >= mean_fl ))
			{
				//find_2seg(RL,blocks,BL,cstart,cend,log_sizes,frag_penalty,read_density,length_pdf,&best_out);
				find_3seg(RL,BL,cstart,cend,log_sizes,frag_penalty,read_density,fragments_in_window,length_pdf,newfosmidlist,&nfosmids);
				
			}
			if (nf == nfosmids)
			{
				// copy the fragments in the cluster to the new list AS IS 
				for (nf=i-fragments_in_window+1;nf<=i;nf++) 
				{
					newfosmidlist[nfosmids].firstblock = fosmidlist[nf].firstblock; newfosmidlist[nfosmids].lastblock = fosmidlist[nf].lastblock; 
					newfosmidlist[nfosmids].start = fosmidlist[nf].start; newfosmidlist[nfosmids].end = fosmidlist[nf].end;
					newfosmidlist[nfosmids].firstread = fosmidlist[nf].firstread; newfosmidlist[nfosmids].lastread = fosmidlist[nf].lastread;
					newfosmidlist[nfosmids].mean  = fosmidlist[nf].mean; newfosmidlist[nfosmids].reads_window = fosmidlist[nf].reads_window; 
					newfosmidlist[nfosmids].ws =  fosmidlist[nf].ws; //newfosmidlist[nfosmids].score_window =fosmidlist[nf].score_window; 
					nfosmids++; 
				}
			}
			else 
			{
				fprintf(stdout,"modified original fragments %d..%d  new %d newlistlength %d\n",i,fragments_in_window,nfosmids-nf,nfosmids);
			}

			k = cend;
			prev = RL[k].previous[4]; 
			while (1)
			{
				l = (k-BL[prev+1]+1);
				if (l >=5) fprintf(stdout,"prev-4 %d-%d length:%d mean:%0.4f\n",BL[prev+1]*BLOCK_SIZE,k*BLOCK_SIZE,l*BLOCK_SIZE,RL[k].pmean[4]*1000/BLOCK_SIZE);
				if (prev <= 0 || BL[prev+1] <= cstart) break; 
				k = BL[prev]; prev = RL[k].previous[4];
			}

			fprintf(stdout,"\n-------- frags:%d ------ W:%d --------- %d ------------%d...%d\n\n\n\n\n",fragments_in_window,total_length,distance_next,cstart,cend);
			fragments_in_window = 0; total_length = 0; curr_cluster = i+1; clusters++; 
			cstart = fosmidlist[i+1].firstblock; cend = fosmidlist[i+1].lastblock;
		}
	}

	fprintf(stdout,"printing new fragments to file \n");
	for (i=0;i<nfosmids;i++) 
	{
		// print intervals for simulations, allow overlapping fragments... don't print middle fragments 
		if (newfosmidlist[i].middle == '0') fprintf(stdout,"BEDFILE %s\t%d\t%d\t%0.0f_%f\n",readlist[newfosmidlist[i].firstread]->chrom,newfosmidlist[i].start,newfosmidlist[i].end,newfosmidlist[i].ws,newfosmidlist[i].mean);	
		else if (newfosmidlist[i].middle == '2') fprintf(stdout,"BEDFILE %s\t%d\t%d\t%d_%f_OV\n",readlist[newfosmidlist[i].firstread]->chrom,newfosmidlist[i].start1,newfosmidlist[i].end1,newfosmidlist[i].end1-newfosmidlist[i].start1,newfosmidlist[i].mean);	
		if (newfosmidlist[i].middle == '1') generate_single_fragment(readlist,flist,varlist,&newfosmidlist[i],"1");
		else generate_single_fragment(readlist,flist,varlist,&newfosmidlist[i],"0");
	}

	return clusters;
}

int dynamic_programming_loop(struct BLOCK_READ* RL,int blocks,int* BL,int validblocks,double read_density,double block_penalty)
{
        int i = 0,j=0,k=0,eb=0,reads=0,neb=0;
        double ws=0,ws_corrected=0,flag=0;
        double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;
        double score0,lambda,score_partial,sum_partial,score1,frag_density;
        double log_sizes[BLOCK_SIZE]; for (i=1;i<BLOCK_SIZE;i++) log_sizes[i] = log(i)-log(BLOCK_SIZE); log_sizes[0]= -10000; // dummy
        double logRD = log(read_density);

        double block_open_penalty_1 = (block_penalty);
        double block_open_penalty_2 = 2*(block_penalty)/3;
        double fragment_p = block_open_penalty_1; double fragment_p_1 = log(1.0-exp(fragment_p)); 
	double gap_nb=0; // gap between adjacent bins that are valid 
	fprintf(stderr,"DP loop, rd %f block-open-penalty %f fragment_p %f \n",read_density,block_penalty,fragment_p);

	for (i=0;i<blocks;i++)
	{
        	for (j=0;j<5;j++) {RL[i].bscore[j] = -100000; RL[i].previous[j] = -1; }
		RL[i].pruned = '0'; 
	}

	i=0; // important to initialize
        score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size]-log(RL[BL[i]].GC_corr))-read_density*(float)RL[BL[i]].adjusted_size/(BLOCK_SIZE*RL[BL[i]].GC_corr);
        RL[BL[i]].bscore[0] = score0;  RL[BL[i]].bscore[4] =  RL[BL[i]].bscore[0]; RL[BL[i]].bgscore = score0; 
        RL[BL[i]].score_partial = RL[BL[i]].reads*(log_sizes[(int)RL[BL[i]].adjusted_size] -log(RL[BL[i]].GC_corr));
        RL[BL[i]].W = (float)RL[BL[i]].adjusted_size/(RL[BL[i]].GC_corr*BLOCK_SIZE);
		
	//fprintf(stdout,"i %d score best %f reads %f %f\n",BL[i],RL[BL[i]].bscore[0],RL[BL[i]].reads,score0);
	// keep track of # non-pruned blocks before block 'i', when we are going backwards from i---j, if there are no more non-pruned blocks, we can stop 

        for (i=1;i<validblocks;i++)
        {
                score0 = RL[BL[i]].reads*(logRD + log_sizes[(int)RL[BL[i]].adjusted_size]-log(RL[BL[i]].GC_corr))-read_density*(float)RL[BL[i]].adjusted_size/(BLOCK_SIZE*RL[BL[i]].GC_corr);
		RL[BL[i]].bgscore = score0;
                // score of empty segment...  no penalty for non-mappability ?? 

                RL[BL[i]].bscore[0] = RL[BL[i-1]].bscore[0] + score0; // + fragment_p_1 is ln(1-p) where 'p' is probability of selecting bin as fragment start 
                RL[BL[i]].bscore[4] = RL[BL[i]].bscore[0];
                // does not matter if previous block was '1' or '0' (no penalty for switching from 1 -> 0), only penalty from 0 -> 1 for opening new block 
                RL[BL[i]].previous[0] = i-1; RL[BL[i]].previous_cluster = RL[BL[i-1]].previous_cluster; RL[BL[i]].score_window = score0;
                RL[BL[i]].previous[4] = i-1;
				

                // mean is the average # of reads per 100 bp window = \lambda
                RL[BL[i]].score_partial = RL[BL[i]].reads*(log_sizes[(int)RL[BL[i]].adjusted_size] -log(RL[BL[i]].GC_corr));
                RL[BL[i]].W = (float)RL[BL[i]].adjusted_size/(RL[BL[i]].GC_corr*BLOCK_SIZE);

                variance = RL[BL[i]].reads*RL[BL[i]].reads*(float)RL[BL[i]].adjusted_size/BLOCK_SIZE;
                score_partial = RL[BL[i]].score_partial; mean = RL[BL[i]].reads; W = RL[BL[i]].W; reads = RL[BL[i]].reads;
	
                j = i-1; eb = 0; ws =0; ws_corrected=0; neb =0;
                // consider segments that start at 'j' and ends at 'i'  (j decreasing) 
                while ( (BL[i]-BL[j]) < (int)(MAX_FRAG_LENGTH/BLOCK_SIZE) && BL[j] > 0 && eb*BLOCK_SIZE < max_empty_stretch) // segfault if we set BL[j] >= 0 
                {
			if ((BL[j+1]-BL[j])*BLOCK_SIZE > 10000) break;  // do not link blocks that are separated by long gap of low mappability bins 
                        if (RL[BL[j]].reads > 0) { neb++; eb = 0; }  else eb++;  

			reads += RL[BL[j]].reads; mean += RL[BL[j]].reads; W += RL[BL[j]].W; score_partial += RL[BL[j]].score_partial;
                        variance += RL[BL[j]].reads*RL[BL[j]].reads*(float)RL[BL[j]].adjusted_size/BLOCK_SIZE;
                        ws_corrected += (float)RL[BL[j]].adjusted_size; ws += BLOCK_SIZE;

                        // default is 10 reads, 1 kb and 0.1= 1 read per 1kb | threshold should be based on background read 2-5 fold density...
			// removed filter on min-reads 10/28/16, reads >= MIN_READS_ gives 0 fragments
                        if (RL[BL[j-1]].pruned == '0' && ws_corrected >= MIN_FRAG_LENGTH)// && mean*BLOCK_SIZE/W >= MIN_READ_DENSITY) // extra conditions to consider it a valid segment...
                        {
				blockscore = score_partial + reads*log(mean/W) - mean;  // score of the full block 

                                if (j ==0) score0 = blockscore + block_open_penalty_1; // +prior_size + gap_length_penalty; 
                                else score0 = RL[BL[j-1]].bscore[0] + blockscore + block_open_penalty_1; // + prior_size + gap_length_penalty;
                                if (score0 > RL[BL[i]].bscore[0])
                                {
                                        RL[BL[i]].bscore[0] = score0; RL[BL[i]].previous[0] = (j-1); RL[BL[i]].reads_window = reads;
                                        RL[BL[i]].score_window = blockscore;    RL[BL[i]].previous_cluster = i;
                                        RL[BL[i]].mean = mean/W; RL[BL[i]].variance = variance/W - RL[BL[i]].mean*RL[BL[i]].mean;
                                }
                                RL[BL[j]].curr_score = RL[BL[j-1]].bscore[0] + blockscore;
				
                                if (j ==0) score0 = blockscore + block_open_penalty_2;
                                else score0 = RL[BL[j-1]].bscore[4] + blockscore + block_open_penalty_2; // + prior_size + gap_length_penalty;
                                if (score0 > RL[BL[i]].bscore[4])
                                {
                                        RL[BL[i]].bscore[4] = score0; RL[BL[i]].previous[4] = (j-1); RL[BL[i]].pmean[4] = mean/W;
                                }
				if (j > 0) RL[BL[j-1]].tscore = blockscore; 
				
                        }
			else if (j > 0) RL[BL[j-1]].tscore = 1; 
                        j--;
                }

        	 k=j; j=i-1; while (j > k) // prune blocks using Kirill's criteria, sub-optimal breakpoints
                {
                        if (j > 0 && RL[BL[j-1]].tscore < 0 &&  (RL[BL[j-1]].tscore + RL[BL[j-1]].bscore[0]  < RL[BL[i]].bscore[0]) && (RL[BL[j-1]].tscore + RL[BL[j-1]].bscore[4]  < RL[BL[i]].bscore[4]) )
                        {
                                if (RL[BL[j-1]].pruned == '0') RL[BL[j-1]].pruned = '1'; // fprintf(stdout,"pruning bp %d... %d \n",BL[j-1],BL[i]); 
                        }
                        j--;
                }
                if (i%500000 ==0) fprintf(stderr,"done for block %d %d\n",i,blocks);
        }
}


// function called on read list from single chromosome....
int process_reads_chrom(struct alignedread** readlist,int s,int e,FRAGMENT* flist,VARIANT* varlist,REFLIST* reflist,REFLIST* genomemask )
{
	find_matepair(readlist,s,e,2000); fprintf(stderr,"start s:%d end read e:%d \n",s,e);//distance_nextread_dist(readlist,s, e); 

	int reads=0,i=0,j=0,k=0,empty=0,non_empty=0,eb=0,t=0;
	double blockscore=0,prior_size=0,mean=0,variance=0,W=0,mean_pb;

	int blocks = readlist[e-1]->position/BLOCK_SIZE+1; if (blocks < 100) return -1; // small number of blocks 
	struct BLOCK_READ* RL = calloc((blocks+1),sizeof(struct BLOCK_READ)); 
	if (SINGLE_READS ==1) { for (i=0;i<blocks;i++) RL[i].subcounts = calloc(BS1*2,sizeof(uint8_t)); }// added 10/21/16
        else { for (i=0;i<blocks;i++) RL[i].subcounts = calloc(BS1,sizeof(uint8_t)); }

	init_blocklist(RL,blocks,readlist,s,e);

	int validblocks = calculate_block_stats(RL,blocks,reflist,genomemask); // calculate GC and mappability of each window, smoothed read-depth
	if (validblocks < 100) { free(RL); return -1; } // small number of blocks 
	int* BL = calloc(validblocks,sizeof(int)); j=0;
	for (i=0;i<blocks;i++) 
	{
		if (RL[i].valid == '0') continue; BL[j] = i; j++; 
	}
	
	double estimates[3];  estimate_lambda_b(RL,blocks,estimates); // use exponential distribution to estimate mean of background read density
	double read_density = (double)BLOCK_SIZE/(estimates[1]+0.0001); read_density *= 1.5;
	if (read_density < 1.0e-8 || read_density > 1)  read_density = 0.01; 

	double block_open_penalty = -7, length_penalty = -4.61; // -4.61 is penalty for segment length = 1/100 uniform distribution over 1-100 kb 
	//double block_open_penalty = 0.5*log(1.0/validblocks), length_penalty = 0.5*log(1.0/100)  + 0.5*log(1.0/100); // length penalty also has sum over density penalty... 
	double block_penalty = block_open_penalty + length_penalty; 
	double lf = 0.001;

	fprintf(stderr,"block open penalty = %0.2f length penalty %0.2f bg-density %f\n",block_open_penalty,length_penalty,read_density); 
	double length_pdf[300]; // 1 kb intervals counts, [0] is sum over all fragments 
	double length_pdf_log[300]; 
	for (i=0;i<300;i++) length_pdf[i] = 0; 

	int refflag = 0;	if (reflist->current >=0) refflag = 1;
	struct FOSMID* fosmidlist = calloc(100000,sizeof(struct FOSMID)); 
	double new_read_density=read_density; int fosmids=0; int iter =0; double mean_fl;
	//fprintf(stdout,"\nDP iter reads %d bg-read-density (per block) %f \n\n",e-s+1,read_density);
	iter=0; while (iter++ < 10) 
	{
		// call function to DP backtracking to find optimal segmentation
		//fprintf(stderr,"DP iter reads %d...%d %d bg-read-density %f %f ",readlist[s]->position,readlist[e-1]->position,e-s+1,read_density_global,read_density);
		dynamic_programming_loop(RL,blocks,BL,validblocks,read_density,block_penalty);
		fosmids =  find_best_segmentation(RL,blocks,BL,validblocks,fosmidlist,refflag,&new_read_density,length_pdf,block_penalty,&mean_fl);
		lf = ((float)fosmids+1.0)/validblocks; 
		block_open_penalty = log(lf)-log(1.0-lf); block_penalty = block_open_penalty + length_penalty;
		if ( fabsf(read_density-new_read_density) < 0.0002) 
		{		
			fprintf(stderr,"final iter convergence %f %f # of iters %d \n",read_density,new_read_density,iter);
			read_density = new_read_density;  
			break; 
		}
		else fprintf(stderr,"iter %d newblockpenalty %f\n",iter+1,block_open_penalty);
		read_density = new_read_density;
	}
	for (i=0;i<300;i++) 
	{
		if (length_pdf[i]  < 2.001) length_pdf[i] = 0.1; length_pdf_log[i] = log(length_pdf[i]) -log(length_pdf[0]); 
	}

	// print list of fragments identified, print fragments that can be merged and sum of lengths...
	qsort(fosmidlist,fosmids,sizeof(struct FOSMID),cmp_fosmid);
	int clusters = print_fosmid_list(fosmidlist, fosmids, readlist, flist,varlist,RL,blocks,BL,block_penalty,read_density,length_pdf_log,mean_fl);
	fprintf(stdout,"statistics for chrom %s fragments %d clusters %d blocks %d validblocks %d %0.4f %0.5f L %0.1f\n",readlist[s]->chrom,fosmids,clusters,blocks,validblocks,(float)validblocks/blocks,read_density*1000/BLOCK_SIZE,mean_fl);
	fprintf(stderr,"statistics for chrom %s fragments %d clusters %d blocks %d validblocks %d %0.4f %0.5f L %0.1f\n",readlist[s]->chrom,fosmids,clusters,blocks,validblocks,(float)validblocks/blocks,read_density*1000/BLOCK_SIZE,mean_fl);
	for (i=0;i<110;i++) fprintf(stderr,"#"); fprintf(stderr,"\n\n");

	for (i=0;i<blocks;i++) free(RL[i].subcounts); 
	free(fosmidlist); free(RL); free(BL);
	return 1;
}

//calculate_GC_stats(RL,blocks,fosmidlist,fosmids); // this is not used to adjust read-depth but can be used in 2nd-pass
