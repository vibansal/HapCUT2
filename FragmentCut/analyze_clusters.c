// likelihood of segment is e^(-l)l^r/r! summed over bins here 'l' is the mean read density (can be user provided) | r! is constant and can be ignored...
// L(s,e,l) = r_i. log(l.correction) - l.correction
// RL[i].W = adjusted_size(mappability)/(GC)*BLOCK_SIZE, ideally this should be 1 but can be less or more 


// this function is not correct since it does not use GC.correction 
// these functions not used elsewhere, try with counts ??
int use_counts=0; // not use smoothed read counts but full counts

// k! is missing -> we get both positive and negative likelihood values...
// RL[i].W incorporates GC correction 
double segment_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* LOG_sizes,double PARAMS[]) // first block .... last block
{
	int i=0; double mean=0,W=0,score_partial=0,reads=0,segment_score=0,sum_partial=0; 
	for (i=fb;i<lb;i++)
	{
		if (RL[i].valid != '0') 
		{ 
			W += RL[i].W; 
			if (use_counts ==1) reads += RL[i].counts;  else reads += RL[i].reads; 
			segment_score += RL[i].score_partial; 
		} 
	}
	mean = reads/W; segment_score += reads*log(mean)-mean*W; // same as below
	PARAMS[0] = reads; PARAMS[1] = W; PARAMS[2] = mean; 
	return segment_score; 
	//fprintf(stderr,"segscore %f mean %f reads %f W %f \n",segment_score,mean,reads,W);
}

// score of segment assuming no fragment overlaps
double background_likelihood(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* LOG_sizes,double read_density)
{
	double LOG_RD = log(read_density);
	int i=0; double segment_score=0, reads=0,W=0;
	for (i=fb;i<lb;i++)
	{
		if (RL[i].valid != '0' && use_counts ==1) segment_score += RL[i].counts*(LOG_RD + LOG_sizes[(int)RL[i].adjusted_size]-(RL[i].GC_corr_log))-read_density*(float)RL[i].adjusted_size/(BLOCK_SIZE*RL[i].GC_corr);
		else if (RL[i].valid != '0' && use_counts ==0) segment_score += RL[i].reads*(LOG_RD + LOG_sizes[(int)RL[i].adjusted_size]-(RL[i].GC_corr_log))-read_density*(float)RL[i].adjusted_size/(BLOCK_SIZE*RL[i].GC_corr);
		//if (RL[i].valid != '0') { W += RL[i].W; reads += RL[i].reads;  segment_score += RL[i].score_partial; } 
	}
	//segment_score += reads*log(read_density)-read_density*W;
	return segment_score;
}

// expected  sum of two read-density | union instead of sum | optimized oct29, 2016, big difference in speed...
double sum_coverages(double m1,double m2)
{
	int i=0,s=BLOCK_SIZE/BS1; double ne=0;
	float a = m1/s; if (a > 1) a= 1; 
	float b = m2/s; if (b > 1) b = 1; 
	float un = (a + b-a*b)*s; return un; 
}


// given start and end as index into RL[], add interval to list of fragments 
int add_fragment_list(struct BLOCK_READ* RL,struct FOSMID* fosmidlist,int fosmids,int f,int l,double mean)
{
        //while (RL[f].valid == '0') f++; while (RL[l].valid == '0' && l > 0) l--;  // is this step needed

        int k=0;
        int first = RL[f].firstread; k=1; while(first < 0) { first = RL[f+k].firstread; k++; }
        int last = RL[l].lastread; k= 1;  while(last < 0) { last = RL[l-k].lastread; k++; }

        fosmidlist[fosmids].firstblock = f; fosmidlist[fosmids].lastblock = l;
        fosmidlist[fosmids].start = f*BLOCK_SIZE; fosmidlist[fosmids].end = (l+1)*BLOCK_SIZE; fosmidlist[fosmids].ws =  (l+1-f)*BLOCK_SIZE;
        fosmidlist[fosmids].firstread = first; fosmidlist[fosmids].lastread = last;
        fosmidlist[fosmids].mean  = mean;
        fosmidlist[fosmids].reads_window =0; for (k=f;k<l;k++) fosmidlist[fosmids].reads_window += RL[k].reads;
        //fosmidlist[fosmids].score_window = RL[BL[i]].score_window; // fv, lv, score_window variables not needed
}

// N=1, N=2,non-overlap/overlap
// best segmentation with 2 overlapping fragments, 1 single fragment or with 2 non-overlapping fragments
double find_2seg(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double frag_penalty,double read_density,int frags,double length_pdf_LOG[],int* ret,struct FOSMID* fosmidlist,int* fsize)
{
	double mean1[3]={0,0,0},mean2[3]={0,0,0},mean3[3]={0,0,0}; double logRD = log(read_density);
	int i=0,j=0,k=0;
	double l1,l2,l3,l2_1,l,l2_empty;
	double bestseg[7] = {-10000000,0,0,0,0,0,0};
	double best2seg[7] = {-10000000,0,0,0,0,0,0};
	double best3seg[7] = {-10000000,0,0,0,0,0,0};
        double best1seg[7] = {-10000000,0,0,0,0,0,0};

	double middlem = 0;
	int min_sl = (int)(2000/BLOCK_SIZE); // minimum segment length units of 100 bp 
	int s1 =0,s2=0; double p1=0,p2=0;

	// also consider the possibility of 3-segments independent... compromise when 2 overlapping frags don't work...

	for (i=fb+min_sl;i<lb-min_sl*2;i+=2)
	{
		if (RL[i].valid == '0') continue; 
		j = i+1; 
		l1 = 	segment_likelihood(RL,BL,fb,i,log_sizes,mean1); 
		l2_empty = background_likelihood(RL,BL,i,j,log_sizes,read_density);  // score of segment with no fragment 
		l2 = 	segment_likelihood(RL,BL,i,j,log_sizes,mean2);  l3 = 	segment_likelihood(RL,BL,j,lb,log_sizes,mean3);
		// further speedup by calculating updates

		while (j<lb-min_sl)
		{
			if (RL[j].valid != '0')
			{
				middlem= sum_coverages(mean1[2],mean3[2]);
				l2_1 = l2 - (mean2[0]*log(mean2[2]) - mean2[0]) + mean2[0]*log(middlem)-mean2[0]*(middlem)/mean2[2]; // using a + b as mean for middle segment
				l = l1+l3 + l2_1 + 2*frag_penalty;
				s1 = (int)(((j-fb)*BLOCK_SIZE)/1000)+1; s2 = (int)(((lb-i)*BLOCK_SIZE)/1000)+1; 
				if (s1 < 200) p1 = length_pdf_LOG[s1]; else p1= 0; 
				if (s2 < 200) p2 = length_pdf_LOG[s2]; else p2= 0; 
				l += (p1)+(p2); // add fragment penalty

				if (mean2[2] > mean1[2] && mean2[2] > mean3[2] && l > bestseg[0] && j-i >= min_sl) 
				{ 
					bestseg[0] = l; bestseg[1] = i; bestseg[2]= j; bestseg[3] = mean1[2]; bestseg[4] = mean3[2];bestseg[6] = middlem; bestseg[5] = mean2[2]; 
				} 
				//fprintf(stdout,"triple-seg %d...%d %f %f %f\n",i,j,l1,l2_1,l3);

				// middle segment is background, not part of fragment 
				l = l1 + l2_empty + l3 + 2*frag_penalty; 
				s1 = (int)(((i-fb)*BLOCK_SIZE)/1000)+1; s2 = (int)(((lb-j)*BLOCK_SIZE)/1000)+1; 
				if (s1 < 200) p1 = length_pdf_LOG[s1]; else p1= 0; 
				if (s2 < 200) p2 = length_pdf_LOG[s2]; else p2= 0; 
				l += (p1)+(p2); // add fragment penalty
				if (l > best2seg[0]) // two non-overlapping fragments 
				{
					best2seg[0] = l; best2seg[1] = i; best2seg[2]= j; best2seg[3] = mean1[2]; best2seg[4] = mean3[2]; best2seg[5] = mean2[2]; 
				}

				// three non-overlapping segments
				s1 = (int)(((j-i)*BLOCK_SIZE)/1000)+1; 
				if (s1 < 200) p1 = length_pdf_LOG[s1]; else p1= 0; 
				l += l2 -l2_empty + frag_penalty + (p1); 
				if (l > best3seg[0]) // two non-overlapping fragments 
				{
					best3seg[0] = l; best3seg[1] = i; best3seg[2]= j; best3seg[3] = mean1[2]; best3seg[4] = mean3[2]; best3seg[5] = mean2[2]; 
				}

				// update l2,l3, mean2[x], mean3[x]  |  j-th bin is being added to segment2 and removed from segment3
				l2 -= mean2[0]*log(mean2[2])-mean2[0]; 
				mean2[1] += RL[j].W; mean2[0] += RL[j].reads;  mean2[2] = mean2[0]/mean2[1]; 
				l2 += log(RL[j].W)*RL[j].reads; l2 += mean2[0]*log(mean2[2])-mean2[0];

				l3 -= mean3[0]*log(mean3[2])-mean3[0];
				mean3[1] -= RL[j].W; mean3[0] -= RL[j].reads;  mean3[2] = mean3[0]/mean3[1]; 
				l3 -= log(RL[j].W)*RL[j].reads; l3 += mean3[0]*log(mean3[2])-mean3[0];
	
				l2_empty += RL[j].reads*(logRD + log_sizes[(int)RL[j].adjusted_size]-(RL[j].GC_corr_log))-read_density*(float)RL[j].adjusted_size/(BLOCK_SIZE*RL[j].GC_corr);
			}
			j +=1; 
		}
			
	}

	// N=1
	s1 = (int)(((lb-fb)*BLOCK_SIZE)/1000)+1;  if (s1 < 200) p1 = length_pdf_LOG[s1]; else p1= 0; 
	fprintf(stdout,"%d s1 %d p1 %f \n",(lb-fb)*BLOCK_SIZE,s1,p1);
	double ls = segment_likelihood( RL,BL,fb,lb,log_sizes,mean2)+frag_penalty; fprintf(stdout,"single-fragment %0.2f %d...%d LL %0.2f LL-p %0.2f\n",mean2[2],fb,lb,ls,ls+(p1));
	ls += (p1); best1seg[0] = ls; best1seg[3] = mean2[2]; 

	// N=2, no overlap 
	fprintf(stdout,"two-NO-fragments %d....%d....%d....%d LL %0.2f | ",fb,(int)best2seg[1],(int)best2seg[2],lb,best2seg[0]);
	fprintf(stdout,"intervals %0.0f:%0.2f %0.0f:%0.2f %0.0f:%0.2f \n",BLOCK_SIZE*(best2seg[1]-fb),best2seg[3],BLOCK_SIZE*(best2seg[2]-best2seg[1]),best2seg[5],BLOCK_SIZE*(lb-best2seg[2]),best2seg[4]);
	// N=2, overlapping by at least 10*BLOCK_SIZE
	fprintf(stdout,"two-OV-fragments %d....%d....%d....%d LL %0.2f | ",fb,(int)bestseg[1],(int)bestseg[2],lb,bestseg[0]);
	fprintf(stdout,"3-NO-fragments LL %0.2f | ",best3seg[0]);
	fprintf(stdout,"intervals %0.0f:%0.2f %0.0f:%0.2f:middle:%0.2f %0.0f:%0.2f \n",BLOCK_SIZE*(bestseg[1]-fb),bestseg[3],BLOCK_SIZE*(bestseg[2]-bestseg[1]),bestseg[6],bestseg[5],BLOCK_SIZE*(lb-bestseg[2]),bestseg[4]);
	int fosmids = *fsize;
	double chi2d = 3.32; 

	if (best3seg[0] > best2seg[0] + chi2d && best3seg[0] > bestseg[0] + 2*chi2d && best3seg[0] > best1seg[0] + 2*chi2d) // 3 segmentation
	{
		*ret = -1; return best3seg[0];  // no change to list
	}
	else if (best2seg[0] > bestseg[0]-chi2d && best2seg[0] > best1seg[0] + chi2d) // 2 non-overlapping fragments 
	{ 
		if (frags != 2) 
		{
			fprintf(stdout,"new fragment %d....%d %d %f \n",fb*BLOCK_SIZE,(int)best2seg[1]*BLOCK_SIZE,(int)best2seg[1]*BLOCK_SIZE-fb*BLOCK_SIZE,best2seg[3]);
			fprintf(stdout,"new fragment %d....%d %d %f \n",(int)best2seg[2]*BLOCK_SIZE,lb*BLOCK_SIZE,lb*BLOCK_SIZE-(int)best2seg[2]*BLOCK_SIZE,best2seg[4]);
			add_fragment_list(RL,fosmidlist,fosmids,fb,best2seg[1],best2seg[3]); fosmids++;
			add_fragment_list(RL,fosmidlist,fosmids,best2seg[2],lb,best2seg[4]); fosmids++; 
			*fsize +=2; 
		}
		*ret =1; return best2seg[0]; 
	}
	else if (bestseg[0] > best2seg[0] + chi2d && bestseg[0] > best1seg[0] + chi2d) // 2 overlapping fragments 
	{ 
		fprintf(stdout,"new fragment %d....%d %d %f \n",fb*BLOCK_SIZE,(int)bestseg[1]*BLOCK_SIZE,(int)bestseg[1]*BLOCK_SIZE-fb*BLOCK_SIZE,bestseg[3]);
		fprintf(stdout,"mid fragment %d....%d %d %f \n",((int)bestseg[1]+1)*BLOCK_SIZE,(int)bestseg[2]*BLOCK_SIZE,(int)bestseg[2]*BLOCK_SIZE-(int)bestseg[1]*BLOCK_SIZE,bestseg[6]);
		fprintf(stdout,"new fragment %d....%d %d %f \n",(int)bestseg[2]*BLOCK_SIZE,lb*BLOCK_SIZE,lb*BLOCK_SIZE-(int)bestseg[2]*BLOCK_SIZE,bestseg[4]);
		add_fragment_list(RL,fosmidlist,fosmids,fb,bestseg[1],bestseg[3]); 
		fosmidlist[fosmids].start1 = fb*BLOCK_SIZE;	fosmidlist[fosmids].end1 = (int)bestseg[2]*BLOCK_SIZE; fosmidlist[fosmids].middle = '2'; fosmids++;
		add_fragment_list(RL,fosmidlist,fosmids,bestseg[1]+1,bestseg[2]-1,bestseg[6]); fosmidlist[fosmids].middle ='1';  fosmids++; 
		add_fragment_list(RL,fosmidlist,fosmids,bestseg[2],lb,bestseg[4]); 
		fosmidlist[fosmids].end1 = lb*BLOCK_SIZE;	fosmidlist[fosmids].start1 = (int)bestseg[1]*BLOCK_SIZE; fosmidlist[fosmids].middle = '2';  fosmids++;
		*fsize += 3; 
		*ret = 2; return bestseg[0]; 
	} 
	else if (frags > 1 && frags < 3 && best1seg[0] > best2seg[0]-chi2d && best1seg[0] > bestseg[0]-chi2d ) // single fragment 
	{
		fprintf(stdout,"new fragment %d....%d %d %f \n",fb*BLOCK_SIZE,lb*BLOCK_SIZE,lb*BLOCK_SIZE-fb*BLOCK_SIZE,best1seg[3]);
		add_fragment_list(RL,fosmidlist,fosmids,fb,lb,best1seg[3]); fosmids++; *fsize +=1; 
		 *ret = 0; return best1seg[0];
	}
	else // nochange to list 
	{
		*ret = -1; return -1; 
	}
}


// 2 + 1 segmentation | N=3,non-overlap, N=3,1+2 should work...
double find_3seg(struct BLOCK_READ* RL, int* BL,int fb,int lb,double* log_sizes,double frag_penalty,double read_density,int frags,double length_pdf_log[],struct FOSMID* flist,int* fsize)
{
	// first find list of all potential breakpoints using weak penalty...  int* bkpts = calloc(sizeof(int),100); 
	int i=0; int k = lb; int prev = RL[lb].previous[0];
	int* bkpts = calloc(sizeof(int),frags*4); int b = 0; 

	while (1)
	{
		if (k-BL[prev] >=5 && BL[prev+1] > fb) bkpts[b++] = BL[prev]; 
		if (prev <= 0 || BL[prev+1] <= fb) break;
		k = BL[prev]; prev = RL[k].previous[0];
	}
	for (i=0;i<b;i++) fprintf(stdout,"bkpts %d:%d ",i,bkpts[i]); fprintf(stdout,"b %d\n",b);

	double mean[3]={0,0,0}; double ls1,ls2;
	int best=0;

	if (frags <=4) 
	{
		ls2 = find_2seg(RL,BL,fb,lb,log_sizes,frag_penalty,read_density,frags,length_pdf_log,&best,flist,fsize);
		if (best ==2) fprintf(stdout,"best-overlap %0.2f \n",ls2);
		else if (best ==1) fprintf(stdout,"best-2frags %0.2f \n",ls2);
		else if (best ==0 && frags > 1) fprintf(stdout,"merge-pair %0.2f \n",ls2);
		// if the 2-way segmentation is much better, don't call the next loop
	}

	/*
	if (frags >=3544) // for doing 1+2 segmentation, disable for now, uses prior list of breakpoints to determine single fragment interval from start or end... 
	{
		fprintf(stdout,"\n");
		ls2 = find_2seg(RL,BL,fb,bkpts[0],log_sizes,frag_penalty,read_density,frags,length_pdf_log,&best,flist,fsize); 
		ls1 = segment_likelihood( RL,BL,bkpts[0],lb,log_sizes,mean)+frag_penalty;
		fprintf(stdout,"------ 2overlap+1 ll %0.2f %0.2f sum %0.2f-----------\n \n",ls1,ls2,ls1+ls2); 

		fprintf(stdout,"\n");
		ls1 = segment_likelihood( RL,BL,fb,bkpts[b-1],log_sizes,mean)+frag_penalty;
		ls2 = find_2seg(RL,BL,bkpts[b-1],lb,log_sizes,frag_penalty,read_density,frags,length_pdf_log,&best,flist,fsize); 
		fprintf(stdout,"-----1+2overlap ll %0.2f %0.2f sum %0.2f---------------\n \n",ls1,ls2,ls1+ls2); 
	}	
	*/
	free(bkpts);

}
