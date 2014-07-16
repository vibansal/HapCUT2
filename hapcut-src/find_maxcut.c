

/****** CODE TO FIND MAX CUT FOR EACH COMPONENT *************/

void evaluate_cut_component(struct fragment* Flist,struct SNPfrags* snpfrag,struct BLOCK* clist, int k,int* slist, char* HAP1,int iter);
float compute_goodcut(struct SNPfrags* snpfrag,char* hap,int* slist,int N,struct fragment* Flist, int algo);

/**************** DETERMINISTIC MAX_CUT MEC IMPLEMENTATION *********************************************************//////

void evaluate_cut_component(struct fragment* Flist,struct SNPfrags* snpfrag,struct BLOCK* clist, int k,int* slist, char* HAP1,int iter)
{
	int i=0,j=0,t=0; int f=0; 
	float cutvalue; float newscore;
	/*
	i=0;for (j=clist[k].offset;j<clist[k].offset+clist[k].length;j++) 
	{
		if (snpfrag[clist[k].offset].component == snpfrag[j].component) { slist[i] = j; i++; } 
	}*/
	// slist is copied from clist[k].slist
	for (t=0;t<clist[k].phased;t++) slist[t] = clist[k].slist[t]; // new way to determine slist

	// not required but we do it to avoid errors
	clist[k].MEC =0; 
	for (i=0;i<clist[k].frags;i++) 
	{
		update_fragscore(Flist,clist[k].flist[i],HAP1); clist[k].MEC += Flist[clist[k].flist[i]].currscore;
	}

	clist[k].bestMEC = clist[k].MEC; newscore= clist[k].MEC; 
	// evaluate the impact of flipping of each SNP on MEC score
	// this loop is essentially a N^2 time hog, do we need it, april 17 2012
	for (t=0;t<clist[k].phased;t++)
	{
		if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0'; 	else if (HAP1[slist[t]] == '0') HAP1[slist[t]] = '1'; 
		
		for (i=0;i<snpfrag[slist[t]].frags;i++) 
		{
			f = snpfrag[slist[t]].flist[i]; 	newscore -= Flist[f].currscore;
			update_fragscore(Flist,f,HAP1); 	newscore += Flist[f].currscore;
		}
		/*
		fprintf(stdout,"%f ",clist[k].MEC);
		fprintf(stdout,"flist %d %d %f %f \n",slist[t],snpfrag[slist[t]].frags,newscore,clist[k].MEC);
		*/

		if (newscore > clist[k].bestMEC)  // newMEC is not better than original MEC so we need to revert back to oldMEC for each fragment 
		{
			if (HAP1[slist[t]] == '1') HAP1[slist[t]] = '0';		else if ( HAP1[slist[t]] == '0') HAP1[slist[t]] = '1';
			//clist[k].MEC = clist[k].bestMEC; 
			for (i=0;i<snpfrag[slist[t]].frags;i++) 
			{
				f = snpfrag[slist[t]].flist[i]; 	newscore -= Flist[f].currscore;
				update_fragscore(Flist,f,HAP1); 	newscore += Flist[f].currscore;
			}
		}
		else
		{
			clist[k].MEC = newscore; clist[k].bestMEC = newscore;
		}
	}
	
	newscore = clist[k].MEC; 
	clist[k].MEC =0; for (i=0;i<clist[k].frags;i++) 
	{
		update_fragscore(Flist,clist[k].flist[i],HAP1); clist[k].MEC += Flist[clist[k].flist[i]].currscore;
	}
	clist[k].bestMEC = clist[k].MEC;
	//if (fabsf(newscore-clist[k].MEC) >= 0.001) fprintf(stdout,"old %f %f MECcheck\n",newscore,clist[k].MEC);

	cutvalue =10; if (clist[k].MEC > 0) cutvalue = compute_goodcut(snpfrag,HAP1,slist,clist[k].phased,Flist,MINCUTALGO);
	// flip the subset of columns in slist with positive value 
	if (cutvalue <= 3 || MINCUTALGO ==2) 
	{	 //getchar();
		//printf("code reached here \n");
		for (i=0;i<clist[k].phased;i++) 
		{
			if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
			else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
		}
		clist[k].bestMEC = clist[k].MEC; 
		clist[k].MEC =0; 
		for (i=0;i<clist[k].frags;i++) 
		{
			update_fragscore(Flist,clist[k].flist[i],HAP1); clist[k].MEC += Flist[clist[k].flist[i]].currscore;
		}
		if (clist[k].MEC > clist[k].bestMEC) 
		{
			for (i=0;i<clist[k].phased;i++)
			{
				if (slist[i] > 0 && HAP1[slist[i]] == '1') HAP1[slist[i]] = '0';
				else if (slist[i] > 0 && HAP1[slist[i]] == '0') HAP1[slist[i]] = '1';
			} 
			clist[k].MEC =0; 
			for (i=0;i<clist[k].frags;i++) 
			{
				update_fragscore(Flist,clist[k].flist[i],HAP1); clist[k].MEC += Flist[clist[k].flist[i]].currscore;
			}
		}
		else clist[k].bestMEC = clist[k].MEC;
	} 
	if (iter > 0 && clist[k].MEC > 0) fprintf(stdout,"component %d offset %d length %d phased %d  calls %d MEC %0.1f cutvalue %f bestMEC %0.2f\n",k,clist[k].offset,clist[k].length,clist[k].phased,clist[k].calls,clist[k].MEC,cutvalue,clist[k].bestMEC);

	// print out cut information
	if (iter >= 2 && clist[k].MEC > clist[k].lastMEC && iter%2==0 && clist[k].phased < 50 && clist[k].MEC > 0 ) 
	{
		/*
		for (i=0;i<clist[k].length;i++) 
		{ 
			if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component) fprintf(stdout,"%c",HAP1[clist[k].offset+i]); else fprintf(stdout,"-"); 
		} fprintf(stdout,"\n");
		for (i=0;i<clist[k].length;i++) 
		{ 
			if (snpfrag[clist[k].offset+i].component == snpfrag[clist[k].offset].component && slist[i] > 0) fprintf(stdout,"c"); else fprintf(stdout,"-"); 
		} fprintf(stdout,"\n");
		*/
	}
}

/********* THIS IS THE MAIN FUNCTION FOR HAPLOTYPE ASSEMBLY USING ITERATIVE MAX CUT computations **************/

float compute_goodcut(struct SNPfrags* snpfrag,char* hap,int* slist,int N,struct fragment* Flist, int algo)
{
	// given a haplotype 'hap' and a fragment matrix, find a cut with positive score 
	int totaledges=0,i=0,j=0,k=0,l=0;  int wf = 0; //if (drand48() < 0.5) wf=1;
	float W = 0;

	/* CODE TO set up the read-haplotype consistency graph */
	for (i=0;i<N;i++)
	{
		snpfrag[slist[i]].tedges=0; k=-1;
		// edges contain duplicates in sorted order, but tedges is unique count of edges
		for (j=0;j<snpfrag[slist[i]].edges;j++) 
		{
			if (k != snpfrag[slist[i]].elist[j].snp) { snpfrag[slist[i]].tedges++; k = snpfrag[slist[i]].elist[j].snp; } 
		}
	}
	for (i=0;i<N;i++)
	{
		snpfrag[slist[i]].tedges=0; k=-1;
		for (j=0;j<snpfrag[slist[i]].edges;j++) 
		{
			if (k != snpfrag[slist[i]].elist[j].snp) 
			{ 
				snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].snp = snpfrag[slist[i]].elist[j].snp; 
				k = snpfrag[slist[i]].elist[j].snp; 
				W = (float)edge_weight(hap,slist[i],k,snpfrag[slist[i]].elist[j].p,Flist,snpfrag[slist[i]].elist[j].frag);
				if (wf ==0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls-1; //(fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1);	
				snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges].w = W; 
				snpfrag[slist[i]].tedges++;  totaledges++;
			} 
			else if (k == snpfrag[slist[i]].elist[j].snp) 
			{
				W = (float)edge_weight(hap,slist[i],k,snpfrag[slist[i]].elist[j].p,Flist,snpfrag[slist[i]].elist[j].frag);
				if (wf ==0) W /= Flist[snpfrag[slist[i]].elist[j].frag].calls -1; //(fraglength(Flist,snpfrag[slist[i]].elist[j].frag)-1); 
				snpfrag[slist[i]].telist[snpfrag[slist[i]].tedges-1].w += W; 
			}
		}
	}
	/* CODE TO set up the read-haplotype consistency graph */

	/* CODE TO find 'K' biggest edges in MEC graph  */
	int K = 5; int smallest = 0; float smallw=1000;
	if (totaledges/2 < K) K = totaledges/2;
	EDGE* edgelist = (EDGE*)malloc(sizeof(EDGE)*K); j=0; i= 0; k= 0;
	for (i=0;i<N;i++)
	{
		for (j=0;j<snpfrag[slist[i]].tedges;j++)
		{
			if (k < K)
			{
				edgelist[k].s = slist[i]; edgelist[k].t = snpfrag[slist[i]].telist[j].snp;
				edgelist[k].w = snpfrag[slist[i]].telist[j].w;
				if (edgelist[k].w < smallw) { smallest = k; smallw = edgelist[k].w; }
				k++;
			}
			else
			{
				if (snpfrag[slist[i]].telist[j].w > smallw)
				{
					edgelist[smallest].s = slist[i]; edgelist[smallest].t = snpfrag[slist[i]].telist[j].snp;
					edgelist[smallest].w = snpfrag[slist[i]].telist[j].w;
					smallw = 1000;
					for (l=0;l<K;l++)
					{
						if (edgelist[l].w < smallw) { smallest = l; smallw = edgelist[l].w; }
					}
				}
			}
		}
	}
	/* CODE TO find 'K' biggest edges in MEC graph  */

	// edge contraction algorithm: merge vertices until only two nodes left or total edge weight of graph is negative  
	int startnode =(int)(drand48()*N); if (startnode ==N) startnode--;  
	int secondnode= -1 ;  // root of 2nd cluster initially not there
	// chose a positive edge to initialize the two clusters and run this algorithm $O(m)$ times for each block 
	// a negative weight cut should have at least one negative edge or if there is no negative weight edge, the edge with lowest weight 

	for (i=0;i<N;i++) { snpfrag[slist[i]].revmap = i; }
	int V = N; 
	float curr_cut=0,best_cut=10000; float oldscore;
	int snp_add;
	int moved =1,c1=0,c2=0; char* bestmincut;
//	int size_small,best_small=0,secondlast=0,last=0;
	int iter=0,maxiter=N/10; 
	if (N/10 < 1) maxiter =1; 
	if (maxiter >= MAXCUT_ITER && MAXCUT_ITER >= 1) maxiter = MAXCUT_ITER;  // added march 13 2013

	int fixheap =0;  PHEAP pheap; pinitheap(&pheap,N); // heap for maxcut 

	/*****************************Maintain two clusters and add each vertex to one of these two ******************/
	bestmincut = (char*)malloc(N); for (i=0;i<N;i++) bestmincut[i] = '0';
	//for (iter=0;iter<totaledges*(int)(log2(totaledges));iter++)
	for (iter=0;iter<maxiter+K;iter++)
	{
		if (iter < K) 
		{ 
			startnode = edgelist[iter].s; secondnode = edgelist[iter].t; 
		} // first K edges are not random but rather the top 'K' edges by weight
		else
		{
			i = (int)(drand48()*totaledges-0.0001); j=0;
			while (i >= snpfrag[slist[j]].tedges) { i -= snpfrag[slist[j]].tedges; j++;} 
			startnode = slist[j]; secondnode = snpfrag[slist[j]].telist[i].snp; 
			if (snpfrag[slist[j]].telist[i].w >=1) continue; 
		}
		for (i=0;i<N;i++) snpfrag[slist[i]].parent = slist[i];

		// new code added for heap based calculation
		for (i=0;i<N;i++) snpfrag[slist[i]].score  = 0;
		pheap.length = N-2; j=0; // heap only has N-2 elements (startnode and secondnode are not there) 
		for (i=0;i<N;i++) 
		{
			if (slist[i] != startnode && slist[i] != secondnode) { pheap.elements[j] = i; snpfrag[slist[i]].heaploc = j; j++; } 
		}
//		for (i=0;i<N;i++) fprintf(stdout,"heaploc %d %d %d-%d\n",slist[i],snpfrag[slist[i]].heaploc,startnode,secondnode);

		// for all vertices linked to startnode -> update score of the neighbors, can be done using list of fragments for snpfrag[startnode]
		for (j=0;j<snpfrag[startnode].tedges;j++) 
		{
			k = snpfrag[snpfrag[startnode].telist[j].snp].revmap; 
			snpfrag[slist[k]].score  += snpfrag[startnode].telist[j].w;
			//if (slist[k] == startnode) fprintf(stdout,"self edge %d %d\n",slist[k],startnode);
		}
		for (j=0;j<snpfrag[secondnode].tedges;j++)
		{
			k = snpfrag[snpfrag[secondnode].telist[j].snp].revmap; 
			snpfrag[slist[k]].score  -= snpfrag[secondnode].telist[j].w;
			//if (slist[k] == secondnode) fprintf(stdout,"self edge %d %d\n",slist[k],secondnode);
		}

		pbuildmaxheap(&pheap,snpfrag,slist); 
		V = N; while (V > 2) // more than two clusters, this loop is O(N^2) and makes the program slow
		{
			// need to replace the code below by a pointer-based heap that stores 'k' at each node which -> slist[k] -> snpfrag[slist[k]].score 
			snp_add = pheap.elements[0]; premovemax(&pheap,snpfrag,slist);  fixheap = 0;
			//if (N < 30) fprintf(stdout,"standard best score %f snp %d %d V %d\n",snpfrag[slist[snp_add]].score,snp_add,slist[snp_add],V);
			if (snpfrag[slist[snp_add]].score > 0) snpfrag[slist[snp_add]].parent = startnode;  
			else if (snpfrag[slist[snp_add]].score < 0 ) snpfrag[slist[snp_add]].parent = secondnode;  
			else if (drand48() < 0.5) snpfrag[slist[snp_add]].parent = startnode; 
			else snpfrag[slist[snp_add]].parent = secondnode;  
			V--;
			for (j=0;j<snpfrag[slist[snp_add]].tedges;j++)
			{
				k = snpfrag[snpfrag[slist[snp_add]].telist[j].snp].revmap;
				if (slist[k] == startnode || slist[k] == secondnode)  continue;
				// serious bug identified here from wheat data, we should discard edges from snp_add to startnode/secondnode since these nodes are nott in the heap 
				if (snpfrag[slist[snp_add]].parent == startnode) 
				{
					oldscore = snpfrag[slist[k]].score; 
					snpfrag[slist[k]].score += snpfrag[slist[snp_add]].telist[j].w; 
				} 
				else if (snpfrag[slist[snp_add]].parent == secondnode) 
				{
					oldscore = snpfrag[slist[k]].score; 
					snpfrag[slist[k]].score -= snpfrag[slist[snp_add]].telist[j].w;
				}
				if (fabsf(oldscore) > fabsf(snpfrag[slist[k]].score)) // score decreased 
				{
					//fixheap = 1;
					//fprintf(stdout,"maxheapify iter %d %d slist %d %d %d\n",iter,V,slist[k],snpfrag[slist[k]].heaploc,slist[snp_add]);
					pmaxHeapify(&pheap,snpfrag[slist[k]].heaploc,snpfrag,slist);
				}
				else // absolute score of node increased 
				{
					//fprintf(stdout,"bubble iter %d %d slist %d %d %d\n",iter,V,slist[k],snpfrag[slist[k]].heaploc,slist[snp_add]);
					pbubbleUp(&pheap,snpfrag[slist[k]].heaploc,snpfrag,slist);
					//fixheap = 1;
				}
			}
			if (fixheap ==1) pbuildmaxheap(&pheap,snpfrag,slist);
		}


		// compute score of the cut computed above 
		for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == startnode) snpfrag[slist[i]].parent =0; else snpfrag[slist[i]].parent = 1;  } 
		c1=0;c2=0; for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == 0) c1++; else c2++;  } 
		curr_cut =0; 
		for(i=0;i<N;i++) 
		{ 
			for (j=0;j<snpfrag[slist[i]].tedges;j++)
			{
				if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent != snpfrag[slist[i]].parent) 		curr_cut += snpfrag[slist[i]].telist[j].w/2; 
			}
		}
		moved =1; while (moved > 0) // any improvement in score of cut  
		{
			moved =0;
			for (i=0;i<N;i++)
			{
				snpfrag[slist[i]].score  = 0;
				for (j=0;j<snpfrag[slist[i]].tedges;j++)
				{
					if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == 0) snpfrag[slist[i]].score += snpfrag[slist[i]].telist[j].w;
					if (snpfrag[snpfrag[slist[i]].telist[j].snp].parent == 1) snpfrag[slist[i]].score -= snpfrag[slist[i]].telist[j].w;
				}
				if (snpfrag[slist[i]].parent ==0 && snpfrag[slist[i]].score < 0 && c1 > 1)
				{
					snpfrag[slist[i]].parent =1; curr_cut += snpfrag[slist[i]].score; moved++; c1--; c2++;
				}
				if (snpfrag[slist[i]].parent ==1 && snpfrag[slist[i]].score > 0 && c2 > 1)
				{
					snpfrag[slist[i]].parent =0; curr_cut -= snpfrag[slist[i]].score; moved++; c2--; c1++;
				}
			}
		} if (c1 ==0 || c2==0) { fprintf(stdout," cut size is 0 red \n"); exit(0); } 
		if (curr_cut < best_cut)
		{
			best_cut = curr_cut; 
			for (i=0;i<N;i++) { if (snpfrag[slist[i]].parent == 1) bestmincut[i] ='1'; else bestmincut[i] = '0'; } 
		}
		//																								if (best_cut < -5) iter = maxiter; 
	}
	for(i=0;i<N;i++) { if (bestmincut[i] == '1') slist[i]= -1*slist[i]-1; }
	free(bestmincut); free(edgelist); free(pheap.elements); 
	return best_cut;
}

