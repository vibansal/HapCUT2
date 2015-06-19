
// code written april 6 2012 //
// for each connected component, calculate the likelihood increase by removing a variant and if good, do it

void find_bestvariant_segment(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,struct BLOCK* clist, int k,char* HAP1,char* HAP2);
int removevariants(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int maxiter,char* HAP1,char* HAP2,struct BLOCK* clist,int components);
float evaluate_variant(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, int v,char* HAP1);

// phase each component individually
void find_bestvariant_segment(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,struct BLOCK* clist, int k,char* HAP1,char* HAP2)
{
	int i=0;
	//	fprintf(stdout,"\ncomponent offset: %d length %d delta scores MEC %f\n",clist[k].offset+1,clist[k].length,clist[k].bestMEC);
	for (i=0;i<clist[k].phased;i++)
	{
		//fprintf(stderr,"flipping variant %d ll %f \n",i,ll);
		//t = clist[k].slist[i]; 
		//delta = evaluate_variant(Flist,fragments,snpfrag,t,HAP1);
		//fprintf(stdout,"%s_%d_%s_%s_%s ",snpfrag[t].chromosome,snpfrag[t].position,snpfrag[t].allele0,snpfrag[t].allele1,snpfrag[t].genotypes);
		//fprintf(stdout," R0 %d %d SNP %d \n",snpfrag[t].R0,snpfrag[t].R1,t);
		//fprintf(stdout,"DELTA %f %f %f %f %f\n",delta,likelihoods[0],likelihoods[1],likelihoods[2],likelihoods[3]);
	}
}

// calculate how much is MEC score reduced by removing a variant 
float calculate_reduction_mecscore(struct fragment* Flist,int f, char* h,int variant)
{
	int j=0,k=0; float good=0,bad=0; float prob =0,prob1=0,origscore=0,newscore =0; float good1=0,bad1=0;
	for (j=0;j<Flist[f].blocks;j++)
	{
		for (k=0;k<Flist[f].list[j].len;k++)
		{
			if (h[Flist[f].list[j].offset+k] == '-') continue;// { fprintf(stdout,"fragment error"); continue;}
			if ((int)Flist[f].list[j].qv[k] -QVoffset < MINQ) continue; 
			prob = QVoffset-(int)Flist[f].list[j].qv[k]; prob /= 10; prob1 = 1.0; prob1 -= pow(10,prob);
			if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good +=prob1; else bad +=prob1;
			if (Flist[f].list[j].offset+k != variant) 
			{
				if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good1 +=prob1; else bad1 +=prob1;
			}
		}
	}
	if (good < bad) origscore = good; else origscore = bad;
	if (good1 < bad1) newscore = good1; else newscore = bad1;
	return origscore - newscore; 
}


// variant can be homozygous reference, homozygous alternate, or het (0/0, 1/1, 0|1, 1|0) or due to copy number polymorphism (0|0|1, 0|1|1) 
// calculate delta in log-likelihood on removal of a variant from connected component
// also calculate the impact on the likelihood if the variant is removed altogether, don't care.... 
float evaluate_variant(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, int v,char* HAP1)
{
	int i=0,j=0,t=0,t1=0,t2=0;
	float Forig1 = 0,Forig2 = 0, Fnew1 = 0, Fnew2=0,delta =0; 
	float Fnew3 =0,Fnew4=0,Fnew5=0,Fnew6=0;
	float prob = 0.0,prob1=0.0;
	float L01 = 0, L00 =0,L11=0,Lnovar=0;
	int R0=0,R1=0;
	snpfrag[v].rMEC = 0; // reduction in MEC score if we remove this variant altogether
	snpfrag[v].G00 = 0; snpfrag[v].G01 = 0; snpfrag[v].G11 = 0; snpfrag[v].A0 = 0; snpfrag[v].A1=0;
	//float LOG2 = log10(0.5); LOG2 =0;

	// also calculate reduction in MEC score if we remove this variant or make it homozygous 

	// should we be using all fragments that cover variant 'v' -> snpfrag[v].flist and frags have that info 
	// this is correct, and fixes the bug in genotype probabilities, april 17 2012 
	for (i=0;i<snpfrag[v].frags;i++)
	{
		j = snpfrag[v].flist[i]; 
		if (Flist[j].blocks ==1 && Flist[j].list[0].len ==1)  // singleton read fragment
		{
			t = Flist[j].list[0].offset; if (HAP1[t] == '-') continue;
			if ((int)Flist[j].list[0].qv[0] -QVoffset < MINQ) continue;
			prob = QVoffset-(int)Flist[j].list[0].qv[0]; prob /= 10; prob1 = log10(1.0 - pow(10,prob));
			if (Flist[j].list[0].hap[0] == '0')
			{
				snpfrag[v].A0++; snpfrag[v].G00 += prob1; snpfrag[v].G11 += prob; snpfrag[v].G01 += log10(0.5);
			}
			else if (Flist[j].list[0].hap[0] == '1')
			{
				snpfrag[v].A1++; snpfrag[v].G00 += prob; snpfrag[v].G11 += prob1; snpfrag[v].G01 += log10(0.5);
			}
			continue;
		}

		snpfrag[v].rMEC += calculate_reduction_mecscore(Flist,j,HAP1,v);

		Forig1=0,Forig2 =0, Fnew1 = 0, Fnew2=0,Fnew3=0,Fnew4=0, Fnew5=0,Fnew6=0;
		for (t1=0;t1<Flist[j].blocks;t1++)
		{
			for (t2=0;t2<Flist[j].list[t1].len;t2++)
			{
				t = Flist[j].list[t1].offset + t2;
				if (HAP1[t] == '-') continue;
				if ((int)Flist[j].list[t1].qv[t2] -QVoffset < MINQ) continue; 
				prob = QVoffset-(int)Flist[j].list[t1].qv[t2]; prob /= 10; prob1 = log10(1.0 - pow(10,prob));

				if (t ==v) // t =v 
				{
					if (Flist[j].list[t1].hap[t2] == HAP1[t]) { Forig1 += prob1;  Forig2 += prob;  }
					else {  Forig1 += prob; Forig2 += prob1; }
					if (Flist[j].list[t1].hap[t2] == '0') { Fnew1 += prob1; Fnew2 += prob1; Fnew3 += prob; Fnew4 += prob; R0++;  } 
					else { Fnew1 += prob; Fnew2 += prob; Fnew3 += prob1; Fnew4 += prob1; R1++;} 
				}
				else 
				{
					if (Flist[j].list[t1].hap[t2] == HAP1[t])
					{
						Forig1 += prob1;  Forig2 += prob; Fnew1 += prob1; Fnew2 += prob; Fnew3 += prob1; Fnew4 += prob;
						Fnew5 += prob1; Fnew6 += prob; 
					}
					else
					{
						Forig1 += prob; Forig2 += prob1; Fnew1 += prob; Fnew2 += prob1; Fnew3 += prob; Fnew4 += prob1; 
						Fnew5 += prob; Fnew6 += prob1; 
					}
				}
			}
		}
		if (Fnew1 > Fnew2) L00 += Fnew1 + log10(1.0+pow(10,Fnew2-Fnew1)); else L00 += Fnew2 + log10(1.0+pow(10,Fnew1-Fnew2)); 
		if (Fnew3 > Fnew4) L11 += Fnew3 + log10(1.0+pow(10,Fnew4-Fnew3)); else L11 += Fnew4 + log10(1.0+pow(10,Fnew3-Fnew4)); 
		if (Fnew5 > Fnew6) Lnovar += Fnew5 + log10(1.0+pow(10,Fnew6-Fnew5)); else Lnovar += Fnew6 + log10(1.0+pow(10,Fnew5-Fnew6)); 
		if (Forig1 > Forig2) L01 += Forig1 + log10(1.0+pow(10,Forig2-Forig1)); else L01 += Forig2 + log10(1.0+pow(10,Forig1-Forig2)); 
	}
	snpfrag[v].L00 = L00; snpfrag[v].L01 = L01; snpfrag[v].L11 = L11; snpfrag[v].Lnovar = Lnovar;
	snpfrag[v].R0 = R0; snpfrag[v].R1 = R1;
	delta = L00-L01; 
	if (L11-L01 > delta) delta = L11-L01;
	return delta;
}


// not used for now april 9 2012
int removevariants(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag,int snps,int maxiter,char* HAP1,char* HAP2,struct BLOCK* clist,int components)
{
	int iter=0,k=0;
	float calls=0, miscalls=0,ll=0, bestll =0;
	//int switches =0; int prev = 0;
	/*****************************************************************************************************/
	//for (i=0;i<snps;i++) fprintf(stdout,"first fragment for SNP %d %d \n",i,snpfrag[i].ff);
	mecscore(Flist,fragments,HAP1,&ll,&calls,&miscalls); bestll = ll; 

	for (iter =0; iter < 1;iter++)
	{
		fprintf(stdout,"MEC score %f %f  LL %f bestLL %f \n",miscalls,calls,ll,bestll);
		fprintf(stderr,"MEC score %f %f LL %f bestLL %f\n",miscalls,calls,ll,bestll);

		for (k=0;k<components;k++) find_bestvariant_segment(Flist,fragments,snpfrag,clist,k,HAP1,HAP2); 				
		//mecscore(Flist,fragments,HAP1,&ll,&calls,&miscalls); if (bestll > ll) bestll = ll;
	}
	return 1;
}
