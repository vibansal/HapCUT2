// CODE STARTED SEPT 10 2007 4pm //  
// april 8 2008 this code used for producing results in ECCB 2008 paper
// author: VIKAS BANSAL (vbansal@scripps.edu) last modified December 23, 2010 

// time ./HAPCUT --fragments /home/vbansal-scripps/Haplotype-assembly/Dec8hairs/NA18507.allfour.matrix.SORTED --variants NA18507.chr6.hetvariants.inputforhapcut --output NA18507.phased --maxiter 60 > na18507.out 

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "common.h"
#include "fragmatrix.h"
#include "pointerheap.h"
#include "annealing.h"
#include "readinputfiles.h"
#include "removefalsehets.c"
//#define MAXBUF 10000

int RANDOM_START=1;
int MINCUTALGO =1;
int QVoffset= 33;
int VCFformat =0;
int MINQ = 10;// additional base quality filter in hapcut added april 18 2012

int maxcut_haplotyping(char* fragmentfile,char* variantfile,int snps,char* outputfile,int maxiter_hapcut)
{
	// IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+
	fprintf(stderr,"calling MAXCUT based haplotype assembly algorithm\n");
	int fragments=0,iter=0,components=0; int i=0,j=0,k=0,t=0,component;	int* slist; 
	float bestscore_mec = 0,calls=0, miscalls=0,ll = 0;
	char buffer[MAXBUF];
	int flag =0;

	/****************************** READ FRAGMENT MATRIX*************************************************/
	struct fragment* Flist; FILE* ff = fopen(fragmentfile,"r"); 
	if (ff == NULL) { fprintf(stderr,"couldn't open fragment file %s\n",fragmentfile); exit(0);}
	fragments =0; while ( fgets(buffer,MAXBUF,ff) != NULL) fragments++; fclose(ff);
	Flist = (struct fragment*)malloc(sizeof(struct fragment)*fragments); 
	flag = read_fragment_matrix(fragmentfile,Flist,fragments);
	if (flag < 0) { fprintf(stderr,"unable to read fragment matrix file %s \n",fragmentfile); return -1; } 

	if (VCFformat ==0) snps = count_variants(variantfile);
	else snps = count_variants_vcf(variantfile);
	if (snps < 0) { fprintf(stderr,"unable to read variant file %s \n",variantfile); return -1; } 
	fprintf(stderr,"fragments %d variants %d\n",fragments,snps); 
	/****************************** READ FRAGMENT MATRIX*************************************************/

	struct SNPfrags* snpfrag = (struct SNPfrags*)malloc(sizeof(struct SNPfrags)*snps); 
	update_snpfrags(Flist,fragments,snpfrag,snps,&components);
	//for (i=0;i<snps;i++) fprintf(stdout,"flist %d %d %d \n",i,snpfrag[i].frags,snpfrag[i].flist[0]);
	components = determine_connected_components(Flist,fragments,snpfrag,snps);
	struct BLOCK* clist = (struct BLOCK*)malloc(sizeof(struct BLOCK)*components);

	/*****************************************************************************************************/

	char* HAP1 = (char*)malloc(snps+1); char* besthap_mec = (char*)malloc(snps+1);
	char* HAP2 = (char*)malloc(snps+1);
	struct tm  *ts1;   char       buf[80];	time_t     now;		
	slist = (int*)malloc(sizeof(int)*snps); char fn[1000];  

	if (RANDOM_START ==1)
        {
                fprintf(stdout,"starting from a completely random solution SOLUTION \n");
                for (i=0;i<snps;i++)
                {
                        if (snpfrag[i].frags ==0) { HAP1[i] = '-'; HAP2[i] = '-'; }
                        else
                        {
                                if (drand48() < 0.5) { HAP1[i] = '0'; HAP2[i] = '1'; }  else  {HAP1[i] = '1'; HAP2[i] = '0'; }
                        }
                }
        }

	if (VCFformat ==0) read_variantfile(variantfile,snpfrag,snps); else read_vcffile(variantfile,snpfrag,snps);
	generate_clist_structure(Flist,fragments,snpfrag,snps,components,clist); //return 1; 

	/*****************************************************************************************************/
	int maxiter = 10;
	annealing_haplotyping(Flist,fragments,snpfrag,snps,maxiter,HAP1,HAP2,clist,components,slist); 
	for (i=0;i<components;i++) find_bestvariant_segment(Flist,fragments,snpfrag,clist,i,HAP1,HAP2); 
	sprintf(fn,"%s",outputfile);   // newfile for every update to score....
        fprintf(stdout,"OUTPUTTING HAPLOTYPE ASSEMBLY TO FILE %s\n",fn);        fprintf(stderr,"OUTPUTTING HAPLOTYPE ASSEMBLY TO FILE %s\n",fn);
        print_hapfile(clist,components,HAP1,Flist,fragments,snpfrag,variantfile,miscalls,fn);
	return 1;
	//annealing_haplotyping_full(Flist,fragments,snpfrag,snps,maxiter,HAP1,HAP2,0); return 1;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	// input arguments are initial fragment file, variant file with variant information and alleles for each variant 
	// number of iterations total, when to output the solution, file to output solution .....
	time_t ts; time(&ts); srand48((long int)ts);
	if (MINCUTALGO ==2) RANDOM_START=0;
	int i=0; int flag =0;
	char fragfile[10000]; char varfile[10000]; char VCFfile[10000]; char hapfile[10000]; int maxiter =100; 
	strcpy(fragfile,"None"); strcpy(varfile,"None");strcpy(hapfile,"None"); 
	for (i=1;i<argc;i+=2)
	{
		if (argc < 6) break;
		if (strcmp(argv[i],"--fragments") ==0)     {   strcpy(fragfile,argv[i+1]); flag++; }
		else if (strcmp(argv[i],"--variants") ==0)    {     strcpy(varfile,argv[i+1]); flag++; }
		else if (strcmp(argv[i],"--VCF") ==0)    {   strcpy(VCFfile,argv[i+1]); VCFformat = 1; flag++; }
		else if (strcmp(argv[i],"--output") ==0)      {   strcpy(hapfile,argv[i+1]); flag++; }
		else if (strcmp(argv[i],"--maxiter") ==0)        maxiter = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--QVoffset") ==0)        QVoffset = atoi(argv[i+1]);
		else if (strcmp(argv[i],"--qvoffset") ==0)        QVoffset = atoi(argv[i+1]);
	}
	if (flag !=3) 
	{ 
		print_hapcut_options(); return 0;
	}
	else
	{
		fprintf(stderr,"\n\nfragment file: %s\nvariantfile (VCF format):%s\nhaplotypes will be output to file: %s\niterations of maxcut algorithm: %d\nQVoffset: %d\n\n",fragfile,VCFfile,hapfile,maxiter,QVoffset); 	
		maxcut_haplotyping(fragfile,VCFfile,0,hapfile,maxiter); 
	}
	return 1;
}




