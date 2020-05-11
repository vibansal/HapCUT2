/* functions for comparing an aligned sequence read to the set of variants to identify alleles and haplotype-informative reads */

//#include "../realign_pairHMM.c"
#include "nw.c"
#include <assert.h>
#include <stdlib.h>

int MINLEN= 8;
int COMPLEXITY_K = 5; // anchor sequences must have unique kmers of this length
int SHORT_HAP_CUTOFF = 20; // separation between variants to be considered part of the same local haplotype for realignment
extern int VERBOSE;
extern int REALIGN_ALL;
extern int PARSEINDELS;
float TINYLOG = -10000;
//int MIN_QUAL = 10;
int MAX_SNPs_SHORT_HAP = 10; // max number of SNVs in a short haplotype
// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 - pow(10, (b) - (a)))) : ((b) + log10(1.0 - pow(10.0, (a) - (b)))))

int compare_strings(const void* a, const void* b)
{
	const char *ia = *(const char **)a;        //const char *ia = (const char *)a;
	const char *ib = *(const char **)b;  //const char *ib = (const char *)b

	return strcmp(ia, ib);
}

// https://stackoverflow.com/questions/2509679/how-to-generate-a-random-number-from-within-a-range
unsigned int rand_interval(unsigned int min, unsigned int max){
	int r;
	const unsigned int range = 1 + max - min;
	const unsigned int buckets = RAND_MAX / range;
	const unsigned int limit = buckets * range;

	/* Create equal size buckets all in a row, then fire randomly towards
	 * the buckets until you land in one of them. All buckets are equally
	 * likely. If you land off the end of the line of buckets, try again. */
	do{
		r = rand();
	} while (r >= limit);

	return min + (r / buckets);
}

// indel ambiguity needs to be accounted when constructing haplotypes...

// generate extended cigar list, memory should already be allocated..
int parse_cigar(struct alignedread* read,REFLIST* reflist,int* fcigarlist)
{
	int current = reflist->current; if (current < 0) return -1;
	int f=0;
	int i=0,t=0, l1=0,l2=0; int l=0; int m=0; int op;
	for (i=0;i<read->cigs;i++)
	{
		op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4;
		if (op == BAM_CMATCH)
		{
			m=0;
			for (t=0;t<l;t++)
			{
				if (read->sequence[l1+t]  != reflist->sequences[current][read->position+l2+t-1] && read->sequence[l1+t] != reflist->sequences[current][read->position+l2+t-1]-32 && read->sequence[l1+t]  !='N')
				{
					read->mismatches++;
					if (m > 0) fcigarlist[f++] = (m<<4)+7; //  m=
					fcigarlist[f++] = 16 + 8; // mismatch 1X
					m=0;
				}
				else m++;
			}
			if (m > 0) fcigarlist[f++] = (m<<4)+7;
		}
		else fcigarlist[f++] = read->cigarlist[i];  // copy as it is for some cigar operations
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)  { l1 +=l; l2 +=l; }
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP) l2 +=l;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)  l1 += l;
	}
	return f;
}

// returns 1 if kmers in seq are unique, else 0
int test_complexity(char* seq, int k){
	//fprintf(stderr,"seq: %s",seq);
	int i = 0, j = 0, q = 0, match = 1;
	int l = strlen(seq);

	for (i = 0; i < l - k + 1; i++){
		for (j = 0; j < l - k + 1; j++){
			if (i == j) continue;
			match = 1;
			for (q = 0; q < k; q++){
				if (seq[i+q] != seq[j+q]){
					match = 0;
					break;
				}
			}

			if (match){
				if (VERBOSE) fprintf(stderr," fail %s\n",seq);
				return 0;
			}
		}
	}
	return 1;
}

// measure complexity of a sequence using the number of distinct kmers, avoid repetitive sequence
// AAAAAAAAA is bad, ACACACAGAGAGAG is also bad... 
/*
int complexity(char* seq,int n)
{
	int k=2;
	if (n >= 10) k = 3; 
	//int max_count = n-k+1;  // min of (n-k+1,4^k)
	int* kmers = calloc(sizeof(short),n-k+1);
	int i=0,j=0;
	int kmer2int=0;
	for (i=0;i<n-k;i++)
	{
		kmer2int=0;
		for (j=i;j<i+k;j++) kmer2int *=4;
		kmers[i]= kmer2int; 
	}
	// sort the kmers list and count number of unique k-mers
	free(kmers);
}
*/

// realign sequence read to haplotype sequences defined by 'k' variants (typically k=1), local realignment for computing quality scores
// ALLELOTYPING 
int realign_HAPs(struct alignedread* read, REFLIST* reflist, int positions[], VARIANT* varlist, int* snplst, int n_variants, FRAGMENT* fragment)
{
	if ( positions[3]-positions[1] < 10 ||  positions[3]-positions[1] > 200 || positions[2]-positions[0] < 10 ||  positions[2]-positions[0] > 200) 
	{
		if  (VERBOSE) fprintf(stderr,"range check failed \n");
		return -1; // check on length in range (15,200)
	}
	if (positions[0] < 0 || positions[2]-1 >= strlen(read->sequence)){
		fprintf(stderr,"BOUNDS CHECK FAIL: read %s, read length %d, start %d, end %d\n", read->readid, read->readlength, positions[0], positions[2]-1);
		exit(1);
	}

	int i=0,j=0,k=0;
	char* subread = malloc(positions[2]-positions[0]+1);

	for (j=positions[0];j<positions[2];j++)   subread[j-positions[0]] = read->sequence[j];
	subread[j-positions[0]]='\0';

	char* refhap = malloc(positions[3]-positions[1]+1);
	// don't forget +1 to strlen for end character
	for (j=positions[1];j<positions[3];j++)  refhap[j-positions[1]] = reflist->sequences[reflist->current][j-1];
	refhap[j-positions[1]]  ='\0';

	char* althap = malloc(positions[3]-positions[1]+1 + 4096); // edit on 11/26/18
	int h=0, s=0, ss=0, max_hap=0, ref_len=0, alt_len=0, total_ref_len=0, total_alt_len=0, n_max_haps = 0, rand_ix = 0;
	double total_score = TINYLOG, max_score = -1000000000;
	int align_qual = 0;
	int* max_haps = (int*) calloc(pow(2,n_variants),sizeof(int));

	double* ref_score_single = malloc(MAX_SNPs_SHORT_HAP * sizeof(double));
	double* alt_score_single = malloc(MAX_SNPs_SHORT_HAP * sizeof(double));
	for (s = 0; s < MAX_SNPs_SHORT_HAP; s++){
		// initialize to very small log value
		ref_score_single[s] = TINYLOG;
		alt_score_single[s] = TINYLOG;
	}

	double altscore=0;
	int r =0;

	//if (VERBOSE) fprintf(stderr,"%s subread\n",subread);
	if (VERBOSE) fprintf(stderr,"read %d:%d:%d ",positions[0],positions[2],positions[2]-positions[0]);
	if (VERBOSE) fprintf(stderr,"on ref %d:%d:%d\n",positions[1],positions[3],positions[3]-positions[1]);
	for (h = 0; h < pow(2,n_variants); h++){

		// we represent a haplotype with integer h
		// in particular the n_variants least significant bits
		// if the i-th bit from the right is 1, computed as h & (pow(2,i)), then the haplotype contains the i-th variant

		total_ref_len = 0;
		total_alt_len = 0;
		for (s = 0; s < n_variants; s++)
		{
			ss = snplst[s];
			total_ref_len += strlen(varlist[ss].allele1);
			if (h & (int)(pow(2,s))) total_alt_len += strlen(varlist[ss].allele2);
			else 	total_alt_len += strlen(varlist[ss].allele1);
		}

		if (VERBOSE) fprintf(stderr,"new alt hap len %d %d\n",positions[3]-positions[1]+1+total_alt_len-total_ref_len,total_alt_len);
		//althap = malloc(positions[3]-positions[1]+1+total_alt_len-total_ref_len);

		j = positions[1]; k = 0;
		for (s = 0; s < n_variants; s++)
		{
			ss = snplst[s];
			ref_len = strlen(varlist[ss].allele1);
			alt_len = strlen(varlist[ss].allele2);

			if (varlist[ss].type != 0) r = varlist[ss].position-1;
			else  r = varlist[ss].position;  // added condition so that the position is correct for indels, 08/20/19

			for (j = j; j<r; j++) 
			{
				althap[k++] = reflist->sequences[reflist->current][j-1];
				//if (VERBOSE) fprintf(stderr,"%d %c ",j,reflist->sequences[reflist->current][j-1]);
			}
			if (h & (int)(pow(2,s)) )
			{
				for (i=0;i<alt_len;i++) althap[k++] = varlist[ss].allele2[i];
			}
			else
			{
				for (i=0;i<ref_len;i++) althap[k++] = varlist[ss].allele1[i];
			}
			j += ref_len;
			//if (VERBOSE) fprintf(stderr,"%d %s ",j,varlist[ss].allele1);
		}

		for (j = j; j < positions[3]; j++) 
		{
			althap[k++] = reflist->sequences[reflist->current][j-1];
			//if (VERBOSE) fprintf(stderr,"%d %c ",j,reflist->sequences[reflist->current][j-1]);
		}
		althap[k] = '\0';

		if (SUM_ALL_ALIGN ==0) altscore = nw(althap,subread,0);
		else if (SUM_ALL_ALIGN ==1) altscore = sum_all_alignments_fast(althap,subread,AP,20);
		else if (SUM_ALL_ALIGN ==2) altscore = sum_all_alignments_logspace(althap,subread,AP,20);  

		if (VERBOSE) fprintf(stderr,"h %d score: %f \n%s haplotype\n%s subread\n",h,altscore,althap,subread);

		// for an index s in the short haplotype,
		// maintain the log sum of scores that have a variant at s
		// and the log sum of scores that have reference at s
		for (s = 0; s < n_variants; s++){

			ss = snplst[s];
			if (h & (int)(pow(2,s))){
				// current haplotype has variant at s-th variant
				alt_score_single[s] = addlogs(alt_score_single[s], altscore);
			}else{
				// current haplotype has reference at s-th variant
				ref_score_single[s] = addlogs(ref_score_single[s], altscore);
			}
		}

		// if this is the current max haplotype alignment score then store it
		if(altscore > max_score){
			max_score = altscore;
			n_max_haps = 0;
			memset(max_haps, 0, pow(2,n_variants)*sizeof(int));
		}
		// in case of a tie we store every max haplotype so we can randomly select one
		if(altscore == max_score){
			max_haps[n_max_haps] = h;
			n_max_haps++;
		}

		// add the alternate score to the total
		total_score = addlogs(total_score, altscore);
		//free(althap);
	}

	// if the max score is sufficiently high then write alleles to haplotype fragments and calculate quality values from scores

	// TODO
	// need to assign each SNVs' score to (sum of haplotypes with variant)/(sum of all haplotypes)
	rand_ix = rand_interval(0,n_max_haps-1);
	max_hap = max_haps[rand_ix];

	for (s = 0; s < n_variants; s++){

		ss = snplst[s];

		if (varlist[ss].type != 0 && !PARSEINDELS){
			continue;
		}
		fragment->alist[fragment->variants].varid = ss;

		if (max_hap & (int)(pow(2,s))){
			align_qual = (int) (-10.0 * (ref_score_single[s] - total_score));
			if(VERBOSE) fprintf(stderr,"varid %d alt-alele %f qv %d %30s\n",ss,ref_score_single[s],align_qual,read->readid);

			if (align_qual < MINQ) continue;

			fragment->alist[fragment->variants].allele = '1';
		}else{
			align_qual = (int) (-10.0 * (alt_score_single[s] - total_score));
			if(VERBOSE) fprintf(stderr,"varid %d ref-alele %f qv %d %30s\n",ss,ref_score_single[s],align_qual,read->readid);

			if (align_qual < MINQ) continue;

			fragment->alist[fragment->variants].allele = '0';
		}

		if (align_qual + QVoffset > 126){
			align_qual = 126 - QVoffset;
		}

		fragment->alist[fragment->variants].qv = (char) (align_qual + QVoffset);
		fragment->variants++;
		varlist[ss].depth++;
		if ((read->flag & 16) == 16) varlist[ss].A1 += 1 << 16;
		else varlist[ss].A1 += 1;
	}
	if (VERBOSE) fprintf(stderr,"**********************************************\n");

	free(subread);free(refhap); free(max_haps);free(ref_score_single); free(alt_score_single);
	free(althap);

	return 0;
}

// offset is distance between variant position and the start of the cigar operation 'f1' on the reference genome (left_anchor_ref)
int find_left_anchor(int* fcigarlist,int fcigs,int f1,REFLIST* reflist,int offset,int* left_anchor_read, int* left_anchor_ref)
{
	char* anchor_seq = malloc(MINLEN+1);
	int op,ol;
	int anchor_start = 0, anchor_end = 0, j = 0;
	int flag = 1;
	//int kmer_test=0;

	// if current cigar is CEQUAL and the variant position is at least MINLEN bases after the cigar start, we are done
	op = fcigarlist[f1]&0xf; ol = fcigarlist[f1]>>4;
	if ((offset >= MINLEN) && op == BAM_CEQUAL && ol >= offset) 
	{
		if (VERBOSE) fprintf(stderr,"done in first pass cigar %d \n",offset);
		if (offset > MINLEN) // too long
		{
			*left_anchor_read += offset-MINLEN; 
			*left_anchor_ref +=  offset-MINLEN;
		}
		flag =0; // don't need to change anchor positions could be too large ??
	}
	else f1--; // go to previous cigar, added aug 2019, seems to be correct since left_anchor_ref corresponds to f1
	// surprisingly, it doesn't affect overall results at all..

	// keep going left until we find a cigar CMATCH of length at least MINLEN 
	while (f1 >= 0 && flag ==1)
	{
		op = fcigarlist[f1]&0xf; ol = fcigarlist[f1]>>4;
		if (op == BAM_CSOFT_CLIP)break;

		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
		{
			*left_anchor_read -= ol;
			*left_anchor_ref -=ol;
		}
		else if (op == BAM_CINS) *left_anchor_read -= ol;
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP ) *left_anchor_ref -= ol;
		if (VERBOSE) fprintf(stderr,"updating leftanchor %d %d f1 %d %d%c\n",*left_anchor_read,*left_anchor_ref,f1,ol,INT_CIGAROP[op]);

		if ( op == BAM_CEQUAL && (ol >= MINLEN || (f1 ==0 && ol >= MINLEN))) // found = cigar operation of length at least MINLEN or it is the first cigar operation, beginning of read
		{
			// don't forget +1 to strlen for end character
			anchor_start = *left_anchor_ref +ol-MINLEN - COMPLEXITY_K;
			anchor_end   = anchor_start + MINLEN + COMPLEXITY_K;
			if (anchor_start < 0)	anchor_start = 0;
			if (anchor_end >= reflist->lengths[reflist->current]) anchor_end = reflist->lengths[reflist->current];
			for (j = anchor_start; j < anchor_end; j++) anchor_seq[j-anchor_start] = reflist->sequences[reflist->current][j-1];
			anchor_seq[j-anchor_start]  ='\0';

			//kmer_test = test_complexity(anchor_seq, COMPLEXITY_K);
			//if (kmer_test >=0)
			{
				flag =0;
				*left_anchor_read += ol-MINLEN; 
				*left_anchor_ref += ol-MINLEN;  // move to the right, we need only MINLEN portion
			}
		}
		f1--;
	}

	if (flag){ // didn't find a left anchor
		if (VERBOSE) fprintf(stderr,"DID NOT FIND left-anchor\n");
		free(anchor_seq);
		return 0;
	}
	else
	{
		if (VERBOSE) fprintf(stderr," found left anchor %d-%d|\n",*left_anchor_read,*left_anchor_ref);
		free(anchor_seq);
		return 1;
	}
}

// rightboundary is the rightmost position of variant, different from left_position_variant for indels
int find_right_anchor(int* fcigarlist,int fcigs,int f2,REFLIST* reflist,int offset,int* right_anchor_read, int* right_anchor_ref,int rightboundary)
{
	char* anchor_seq = malloc(MINLEN+1);
	int anchor_start = 0, anchor_end = 0, j = 0;
	//int kmer_test=0;
	int flag = 1;
	int op=fcigarlist[f2]&0xf, ol=fcigarlist[f2]>>4;
	int delta=0;
	if (VERBOSE) fprintf(stderr,"op %c ol %d offset %d..... right anchor Rboundary: %d\n",INT_CIGAROP[op],ol,offset,rightboundary);

	if (ol-offset >= MINLEN && op == BAM_CEQUAL) //  offset = last_variant_pos-right_anchor_ref, ol-offset > MINLEN = DONE 
	{
		flag =0;
		if (VERBOSE) fprintf(stderr,"right M cigar ol-offset %d offset %d\n",ol-offset,offset);
		if (ol < offset+MINLEN) 
		{
			*right_anchor_read += ol; *right_anchor_ref += ol;
		}
		else 
		{
			*right_anchor_read += offset+MINLEN; *right_anchor_ref += offset+MINLEN;
		}
	}

	int f2_bak = f2; 
	while (f2 < fcigs && flag ==1)
	{
		op = fcigarlist[f2]&0xf;
		ol = fcigarlist[f2]>>4;

		if (op == BAM_CSOFT_CLIP){
			break;
		}else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			*right_anchor_read += ol;
			*right_anchor_ref += ol;
		}
		else if (op == BAM_CINS){
			*right_anchor_read += ol;
		}
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
			*right_anchor_ref += ol;
		}
		if (VERBOSE) fprintf(stderr,"updating rightanchor %d %d f2 %d %d%c\n",*right_anchor_read,*right_anchor_ref,f2,ol,INT_CIGAROP[op]);

		if (f2 > f2_bak && op == BAM_CEQUAL && ol >= MINLEN && *right_anchor_ref - rightboundary >= MINLEN/2+2) 
		// f2 > f2_bak is needed to avoid selecting the first 'f2' which was already checked in the previous if loop
		{		
			// don't forget +1 to strlen for end character
			anchor_start = *right_anchor_ref - ol - COMPLEXITY_K;
			anchor_end = anchor_start + MINLEN + COMPLEXITY_K;
			if (anchor_start < 0)
				anchor_start = 0;
			if (anchor_end >= reflist->lengths[reflist->current])
				anchor_end = reflist->lengths[reflist->current];

			for (j = anchor_start; j < anchor_end; j++){
				anchor_seq[j - anchor_start] = reflist->sequences[reflist->current][j-1];
			}
			anchor_seq[j-anchor_start]  ='\0';

			//kmer_test = test_complexity(anchor_seq, COMPLEXITY_K);
			//if (kmer_test >=0)
			{
				flag =0;
				*right_anchor_read -= ol-MINLEN;
				*right_anchor_ref -= ol-MINLEN;
				delta = *right_anchor_ref-rightboundary-(MINLEN/2+2);
				if (delta < 0)
				{
					if (VERBOSE) fprintf(stderr,"need to move anchor to right by delta %d\n",delta);
					*right_anchor_read -= delta; 
					*right_anchor_ref -= delta;
				}
			}
		}
		f2++;
	}

	free(anchor_seq);
	if (flag){ // didn't find a right anchor
		if (VERBOSE) fprintf(stderr,"DID NOT FIND right-anchor\n");
		return 0;
	}
	else
	{
		if (VERBOSE) fprintf(stderr,"found right-anchor %d-%d\n",*right_anchor_read,*right_anchor_ref);
		return 1;
	}
}

/*
   f1 and f2 are indexes in fcigarlist containing the variant set, f1 is the index for the left of the first variant, f2 is the index of cigar op containing the last variant
   algorithm moves left from f1 (f1--) and right from f2 (f2++) to find the left and right anchor sequences
   we want to find anchor seq of length at least MINLEN
 */
int compare_read_HAPs(struct alignedread* read,VARIANT* varlist,int* snplst, int n_variants, int hap_pos[], int* fcigarlist,int fcigs,int f1, int f2, REFLIST* reflist, FRAGMENT* fragment){

	// left anchor position on read, left pos on ref, right anchor position on read, right position on ref...
	int left_anchor_read = hap_pos[0]; int left_anchor_ref = hap_pos[1]; 
	int right_anchor_read = hap_pos[2]; int right_anchor_ref = hap_pos[3]; 
	int ss,s, i;

	int rightboundary = 0;
	for (s=0;s<n_variants;s++) 
	{
		ss = snplst[s];
		if (varlist[ss].position + varlist[ss].shift > rightboundary) rightboundary = varlist[ss].position + varlist[ss].shift;
		if (VERBOSE) fprintf(stderr,"var %d %d %s %s shift %d\n",ss,varlist[ss].position,varlist[ss].allele1,varlist[ss].allele2,varlist[ss].shift);
	}
	if (VERBOSE)
	{
		//for (i=f1-2;i<=f2+2;i++) fprintf(stderr,"%d%c ",fcigarlist[i]>>4,INT_CIGAROP[fcigarlist[i]&0xf]); 
		fprintf(stderr," cigar f1=%d f2=%d L %d--%d  R %d--%d\n",f1,f2,left_anchor_read,right_anchor_read,left_anchor_ref,right_anchor_ref);
	}
	// call function to find left anchor
	int offsetL = varlist[snplst[0]].position-left_anchor_ref; // hap_pos[1] is the start of the 'f1' cigar operation
	//int offsetR = varlist[snplst[n_variants-1]].position-hap_pos[4];
	int offsetR = rightboundary-hap_pos[4];
	int foundL = find_left_anchor(fcigarlist,fcigs,f1,reflist,offsetL,&left_anchor_read,&left_anchor_ref);
	// call function to find right anchor
	if (VERBOSE) fprintf(stderr,"offset %d lastpos %d nvar %d\n",offsetR,varlist[snplst[n_variants-1]].position,n_variants);
	int foundR = find_right_anchor(fcigarlist,fcigs,f2,reflist,offsetR,&right_anchor_read,&right_anchor_ref,rightboundary);
	if (foundL ==0 || foundR ==0) 
	{
		if (VERBOSE && foundL==0) fprintf(stderr,"no left anchor for block %d\n",varlist[snplst[0]].position);
		if (VERBOSE && foundR==0) fprintf(stderr,"no right anchor for block %d\n",varlist[snplst[0]].position);
		return 0;
	}

	int positions[4] = {left_anchor_read,left_anchor_ref,right_anchor_read,right_anchor_ref}; 
	if (positions[0] < 0)positions[0] = 0;
	if (positions[1] < 0)	positions[1] = 0;
	if (positions[2] > read->readlength-1) positions[2] = read->readlength-1;
	if (positions[3] > reflist->lengths[reflist->current]-1) positions[3] = reflist->lengths[reflist->current]-1;

	for (i = 0; i < n_variants; i++){
		ss = snplst[i];
		assert(varlist[ss].position >= positions[1]);
		assert(varlist[ss].position <= positions[3]);
	}
	ss = snplst[n_variants-1];
	assert(reflist->current >= 0 && (varlist[ss].position > positions[1] && varlist[ss].position < positions[3]));
	if (VERBOSE) fprintf(stderr,"calling realign_HAPS, %d-%d %d-%d\n",positions[0],positions[2],positions[1],positions[3]);
	// positions is a 4-tuple that defines realignment region between read and haplotype
	realign_HAPs(read,reflist,positions,varlist, snplst, n_variants, fragment);

	return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// identify the variants that overlap the read and call the local realignment function (compare_read_HAP) for "clusters of variants" (1 to MAX_SNP) 

int realign_and_extract_variants_read(struct alignedread* read,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,int paired,FRAGMENT* fragment,int chrom,REFLIST* reflist)
{
	int* snplst = malloc(MAX_SNPs_SHORT_HAP*2*sizeof(int));
	int start = read->position; int end = start + read->span; int ss=0,firstvar=0,j=0,ov=0, i=0;
	// int k=0, has_a_SNV = 0;
	j = (int)(start/BSIZE);
	if (j >= chromvars[chrom].blocks) return 0; 
	ss = chromvars[chrom].intervalmap[j];
	if (ss < 0 || ss >= VARIANTS) return 0;

	// check if ss is less than first variant for the chromosome 'chrom', if so assign it to the first variant
	if (ss < chromvars[chrom].first) ss = chromvars[chrom].first;
	if (varlist[ss].position <= end)
	{
		while(ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position < start) ss++;
		firstvar = ss;
		while (ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position <= end) // BUG fixed 12/03/18, VARIANTS-1 => VARIANTS
		{
			ov++; ss++;
		}
	}
	if ((paired ==0 && ov < 2 && SINGLEREADS ==0) || (paired ==0 && ov < 1 && SINGLEREADS ==1) || (paired ==1 && ov < 1)) return 0;
	ss = firstvar; // use variable firstvar to store first variant that overlaps this read

	int fcigs= parse_cigar(read,reflist,fcigarlist);
	if (VERBOSE)
	{
		//for (i=0;i<fcigs;i++) fprintf(stdout,"%d%c ",fcigarlist[i]>>4,INT_CIGAROP[fcigarlist[i]&0xf]); fprintf(stdout," readstart %d\n",read->position);
	}

	int l1=0,l2=read->position; // l1 is advance on read, l2 is advance on reference genome
	int hap_pos[5] = {-1,-1,-1,-1,-1};
	// hap_pos[0] is advance on read for first variant in haplotype
	// hap_pos[1] is advance on ref for first variant in haplotype
	// hap_pos[2] is advance on read for last variant in haplotype
	// hap_pos[3] is advance on ref for last variant in haplotype
        // hap_pos[4] is ...

	int op=0,ol=0;
	int n_variants = 0; int n_snvs=0;
	int prev_snp_position = -1;
	int f1=0,f2=0;
	int len_a1 = 0, len_a2 = 0;
	int left_on_read = 0;

	for (i=0;i<fcigs;i++)
	{
		op = fcigarlist[i]&0xf; ol = fcigarlist[i]>>4;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)left_on_read += ol;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)left_on_read += ol;
	}

	if(VERBOSE) fprintf(stderr,"\n\n.........processing read start-pos=%d...%d overlapping %d variants, firstvar=%d \n",read->position,end,ov,firstvar+1);
	for (i=0;i<fcigs;i++)
	{
		while (ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position < l2) ss++;
		if (ss > chromvars[chrom].last || ss >= VARIANTS) break;
		op = fcigarlist[i]&0xf; ol = fcigarlist[i]>>4;

		//if (VERBOSE) fprintf(stderr,"read %d %d%c l1 %d l2 %d \n",i,ol,INT_CIGAROP[op],l1,l2);

		// BUG: encountered problem where insertion of size I, with a SNV at pos z < l2+ol,
		// would accidenally try to be handled even though the insertion doesn't actually reach the SNV
		// currently restricting to non-insertions
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) //|| op == BAM_CINS)
		{
			len_a1 = strlen(varlist[ss].allele1);
			len_a2 = strlen(varlist[ss].allele2);
			while (ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position >= l2 && varlist[ss].position < l2 + ol
					&& left_on_read > len_a1 + MINLEN && left_on_read > len_a2 + MINLEN){ // so that the read is long enough to span an indel

				if (varlist[ss].heterozygous == '1' || HOMOZYGOUS ==1){
					//fprintf(stderr,"%d %s",varlist[ss].position, varlist[ss].allele1, varlist[ss].allele2);
					// If this variant is far away from the last variant, then analyze the cluster of variants seen up til now
					if (n_variants > 0 && ((varlist[ss].position - prev_snp_position > SHORT_HAP_CUTOFF) || (n_variants >= MAX_SNPs_SHORT_HAP))){

						// if we aren't parsing INDELs, make sure this short haplotype has at least one SNV
						// call allelotyping function
						if (n_snvs > 0 || PARSEINDELS) compare_read_HAPs(read,varlist,snplst,n_variants,hap_pos,fcigarlist,fcigs,f1,f2,reflist,fragment);
						// empty the current cluster of variants
						n_variants = 0; n_snvs=0;
						for (j=0; j < 5; j++)hap_pos[j] = -1;
					}

					if (n_variants == 0){
						hap_pos[0] = l1;
						hap_pos[1] = l2;
						f1 = i;
					}

					// currently just adding the indel value to the ends of the boundaries
					// the ends might not match up perfectly, ideally we want to maintain a position we have to reach,
					// given the (possibly) large indels we've seen, and continue
					// parsing the CIGAR in the normal way until we've well passed that position
					len_a1 = strlen(varlist[ss].allele1);
					len_a2 = strlen(varlist[ss].allele2);
					if (len_a1 > len_a2){ // deletion
						if (l1 + len_a1 > hap_pos[2]) hap_pos[2] = l1 + len_a1;
						if (l2 + len_a1 > hap_pos[3]) hap_pos[3] = l2 + len_a1;
						prev_snp_position = varlist[ss].position + len_a1;
					}else{  // insertion
						if (l1 + len_a2 > hap_pos[2]) hap_pos[2] = l1 + len_a2;
						if (l2 + len_a2 > hap_pos[3]) hap_pos[3] = l2 + len_a2;
						prev_snp_position = varlist[ss].position + len_a2;
					}

					f2 = i; hap_pos[4] = l2; 
					// add variant to the list
					snplst[n_variants] = ss; // no check on length of this list, not needed
					calculate_rightshift(varlist, ss, reflist); // no need to do it multiple times..
					if (varlist[ss].type ==0) n_snvs++;
					n_variants++;
				}
				ss++;
				if (ss < VARIANTS && ss <= chromvars[chrom].last){
					len_a1 = strlen(varlist[ss].allele1);
					len_a2 = strlen(varlist[ss].allele2);
				}
			}
		}

		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			l1 += ol;
			l2 += ol;
			left_on_read -= ol;
		}
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP){
			l1 += ol;
			left_on_read -= ol;
		}else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
			l2 += ol;
		}
	}
	// might have a straggler variant left over
	if(n_variants > 0 && (n_snvs > 0 || PARSEINDELS)) compare_read_HAPs(read,varlist,snplst,n_variants,hap_pos,fcigarlist,fcigs,f1,f2,reflist,fragment);
	free(snplst);

	return 0;
}
