/* functions for comparing an aligned sequence read to the set of variants to identify alleles and haplotype-informative reads */

#include "nw.c"
#include <assert.h>
#include <stdlib.h>

int MINLEN= 10;
int COMPLEXITY_K = 5; // anchor sequences must have unique kmers of this length
// find the variants that are covered by the read and determine the alleles at each of those variants
int SHORT_HAP_CUTOFF = 20;
extern int VERBOSE;
extern int REALIGN_ALL;
extern int PARSEINDELS;
float TINYLOG = -10000;
int MIN_QUAL = 10;
int MAX_SNPs_SHORT_HAP = 10; // max number of SNVs in a short haplotype
// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 - pow(10, (b) - (a)))) : ((b) + log10(1.0 - pow(10.0, (a) - (b)))))

// credit to Jeff Mercado at stack overflow for the original version of this function found here: https://stackoverflow.com/questions/5370753/using-stdlibs-qsort-to-sort-an-array-of-strings
int compare_strings(const void* a, const void* b)
{
    const char *ia = *(const char **)a;        //const char *ia = (const char *)a;
    const char *ib = *(const char **)b;  //const char *ib = (const char *)b

    return strcmp(ia, ib);
}

// full credit for this function goes to theJPster at stackoverflow
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
	//fprintf(stdout,"%c \n",reflist->sequences[current][10000]);	//return -1;

	int f=0;
	int i=0,t=0, l1=0,l2=0; int l=0; int m=0; int op;
	for (i=0;i<read->cigs;i++)
	{
		op = read->cigarlist[i]&0xf; l = read->cigarlist[i]>>4;
                if (op == BAM_CMATCH)
		{
			//fprintf(stdout,"%.*s | ",l,read->sequence+l1); fprintf(stdout,"%.*s\n",l,reflist->sequences[current]+read->position+l2-1);
			m=0;
			for (t=0;t<l;t++)
                        {
                                if (read->sequence[l1+t]  != reflist->sequences[current][read->position+l2+t-1] && read->sequence[l1+t] != reflist->sequences[current][read->position+l2+t-1]-32 && read->sequence[l1+t]  !='N')
				{
					read->mismatches++;
					if (m > 0) fcigarlist[f++] = (m<<4)+7; //  m=
					fcigarlist[f++] = 24; // mismatch 1X
					m=0;
				}
				else m++;
                        }
			if (m > 0) fcigarlist[f++] = (m<<4)+7;
		}
		else fcigarlist[f++] = read->cigarlist[i];  // copy as it is for some cigar operations

                if (op == BAM_CMATCH || op >= 7) { l1 +=l; l2 +=l; }
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
            if (i == j)
                continue;

            match = 1;
            for (q = 0; q < k; q++){
                if (seq[i+q] != seq[j+q]){
                    match = 0;
                    break;
                }
            }

            if (match){
                //fprintf(stderr," fail\n");
                return 0;
            }
        }
    }
    //fprintf(stderr," pass\n");

    return 1;
}

// make sure that the code works for indels/complex-vars |  illumina reads |  low-complexity regions | GENERIC CODE

int realign_HAPs(struct alignedread* read, REFLIST* reflist, int positions[], VARIANT* varlist, int* snplst, int n_snps, FRAGMENT* fragment)
{
	if ( positions[3]-positions[1] < 15 ||  positions[3]-positions[1] > 200) return -1;
	if ( positions[2]-positions[0] < 15 ||  positions[2]-positions[0] > 200) return -1;

	int i=0,j=0,k=0;

	char* subread = malloc(positions[2]-positions[0]+1);

	if (positions[0] < 0 || positions[2]-1 >= strlen(read->sequence)){
		fprintf(stderr,"BOUNDS CHECK FAIL: read %s, read length %d, start %d, end %d\n", read->readid, read->readlength, positions[0], positions[2]-1);
		exit(1);
	}

	for (j=positions[0];j<positions[2];j++){
        subread[j-positions[0]] = read->sequence[j];
    }
    subread[j-positions[0]]='\0';

	char* refhap = malloc(positions[3]-positions[1]+1);

    // don't forget +1 to strlen for end character
    for (j=positions[1];j<positions[3];j++){
        refhap[j-positions[1]] = reflist->sequences[reflist->current][j-1];
    }
    refhap[j-positions[1]]  ='\0';

	char* althap;
	int h=0, s=0, ss=0, max_hap=0, ref_len=0, alt_len=0, total_ref_len=0, total_alt_len=0, n_max_haps = 0, rand_ix = 0;
	double total_score = TINYLOG, max_score = -1000000000;
    int align_qual = 0;
	int* max_haps = (int*) calloc(pow(2,n_snps),sizeof(int));

	double* ref_score_single = malloc(MAX_SNPs_SHORT_HAP * sizeof(double));
	double* alt_score_single = malloc(MAX_SNPs_SHORT_HAP * sizeof(double));
	for (s = 0; s < MAX_SNPs_SHORT_HAP; s++){
		// initialize to very small log value
		ref_score_single[s] = TINYLOG;
		alt_score_single[s] = TINYLOG;
	}

    if (VERBOSE) fprintf(stderr,subread);
    if (VERBOSE) fprintf(stderr,"\n");
	if (VERBOSE) fprintf(stderr,"read %d:%d:%d ",positions[0],positions[2],positions[2]-positions[0]);
	if (VERBOSE) fprintf(stderr,"on ref %d:%d:%d\n",positions[1],positions[3],positions[3]-positions[1]);
	for (h = 0; h < pow(2,n_snps); h++){

		// we represent a haplotype with integer h
		// in particular the n_snps least significant bits
		// if the i-th bit from the right is 1, computed as h & (pow(2,i)), then the haplotype contains the i-th variant

		total_ref_len = 0;
		total_alt_len = 0;
		for (s = 0; s < n_snps; s++){
			ss = snplst[s];
			total_ref_len += strlen(varlist[ss].allele1);
			if (h & (int)(pow(2,s))){
				total_alt_len += strlen(varlist[ss].allele2);
			}else{
				total_alt_len += strlen(varlist[ss].allele1);
			}
		}

		althap = malloc(positions[3]-positions[1]+1+total_alt_len-total_ref_len);

		j = positions[1];
		k = 0;

		for (s = 0; s < n_snps; s++){

			ss = snplst[s];
			ref_len = strlen(varlist[ss].allele1);
			alt_len = strlen(varlist[ss].allele2);

			for (j = j; j<varlist[ss].position; j++)
				althap[k++] = reflist->sequences[reflist->current][j-1];
			if (h & (int)(pow(2,s))){
				for (i=0;i<alt_len;i++) althap[k++] = varlist[ss].allele2[i];
			}else{
				for (i=0;i<ref_len;i++) althap[k++] = varlist[ss].allele1[i];
			}

			j += ref_len;
		}

		for (j = j; j < positions[3]; j++){
			althap[k++] = reflist->sequences[reflist->current][j-1];
		}

		althap[k] = '\0';

		if (VERBOSE) fprintf(stdout,"hap%d %s ",h,althap);

		double altscore = nw(althap,subread,0);

		if (VERBOSE) fprintf(stdout,"score: %f \n",altscore);

		// for an index s in the short haplotype,
		// maintain the log sum of scores that have a variant at s
		// and the log sum of scores that have reference at s
		for (s = 0; s < n_snps; s++){

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
			memset(max_haps, 0, pow(2,n_snps)*sizeof(int));
		}
		// in case of a tie we store every max haplotype so we can randomly select one
		if(altscore == max_score){
			max_haps[n_max_haps] = h;
			n_max_haps++;
		}

		// add the reference score and alternate score to the total
		//total_score = addlogs(total_score, refscore);
		total_score = addlogs(total_score, altscore);

        free(althap);
	}
	if (VERBOSE) fprintf(stdout,"**********************************************\n");


	// if the max score is sufficiently high then write alleles to haplotype fragments and calculate quality values from scores

	// TODO
	// need to assign each SNVs' score to (sum of haplotypes with variant)/(sum of all haplotypes)
	rand_ix = rand_interval(0,n_max_haps-1);
	max_hap = max_haps[rand_ix];

	for (s = 0; s < n_snps; s++){

		ss = snplst[s];

        if (varlist[ss].type != 0 && !PARSEINDELS){
            continue;
        }

		fragment->alist[fragment->variants].varid = ss;

		if (max_hap & (int)(pow(2,s))){
			align_qual = (int) (-10.0 * (ref_score_single[s] - total_score));

			if (align_qual < MIN_QUAL) continue;

			fragment->alist[fragment->variants].allele = '1';
		}else{
			align_qual = (int) (-10.0 * (alt_score_single[s] - total_score));

			if (align_qual < MIN_QUAL) continue;

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

	free(subread);free(refhap); free(max_haps);free(ref_score_single); free(alt_score_single);

    return 0;
}

int compare_read_HAPs(struct alignedread* read,VARIANT* varlist,int* snplst, int n_snps, int hap_pos[], int* fcigarlist,int fcigs,int f1, int f2, REFLIST* reflist, FRAGMENT* fragment){
	int op = fcigarlist[f1+1]&0xf; int ol = fcigarlist[f1+1]>>4;
	int offset = varlist[snplst[0]].position-hap_pos[1];
    char* anchor_seq = malloc(MINLEN+1);
    int anchor_start = 0, anchor_end = 0, j = 0;

	int positions[4] = {hap_pos[0],hap_pos[1],hap_pos[2],hap_pos[3]}; // left anchor position on read, left pos on ref, right position on read, right position on ref...

	int flag = 1;
	if (offset >= MINLEN && (fcigarlist[f1]&0xf) == BAM_CMATCH)
	{
		flag =0; // don't need to change positions[0] or positions[1], could be too large ??
	}

	while (f1 >= 0 && flag ==1){

		op = fcigarlist[f1]&0xf; ol = fcigarlist[f1]>>4;
        if (VERBOSE) fprintf(stdout,"%d%c | ",ol,INT_CIGAROP[op]);
		if (op == BAM_CSOFT_CLIP){
			break;
		}else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			positions[0] -= ol;
			positions[1] -=ol;
		}else if (op == BAM_CINS){
			positions[0] -= ol;
		}else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
			positions[1] -= ol;
		}

		if (op == BAM_CEQUAL && ol >= MINLEN){

            // don't forget +1 to strlen for end character
            anchor_start = positions[1]+ol-MINLEN - COMPLEXITY_K;
            anchor_end   = anchor_start + MINLEN + COMPLEXITY_K;
            if (anchor_start < 0)
                anchor_start = 0;
            if (anchor_end >= reflist->lengths[reflist->current])
                anchor_end = reflist->lengths[reflist->current];

            for (j = anchor_start; j < anchor_end; j++){
                anchor_seq[j-anchor_start] = reflist->sequences[reflist->current][j-1];
            }
            anchor_seq[j-anchor_start]  ='\0';

            if (test_complexity(anchor_seq, COMPLEXITY_K)){
    			flag =0;
                positions[0] += ol-MINLEN;
    			positions[1] += ol-MINLEN;
            }
		}

		f1--;
	}

    if (flag){ // didn't find a left anchor
        free(anchor_seq);
        return 0;
    }

	if (VERBOSE) fprintf(stdout," left |");

	flag = 1;
	if ((hap_pos[3] +(fcigarlist[f2]>>4)-varlist[snplst[n_snps-1]].position >= MINLEN) && (fcigarlist[f2]&0xf) == BAM_CMATCH)
	{
		flag =0;
		positions[2] += fcigarlist[f2]>>4;
		positions[3] += fcigarlist[f2]>>4;
	}

	int f2_bak = f2;
	while (f2 < fcigs && flag ==1){

		op = fcigarlist[f2]&0xf;
		ol = fcigarlist[f2]>>4;
		if (VERBOSE) fprintf(stdout,"%d%c | ",ol,INT_CIGAROP[op]);

		if (op == BAM_CSOFT_CLIP){
			break;
		}else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			positions[2] += ol;
			positions[3] += ol;
		}
		else if (op == BAM_CINS){
			positions[2] += ol;
		}
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP){
			positions[3] += ol;
		}

		if (f2 > f2_bak && op == BAM_CEQUAL && ol >= MINLEN){

            // don't forget +1 to strlen for end character
            anchor_start = positions[3] - ol - COMPLEXITY_K;
            anchor_end = anchor_start + MINLEN + COMPLEXITY_K;
            if (anchor_start < 0)
                anchor_start = 0;
            if (anchor_end >= reflist->lengths[reflist->current])
                anchor_end = reflist->lengths[reflist->current];

            for (j = anchor_start; j < anchor_end; j++){
                anchor_seq[j - anchor_start] = reflist->sequences[reflist->current][j-1];
            }
            anchor_seq[j-anchor_start]  ='\0';

            if (test_complexity(anchor_seq, COMPLEXITY_K)){
    			flag =0;
                positions[2] -= ol-MINLEN;
    			positions[3] -= ol-MINLEN;
            }
		}
		f2++;
	}

    if (flag){ // didn't find a right anchor
        free(anchor_seq);
        return 0;
    }

	if (VERBOSE) fprintf(stdout,"found right-anchor\n");

	if (positions[0] < 0){
		positions[0] = 0;
	}
	if (positions[1] < 0){
		positions[1] = 0;
	}
	if (positions[2] > read->readlength-1){
		positions[2] = read->readlength-1;
	}
	if (positions[3] > reflist->lengths[reflist->current]-1){
		positions[3] = reflist->lengths[reflist->current]-1;
	}
	int ss, i;
	for (i = 0; i < n_snps; i++){
		ss = snplst[i];
		assert(varlist[ss].position >= positions[1]);
		assert(varlist[ss].position <= positions[3]);
	}

	ss = snplst[n_snps-1];
	assert(reflist->current >= 0 && (varlist[ss].position > positions[1] && varlist[ss].position < positions[3]));

	realign_HAPs(read,reflist,positions,varlist, snplst, n_snps, fragment);

    free(anchor_seq);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int realign_and_extract_variants_read(struct alignedread* read,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,int paired,FRAGMENT* fragment,int chrom,REFLIST* reflist){
	//fprintf(stderr,"%s \n",read->readid);

    int* snplst = malloc(MAX_SNPs_SHORT_HAP*sizeof(int));
	int start = read->position; int end = start + read->span; int ss=0,firstvar=0,j=0,ov=0, i=0, k=0, has_a_SNV = 0;
	j = (int)(start/BSIZE);
	if (j >= chromvars[chrom].blocks) return 0; // another BUG april29 2011 found here
	ss = chromvars[chrom].intervalmap[j];
	if (ss < 0 || ss >= VARIANTS) return 0;

	// check if ss is less than first variant for the chromosome 'chrom', if so assign it to the first variant
	if (ss < chromvars[chrom].first) ss = chromvars[chrom].first;

	if (varlist[ss].position <= end)
	{
		while(ss < VARIANTS-1 && ss <= chromvars[chrom].last && varlist[ss].position < start){
            ss++;
        }
        firstvar = ss;
		while (ss < VARIANTS-1 && ss <= chromvars[chrom].last && varlist[ss].position <= end)
		{
			ov++; ss++;
		}
	}
	if ((paired ==0 && ov < 2 && SINGLEREADS ==0) || (paired ==0 && ov < 1 && SINGLEREADS ==1) || (paired ==1 && ov < 1)) return 0;
	ss = firstvar; // use variable firstvar to store first variant that overlaps this read
//	fprintf(stderr,"parsing read %s \n",read->readid);

	int fcigs=0;
	fcigs = parse_cigar(read,reflist,fcigarlist);
	//fprintf(stdout,"getting full cigar %d %d \n",read->cigs,fcigs);
	//for (i=0;i<read->cigs;i++) fprintf(stdout,"%d%c ",read->cigarlist[i]>>4,INT_CIGAROP[read->cigarlist[i]&0xf]); fprintf(stdout,"\n");
	//for (i=0;i<fcigs;i++) fprintf(stdout,"%d%c ",fcigarlist[i]>>4,INT_CIGAROP[fcigarlist[i]&0xf]); fprintf(stdout," readstart %d\n",read->position);
	//exit(0);

	int l1=0,l2=read->position; // l1 is advance on read, l2 is advance on reference genome
	int hap_pos[4] = {-1,-1,-1,-1};
	// hap_pos[0] is advance on read for first variant in haplotype
	// hap_pos[1] is advance on ref for first variant in haplotype
	// hap_pos[2] is advance on read for last variant in haplotype
	// hap_pos[3] is advance on ref for last variant in haplotype

	int op=0,ol=0;
        //fprintf(stdout,"%d\n",fcigs);
	int n_snps = 0;
	int prev_snp_position = -1;
	int f1=0,f2=0;
    int len_a1 = 0, len_a2 = 0;
    int left_on_read = 0;

    for (i=0;i<fcigs;i++)
	{
		op = fcigarlist[i]&0xf;
        ol = fcigarlist[i]>>4;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
			left_on_read += ol;
		}else if (op == BAM_CINS || op == BAM_CSOFT_CLIP)
			left_on_read += ol;
	}

	for (i=0;i<fcigs;i++)
	{
		while (ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position < l2){
			ss++;
		}
		if (ss > chromvars[chrom].last || ss >= VARIANTS){
			break;
		}
		op = fcigarlist[i]&0xf; ol = fcigarlist[i]>>4;

		// BUG: encountered problem where insertion of size I, with a SNV at pos z < l2+ol,
		// would accidenally try to be handled even though the insertion doesn't actually reach the SNV
		// currently restricting to non-insertions
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL){ //|| op == BAM_CINS){

            len_a1 = strlen(varlist[ss].allele1);
            len_a2 = strlen(varlist[ss].allele2);
			while (ss < VARIANTS && ss <= chromvars[chrom].last && varlist[ss].position >= l2 && varlist[ss].position < l2 + ol
                   && left_on_read > len_a1 + MINLEN && left_on_read > len_a2 + MINLEN){ // so that the read is long enough to span an indel

				if (varlist[ss].heterozygous == '1'){
                    //fprintf(stderr,"%d %s",varlist[ss].position, varlist[ss].allele1, varlist[ss].allele2);
					// If this variant is far away from the last variant, then analyze the cluster of variants seen up til now
					if (n_snps > 0 && ((varlist[ss].position - prev_snp_position > SHORT_HAP_CUTOFF) || (n_snps == MAX_SNPs_SHORT_HAP))){

                        // if we aren't parsing INDELs, make sure this short haplotype has at least one SNV
                        has_a_SNV = 0;
                        for (k = 0; k < n_snps; k++){
                            if (varlist[snplst[k]].type == 0){
                                has_a_SNV = 1;
                            }
                        }

                        if (has_a_SNV || PARSEINDELS){
							compare_read_HAPs(read,varlist,snplst,n_snps,hap_pos,fcigarlist,fcigs,f1,f2,reflist,fragment);
                        }

						// empty the current cluster of variants
						n_snps = 0;

						for (j=0; j < 4; j++){
							hap_pos[j] = -1;
						}
					}

					if (n_snps == 0){
						hap_pos[0] = l1;
						hap_pos[1] = l2;
						f1 = i;
					}

                    //hap_pos[2] = l1;
                    //hap_pos[3] = l2;

                    // currently just adding the indel value to the ends of the boundaries
                    // this is not an ideal solution
                    // the ends might not match up perfectly
                    // ideally we want to maintain a position we have to reach,
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

					f2 = i;

					// add variant to the list
					snplst[n_snps] = ss;

					n_snps++;
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

	// might have a straggler SNP left over
	if(n_snps > 0){

        // if we aren't parsing INDELs, make sure this short haplotype has at least one SNV
        has_a_SNV = 0;
        for (k = 0; k < n_snps; k++){
            if (varlist[snplst[k]].type == 0){
                has_a_SNV = 1;
            }
        }

        //
        if (has_a_SNV || PARSEINDELS){
            compare_read_HAPs(read,varlist,snplst,n_snps,hap_pos,fcigarlist,fcigs,f1,f2,reflist,fragment);
        }
	}

    free(snplst);

	return 0;
}
