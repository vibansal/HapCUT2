extern int HIC;

/* functions to calculate likelihoods P(read| haplotype) for sequencing errors and chimeric fragments */

void calculate_fragscore1(struct fragment* Flist, int f, char* h, float* ll) {
    struct fragment* FRAG = &Flist[f];
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    float p0h = 0, p1h = 0, normal_ll=0, htrans_ll=0; // p0 and p1 for h-trans

    for (j = 0; j < FRAG->blocks; j++) {
        for (k = 0; k < FRAG->list[j].len; k++) {

            if (h[FRAG->list[j].offset + k] == '-' || (int) FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
            prob = (QVoffset - (int) FRAG->list[j].qv[k]); prob /= 10;
            prob2 = FRAG->list[j].p1[k];

            if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }

            if (HIC && FRAG->data_type == 1 && FRAG->mate2_ix != -1){ // HiC read
                if (FRAG->list[j].offset + k >= FRAG->mate2_ix){
                    // we are on mate 2 so haplotypes (and probabilities) are flipped, since h-trans
                    if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]){
                        p0h += prob;
                        p1h += prob2;
                    } else {
                        p0h += prob2;
                        p1h += prob;
                    }
                } else {
                    // we are not on mate 2 so the LL calculation is normal
                    if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]){
                        p0h += prob2;
                        p1h += prob;
                    } else {
                        p0h += prob;
                        p1h += prob2;
                    }
                }
            }

        }
    }

    if (HIC && FRAG->data_type == 1 && FRAG->mate2_ix != -1){ // HiC read

        // Hi-C LL calculation accounting for h-trans

        normal_ll = addlogs(p0, p1);   // probability of read assuming h-cis
        htrans_ll = addlogs(p0h, p1h); // probability of read assuming h-trans

        normal_ll += subtractlogs(0,FRAG->htrans_prob);  // multiply the normal read LL calcuation by probability read is normal
        htrans_ll += FRAG->htrans_prob;                  // multiply the h-trans read LL calculation by probability read is h-trans

        *ll = addlogs(normal_ll, htrans_ll);
    }else{
        // normal LL calculation

        *ll = addlogs(p0,p1);

    }
}

void update_fragscore1(struct fragment* Flist,int f, char* h) {
    struct fragment* FRAG = &Flist[f];
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    float p0h = 0, p1h = 0; // p0 and p1 for h-trans interactions
    float normal_ll, htrans_ll;
    FRAG->calls = 0;
    float good = 0, bad = 0;
    int switches = 0;
    int m = 0;
    if (h[FRAG->list[0].offset] == FRAG->list[0].hap[0]) m = 1;
    else m = -1; // initialize
    for (j = 0; j < FRAG->blocks; j++) {
        FRAG->calls += FRAG->list[j].len;
        for (k = 0; k < FRAG->list[j].len; k++) {
            if (h[FRAG->list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) FRAG->list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) FRAG->list[j].qv[k];
            prob /= 10;
            prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
            prob2 = FRAG->list[j].p1[k];

            if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) good += prob1;
            else bad += prob1;
            //if (h[FRAG->list[j].offset+k] == FRAG->list[j].hap[k]) good++; else bad++;
            // this is likelihood based calculation
            if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }

            // this is the likelihood based calculation assuming a Hi-C h-trans interaction
            if (HIC && FRAG->data_type == 1){
                if (FRAG->mate2_ix != -1 && FRAG->list[j].offset + k >= FRAG->mate2_ix){
                    // we are on mate 2 so haplotypes (and probabilities) are flipped, since h-trans
                    if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]){
                        p0h += prob;
                        p1h += prob2;
                    } else {
                        p0h += prob2;
                        p1h += prob;
                    }
                } else {
                    // we are not on mate 2 so the LL calculation is normal
                    if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]){
                        p0h += prob2;
                        p1h += prob;
                    } else {
                        p0h += prob;
                        p1h += prob2;
                    }
                }
            }

            if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k] && m == -1) {
                m = 1;
                switches++;
            } else if (h[FRAG->list[j].offset + k] != FRAG->list[j].hap[k] && m == 1) {
                m = -1;
                switches++;
            }
        }
    }

    if (HIC == 0 || FRAG->data_type == 0){ // normal fragment
        // normal LL calculation
        FRAG->ll = addlogs(p0,p1);
    } else if (HIC && FRAG->data_type == 1){ // HiC fragment 
        // Hi-C LL calculation accounting for h-trans
        normal_ll = addlogs(p0,p1);
        htrans_ll = addlogs(p0h,p1h);

        normal_ll += subtractlogs(0,FRAG->htrans_prob);  // multiply the normal read LL calcuation by probability read is normal
        htrans_ll += FRAG->htrans_prob;                   // multiply the h-trans read LL calculation by probability read is h-trans

        FRAG->ll = addlogs(normal_ll, htrans_ll);
    }

    FRAG->currscore = -1 * FRAG->ll;
}

// compute score of fragment
// don't mutate Flist or other structures.
// return score as return value
// homozygous: 0-based index of a homozygous position. -1 if no homozygous pos

// switch_ix: 0-based index of the switch error being tested, -1 if none
float fragment_ll1(struct fragment* Flist, int f, char* h, int homozygous, int switch_ix) {
    struct fragment* FRAG = &Flist[f];
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, p0h = 0, p1h =0, prob = 0, prob2 = 0;
    //float prob1=0;
    float ll=0;
    int snp_ix, switched, htrans_flipped;

    // normal LL calculation, no Hi-C h-trans
    if (HIC == 0||FRAG->data_type == 0){
        for (j = 0; j < FRAG->blocks; j++) {
            for (k = 0; k < FRAG->list[j].len; k++) {
                snp_ix = FRAG->list[j].offset + k; // index of current position with respect to all SNPs

                // conditions to skip this base
                if (h[snp_ix] == '-') continue;

                prob = QVoffset - (int) FRAG->list[j].qv[k];
                prob /= 10;
                //prob1 = 1.0 - pow(10, prob);
                prob2 = FRAG->list[j].p1[k];

              //  if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) good += prob1;
               // else bad += prob1;
		if (homozygous == snp_ix && switch_ix == 10) continue; // for calculating likelihood if variant is ignored, hack added 03/13/2018

                // this is likelihood based calculation
                switched = (switch_ix != -1 && snp_ix >= switch_ix);
                if ((h[snp_ix] == FRAG->list[j].hap[k]) != switched) { // true if match, or not match but switched
                    p0 += prob2; 
	            if (snp_ix != homozygous) p1 += prob;
		    else p1 += prob2;
                } else {
                    p0 += prob; 
                    if (snp_ix != homozygous) p1 += prob2;
                    else p1 += prob;
                }
            }
        }

        ll = addlogs(p0, p1);
    }else if (HIC && FRAG->data_type == 1){ // Hi-C LL calculation that accounts for Hi-C h-trans

        for (j = 0; j < FRAG->blocks; j++) {
            for (k = 0; k < FRAG->list[j].len; k++) {
                snp_ix = FRAG->list[j].offset + k; // index of current position with respect to all SNPs

                // conditions to skip this base
                if (h[snp_ix] == '-') continue;
                htrans_flipped = (FRAG->mate2_ix != -1 && snp_ix >= FRAG->mate2_ix); // are we flipped due to h-trans?
                prob = QVoffset - (int) FRAG->list[j].qv[k];
                prob /= 10;
                //prob1 = 1.0 - pow(10, prob);
                prob2 = FRAG->list[j].p1[k];

                //if (h[FRAG->list[j].offset + k] == FRAG->list[j].hap[k]) good += prob1;
                //else bad += prob1;

                if (h[snp_ix] == FRAG->list[j].hap[k]) { // true if match, or not match but switched
                    p0 += prob2; 
                    if (snp_ix != homozygous)
                        p1 += prob;
                    else
                        p1 += prob2;
                } else {
                    p0 += prob;
                    if (snp_ix != homozygous)
                        p1 += prob2;
                    else
                        p1 += prob;
                }

                if (((h[snp_ix] == FRAG->list[j].hap[k])&&(!htrans_flipped))
                  ||((h[snp_ix] != FRAG->list[j].hap[k])&&(htrans_flipped))){ // true if match, or not match but switched
                    p0h += prob2;
                    if (snp_ix != homozygous)
                        p1h += prob;
                    else
                        p1h += prob2;
                } else {
                    p0h += prob;
                    if (snp_ix != homozygous)
                        p1h += prob2;
                    else
                        p1h += prob;
                }
            }
        }

        ll = addlogs(addlogs(p0, p1)+subtractlogs(0,FRAG->htrans_prob), addlogs(p0h, p1h)+FRAG->htrans_prob);

    }
    //if (ll > 0.31) fprintf(stdout,"ERROR, likelihood should always be less than 2, frag %d %d %f\n",f,homozygous,ll);
    return ll;
}

