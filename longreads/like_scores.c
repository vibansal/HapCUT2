/* functions to calculate likelihoods P(read| haplotype) for sequencing errors and chimeric fragments */

// need function that calculates best score using a single switch error vs no switch error in likelihood space  
//added 03/03/2015 also calculates minimum score over bitflips + switch errors...

void calculate_fragscore(struct fragment* Flist, int f, char* h, float* mec_ll, float* chimeric_ll) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    float chim_prob = -1000000;
    int bit = 0, bits = 0;

    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10; // log10(e)
            prob1 = 1.0 - pow(10, prob);
            prob2 = Flist[f].list[j].p1[k];
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }
            bit += 1; // counter over the alleles of the fragment, ignore invalid alleles 
        }
    }
    bits = bit;
    bit = 0;

    if (bits > 2) {
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (h[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
                prob = QVoffset - (int) Flist[f].list[j].qv[k];
                prob /= 10; // log10(e)
                prob1 = 1.0 - pow(10, prob);
                prob2 = Flist[f].list[j].p1[k];

                if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                    p0 -= prob2;
                    p0 += prob;
                    p1 -= prob;
                    p1 += prob2;
                } else {
                    p0 -= prob;
                    p0 += prob2;
                    p1 -= prob2;
                    p1 += prob;
                }
                if (bit > 0 && bit < bits - 1) // add the switch error likelihood to chim_prob
                {
                    if (p0 > chim_prob) chim_prob = p0 + log10(1.0 + pow(10, chim_prob - p0));
                    else chim_prob += log10(1.0 + pow(10, p0 - chim_prob));
                    if (p1 > chim_prob) chim_prob = p1 + log10(1.0 + pow(10, chim_prob - p1));
                    else chim_prob += log10(1.0 + pow(10, p1 - chim_prob));
                }
                bit += 1;
            }
        }
        chim_prob -= log10(bits - 2);
    }
    *chimeric_ll = chim_prob;

    if (p0 > p1) *mec_ll = (p0 + log10(1 + pow(10, p1 - p0)));
    else *mec_ll = (p1 + log10(1 + pow(10, p0 - p1)));
    //return mec_prob; 
}

void update_fragscore(struct fragment* Flist, int f, char* h) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    Flist[f].calls = 0;
    float good = 0, bad = 0;
    int switches = 0;
    int m = 0;
    if (h[Flist[f].list[0].offset] == Flist[f].list[0].hap[0]) m = 1;
    else m = -1; // initialize 
    for (j = 0; j < Flist[f].blocks; j++) {
        Flist[f].calls += Flist[f].list[j].len;
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
            prob2 = Flist[f].list[j].p1[k];

            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) good += prob1;
            else bad += prob1;
            //if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
            // this is likelihood based calculation 
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k] && m == -1) {
                m = 1;
                switches++;
            } else if (h[Flist[f].list[j].offset + k] != Flist[f].list[j].hap[k] && m == 1) {
                m = -1;
                switches++;
            }
        }
    }
    if (p0 > p1) Flist[f].ll = (p0 + log10(1 + pow(10, p1 - p0)));
    else Flist[f].ll = (p1 + log10(1 + pow(10, p0 - p1)));

    if (SCORING_FUNCTION == 0) {
        if (good < bad) Flist[f].currscore = good;
        else Flist[f].currscore = bad;
    } else if (SCORING_FUNCTION = 5) Flist[f].currscore = -1 * Flist[f].ll;
    else if (SCORING_FUNCTION == 1) Flist[f].currscore = switches;
    else {
        if (switches == 2 && (good <= 1.0 || bad <= 1.0)) Flist[f].currscore += 1;
        else if (switches == 4 && (good <= 2.0 || bad <= 2.0)) Flist[f].currscore += 2;
        else Flist[f].currscore = switches;
    }
}

// calculate probability of observing k=0,1,2 seq errors in the read 03/04/2015

int calculate_error_probs(struct fragment* Flist, int f, char* h, float perr[], int max) {
    float perror[Flist[f].calls][max];
    int j = 0, k = 0, t = 0;
    float prob = 0, prob1 = 0;
    int bit = 0;
    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;

            if (bit == 0) perror[bit][0] = Flist[f].list[j].p1[k];
            else perror[bit][0] = perror[bit - 1][0] + Flist[f].list[j].p1[k];

            for (t = 1; t < max; t++) {
                if (bit == 0 && t == 1) perror[bit][t] = prob;
                else if (bit >= t - 1) perror[bit][t] = perror[bit - 1][t - 1] + prob;
                if (bit >= t) {
                    perror[bit][t] += log10(1 + pow(10, perror[bit - 1][t] + Flist[f].list[j].p1[k] - perror[bit][t]));
                }
            }
            bit++;
        }
    }
    for (t = 0; t < max; t++) perr[t] = perror[bit - 1][t];
    return bit;
}

// compute score of fragment
// don't mutate Flist or other structures.
// return score as return value
// homozygous: 0-based index of a homozygous position. -1 if no homozygous pos

float simple_fragscore(struct fragment* Flist, int f, char* h, int homozygous) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    float good = 0, bad = 0, ll;

    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
            prob2 = Flist[f].list[j].p1[k];

            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) good += prob1;
            else bad += prob1;

            // this is likelihood based calculation 
            if (Flist[f].list[j].offset + k != homozygous) { // not homozygous
                if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                    p0 += prob2;
                    p1 += prob;
                } else {
                    p0 += prob;
                    p1 += prob2;
                }
            } else { // homozygous at this postion
                if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                    p0 += prob2; // both hap1 and hap2 match
                    p1 += prob2;
                } else {
                    p0 += prob;
                    p1 += prob;
                }
            }
        }
    }
    if (p0 > p1) ll = (p0 + log10(1 + pow(10, p1 - p0)));
    else ll = (p1 + log10(1 + pow(10, p0 - p1)));

    return ll;
}