
#include "like_scores.c" // additional functions for likelihood based calculation

// scoring function =5 -> use -1*LL instead of MEC function 

// function to compute weight of a edge between two variants in the max-cut graph (score of +1 if the phasing agrees with the phasing suggested by the fragment and score of -1 if it does not)
// this function needs to be modified to utilize base-quality scores, so that the weight of an edge is (1-e1)(1-e2) for the two base-calls

float edge_weight(char* hap, int i, int j, char* p, struct fragment* Flist, int f) {
    float q1 = 1, q2 = 1;
    int k = 0, l = 0;
    // new code added so that edges are weighted by quality scores, running time is linear in length of fragment !! 08/15/13 | reduce this
    for (k = 0; k < Flist[f].blocks; k++) {
        for (l = 0; l < Flist[f].list[k].len; l++) {
            if (Flist[f].list[k].offset + l == i) q1 = Flist[f].list[k].pv[l];
            else if (Flist[f].list[k].offset + l == j) q2 = Flist[f].list[k].pv[l];
        }
    }
    float p1 = q1 * q2 + (1 - q1)*(1 - q2);
    float p2 = q1 * (1 - q2) + q2 * (1 - q1);
    if (hap[i] == hap[j] && p[0] == p[1]) return log10(p1 / p2);
    else if (hap[i] != hap[j] && p[0] != p[1]) return log10(p1 / p2);
    else if (hap[i] == hap[j] && p[0] != p[1]) return log10(p2 / p1);
    else if (hap[i] != hap[j] && p[0] == p[1]) return log10(p2 / p1);
    else return 0;
}

// compute MECSCORE of the fragment matrix compared to a haplotype 'h' 

int mecscore(struct fragment* Flist, int fragments, char* h, float* ll, float* calls, float* miscalls) {
    float trueMEC = 0;
    //	fprintf(stderr,"QVoffset is now %d \n",QVoffset);
    int j = 0, k = 0, f = 0;
    *ll = 0;
    float p0, p1;
    *calls = 0;
    *miscalls = 0;
    float prob, prob1, prob2, L;
    float good = 0, bad = 0;
    int switches = 0;
    int bitflips = 0;
    int m = 0;
    int lastswitch = -10;
    int bit = 0;
    Flist[f].vbits = 0;
    for (f = 0; f < fragments; f++) {
        good = bad = 0;
        p0 = p1 = 0; //if (Flist[f].blocks ==1 && Flist[f].list[0].len ==1) continue;
        if (h[Flist[f].list[0].offset] == Flist[f].list[0].hap[0]) m = 1;
        else m = -1; // initialize 
        switches = 0;

        for (j = 0; j < Flist[f].blocks; j++) {
            *calls += Flist[f].list[j].len;
            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (h[Flist[f].list[j].offset + k] == '-') continue;
                if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
                prob = QVoffset - (int) Flist[f].list[j].qv[k];
                prob /= 10; // log10(e)
                prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
                prob2 = Flist[f].list[j].p1[k];

                //if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good +=1; else bad +=1;
                if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) good += prob1;
                else bad += prob1;

                if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k] && m == -1) {
                    if (lastswitch == bit - 1 && bit > 0) {
                        bitflips++;
                        switches--;
                    } else {
                        switches++;
                        lastswitch = bit;
                    }
                    m = 1;
                } else if (h[Flist[f].list[j].offset + k] != Flist[f].list[j].hap[k] && m == 1) {
                    if (lastswitch == bit - 1 && bit > 0) {
                        bitflips++;
                        switches--;
                    } else {
                        switches++;
                        lastswitch = bit;
                    }
                    m = -1;
                }

                prob1 = log10(prob1);
                //printf("prob %s %f %f \n",Flist[f].list[j].qv,prob,prob1); 
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
        if (p0 > p1) L = (p0 + log10(1 + pow(10, p1 - p0)));
        else L = (p1 + log10(1 + pow(10, p0 - p1)));
        *ll += L;
        if (SCORING_FUNCTION == 0) {
            if (good < bad) *miscalls += good;
            else *miscalls += bad;
        } else if (SCORING_FUNCTION == 5) *miscalls += -1 * L;
        else if (SCORING_FUNCTION == 1) *miscalls += switches + bitflips;
        else {
            if (switches == 2 && (good <= 1.0 || bad <= 1.0)) *miscalls += 1;
            else if (switches == 4 && (good <= 2.0 || bad <= 2.0)) *miscalls += 2;
            else *miscalls += switches + bitflips;
        }
        //fprintf(stdout,"good %f bad %f frag %d %f\n",good,bad,f,*calls);
        if (good < bad) trueMEC += good;
        else trueMEC += bad;
    }
    Flist[f].vbits = bit;
    //fprintf(stderr,"true MEC %f \n",trueMEC);
    if (SCORING_FUNCTION == 5) return ((int) trueMEC);
    return *calls;
}

// function to print comparison of fragment to haplotype | format is 45333 000:010:BAF  (fragment block):(hap block):(qscore block)

int print_fragment_MEC(struct fragment* Flist, int f, char* h, FILE* outfile) {
    int j = 0, k = 0, t = 0;
    float prob = 0, prob1 = 0, good = 0, bad = 0, mec = 0;
    float fl[3] = {0, 0, 0};
    float sep = -3; // switch error probability
    int hap = 0;
    int switches = 0;
    int m = 0;
    int mm = 0;
    int bit = 0;
    int bitflips = 0;
    int lastswitch = -10;
    if (h[Flist[f].list[0].offset] == Flist[f].list[0].hap[0]) m = 1;
    else m = -1; // initialize 

    static float expected_counts[5] = {0, 0, 0, 0, 0};
    static int fragments = 0;
    float perr[4];
    int bits = calculate_error_probs(Flist, f, h, perr, 4);
    for (t = 0; t < 4; t++) {
        if (t < bits) expected_counts[t] += pow(10, perr[t]);
    }
    fragments += 1;
    //fprintf(outfile,"frags %d counts %0.2f %0.2f %0.2f %0.2f\n",fragments,expected_counts[0],expected_counts[1],expected_counts[2],expected_counts[3]);

    float lastprob1;

    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob1 = 1.0;
            prob1 -= pow(10, prob);
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) good += prob1;
            else bad += prob1;
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                fl[0] += Flist[f].list[j].p1[k];
                fl[1] += prob;
            } else {
                fl[0] += prob;
                fl[1] += Flist[f].list[j].p1[k];
            }

            // don't allow switch at first base of fragment, only bitflip 
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k] && m == -1) {
                m = 1;
                if (bit == 1 || (bit == bits - 1)) {
                    bitflips++;
                    fl[3] += lastprob1;
                } else if (lastswitch == bit - 1 && bit > 0) {
                    bitflips++;
                    switches--;
                    fl[2] += lastprob1;
                } else {
                    switches++;
                    lastswitch = bit;
                }
            } else if (h[Flist[f].list[j].offset + k] != Flist[f].list[j].hap[k] && m == 1) {
                m = -1;
                if (bit == 1 || (bit == bits - 1)) {
                    bitflips++;
                    fl[2] += lastprob1;
                } else if (lastswitch == bit - 1 && bit > 0) {
                    bitflips++;
                    switches--;
                    fl[2] += lastprob1;
                } else {
                    switches++;
                    lastswitch = bit;
                }
            } else fl[2] += Flist[f].list[j].p1[k];
            bit++;
            lastprob1 = prob;
        }
    }
    fl[2] += switches * sep * log10(bits);
    if (good < bad) {
        mec = good;
        hap = 0;
    } else {
        mec = bad;
        hap = 1;
    }

    if (mec <= 0.001) return 0;
    if (mec <= 1.01) fprintf(outfile, "M1E ");
    else if (switches == 1 && bitflips + switches < mec - 1) fprintf(outfile, "SWE ");
    else fprintf(outfile, "MEC ");

    // only interested in fragments with multiple MECs
    fprintf(outfile, "%0.1f %d %d %d ", mec, switches, bitflips, Flist[f].blocks);
    fprintf(outfile, "EP %0.2f %0.2f %0.2f ", perr[1], perr[2], perr[3]);
    fprintf(outfile, "FL %0.2f %0.2f %0.2f ", fl[0], fl[1], fl[2]);
    fprintf(outfile, "%s %d ", Flist[f].id, Flist[f].list[0].offset);

    for (j = 0; j < Flist[f].blocks; j++) {
        mm = 0;
        //fprintf(outfile,"%d ",Flist[f].list[j].offset); 
        for (k = 0; k < Flist[f].list[j].len; k++) fprintf(outfile, "%c", Flist[f].list[j].hap[k]);
        fprintf(outfile, ":");
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            if (hap == 1 && h[Flist[f].list[j].offset + k] != Flist[f].list[j].hap[k]) mm++;
            else if (hap == 0 && h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) mm++;

            if (hap == 1) fprintf(outfile, "%c", h[Flist[f].list[j].offset + k]);
            else if (h[Flist[f].list[j].offset + k] == '0') fprintf(outfile, "1");
            else if (h[Flist[f].list[j].offset + k] == '1') fprintf(outfile, "0");
            else fprintf(outfile, "-");
        }
        if (mm > 0) {
            fprintf(outfile, ":");
            for (k = 0; k < Flist[f].list[j].len; k++) fprintf(outfile, "%d,", (int) Flist[f].list[j].qv[k] - 33);
        }
        fprintf(outfile, " ");
    }
    fprintf(outfile, "\n");
    return 0;
}

void print_fragmentmatrix_MEC(struct fragment* Flist, int fragments, char* h, char* outfileprefix) {
    fprintf(stderr, "printing fragment matrix along with MEC scores to a file \n");
    char outfile[1024];
    sprintf(outfile, "%s.fragments", outfileprefix);
    FILE* fp = fopen(outfile, "w");
    int i = 0;
    for (i = 0; i < fragments; i++) print_fragment_MEC(Flist, i, h, fp);
    fclose(fp);
}
