
/*
 float p1 = q1*q2+(1-q1)*(1-q2); float p2 = q1*(1-q2)+q2*(1-q1);
        if (hap[i] == hap[j] && p[0] == p[1]) return log10(p1/p2);
        else if (hap[i] != hap[j] && p[0] != p[1]) return log10(p1/p2);
        else if (hap[i] == hap[j] && p[0] != p[1]) return log10(p2/p1);
        else if (hap[i] != hap[j] && p[0] == p[1]) return log10(p2/p1);
 */

// algorithm to pick a neighbor of a 'node' for seeding max-cut algorithm, pick negative wt edges 

int select_neighbor(struct SNPfrags* snpfrag, int* slist, int N, int node, struct fragment* Flist, char* hap) {
    int j = 0, f = 0, i = 0, n = 0, k = 0;
    char allele, allele1;
    float q2, q1, p1, p2, ew, match;

    for (i = 0; i < N; i++) snpfrag[slist[i]].score = 0;
    for (n = 0; n < snpfrag[node].frags; n++) {
        f = snpfrag[node].flist[n]; // index into Flist global 
        for (j = 0; j < Flist[f].blocks; j++) // first find the allele at 'node' in the fragment 'f' 
        {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (Flist[f].list[j].offset + k == node) {
                    allele = Flist[f].list[j].hap[k];
                    q1 = Flist[f].list[j].pv[k];
                }
            }
        }
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                allele1 = Flist[f].list[j].hap[k];
                q2 = Flist[f].list[j].pv[k];
                if (node == Flist[f].list[j].offset + k) continue;
                p1 = q1 * q2 + (1 - q1)*(1 - q2);
                p2 = q1 * (1 - q2) + q2 * (1 - q1);
                ew = log10(p1 / p2);
                if (hap[node] == hap[Flist[f].list[j].offset + k] && allele == allele1) match = 1; //log10(p1/p2);
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele != allele1) match = 1; //log10(p1/p2);
                else if (hap[node] == hap[Flist[f].list[j].offset + k] && allele != allele1) match = -1; //*log10(p1/p2);
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele == allele1) match = -1; //*log10(p1/p2);
                snpfrag[Flist[f].list[j].offset + k].score += match*ew;
            }
        }
    }
    float best = 1000, sbest = 1000;
    int b0 = -1, b1 = -1;
    for (i = 0; i < N; i++) {
        if (snpfrag[slist[i]].score < best) {
            if (b1 < 0) {
                sbest = best;
                b1 = b0;
            }
            best = snpfrag[slist[i]].score;
            b0 = slist[i];
        } else if (snpfrag[slist[i]].score < sbest) {
            sbest = snpfrag[slist[i]].score;
            b1 = slist[i];
        }
    }
    if (b1 < 0 || drand48() < 0.75) return b0;
    else return b1;
}

// function added 10/24/2014

void init_neighbor_scores(struct SNPfrags* snpfrag, int node, struct fragment* Flist, char* hap, float weight) {
    // for all vertices linked to node -> update score of the neighbors, can be done using list of fragments for snpfrag[node]
    int j = 0, f = 0, n = 0, k = 0;
    char allele, allele1;
    float q2, q1, p1, p2, ew, match;
    for (n = 0; n < snpfrag[node].frags; n++) {
        f = snpfrag[node].flist[n]; // index into Flist global 
        for (j = 0; j < Flist[f].blocks; j++) // first find the allele at 'node' in the fragment 'f' 
        {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (Flist[f].list[j].offset + k == node) {
                    allele = Flist[f].list[j].hap[k];
                    q1 = Flist[f].list[j].pv[k];
                }
            }
        }
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                allele1 = Flist[f].list[j].hap[k];
                q2 = Flist[f].list[j].pv[k];
                if (node == Flist[f].list[j].offset + k) continue;
                p1 = q1 * q2 + (1 - q1)*(1 - q2);
                p2 = q1 * (1 - q2) + q2 * (1 - q1);
                ew = log10(p1 / p2);
                if (hap[node] == hap[Flist[f].list[j].offset + k] && allele == allele1) match = 1; //log10(p1/p2);
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele != allele1) match = 1; //log10(p1/p2);
                else if (hap[node] == hap[Flist[f].list[j].offset + k] && allele != allele1) match = -1; //*log10(p1/p2);
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele == allele1) match = -1; //*log10(p1/p2);
                snpfrag[Flist[f].list[j].offset + k].score += match * weight*ew;
            }
        }
    }
}

// function added 10/24/2014

void update_neighbor_scores(struct SNPfrags* snpfrag, int node, int startnode, int secondnode, struct fragment* Flist, char* hap, struct PHEAP* pheap, int* slist) {
    // for all vertices linked to node -> update score of the neighbors, can be done using list of fragments for snpfrag[node]
    int j = 0, f = 0, n = 0, k = 0, newnode = 0;
    char allele, allele1;
    float q2, q1, p1, p2, ew, match;
    double oldscore = 0;
    for (n = 0; n < snpfrag[node].frags; n++) {
        f = snpfrag[node].flist[n]; // index into Flist global 
        for (j = 0; j < Flist[f].blocks; j++) // first find the allele at 'node' in the fragment 'f' 
        {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (Flist[f].list[j].offset + k == node) {
                    allele = Flist[f].list[j].hap[k];
                    q1 = Flist[f].list[j].pv[k];
                }
            }
        }
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                allele1 = Flist[f].list[j].hap[k];
                q2 = Flist[f].list[j].pv[k];
                newnode = Flist[f].list[j].offset + k;
                if (newnode == node || newnode == startnode || newnode == secondnode) continue;

                p1 = q1 * q2 + (1 - q1)*(1 - q2);
                p2 = q1 * (1 - q2) + q2 * (1 - q1);
                ew = log10(p1 / p2);
                if (hap[node] == hap[Flist[f].list[j].offset + k] && allele == allele1) match = 1;
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele != allele1) match = 1;
                else if (hap[node] == hap[Flist[f].list[j].offset + k] && allele != allele1) match = -1;
                else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele == allele1) match = -1;
                if (snpfrag[node].parent == startnode) {
                    oldscore = snpfrag[newnode].score;
                    snpfrag[newnode].score += match*ew;
                } else if (snpfrag[node].parent == secondnode) {
                    oldscore = snpfrag[newnode].score;
                    snpfrag[newnode].score -= match*ew;
                }
                if (fabsf(oldscore) > fabsf(snpfrag[newnode].score)) // score decreased 
                {
                    pmaxHeapify(pheap, snpfrag[newnode].heaploc, snpfrag, slist);
                } else pbubbleUp(pheap, snpfrag[newnode].heaploc, snpfrag, slist);
            }
        }
    }
}

float cut_MEC(struct SNPfrags* snpfrag, struct fragment* Flist, char* hap, int* slist, int N, struct BLOCK* component) {

    int i = 0, j = 0, f = 0;
    float newscore = 0;
    for (i = 0; i < N; i++) {
        // flip the alleles at sites for which parent == 0 (one shore of cut) 
        if (snpfrag[slist[i]].parent == 0 && hap[slist[i]] == '1') hap[slist[i]] = '0';
        else if (snpfrag[slist[i]].parent == 0 && hap[slist[i]] == '0') hap[slist[i]] = '1';
    }

    for (j = 0; j < component->frags; j++) {
        f = component->flist[j];
        update_fragscore(Flist, f, hap);
        newscore += Flist[f].currscore;
    }

    // revert back
    for (i = 0; i < N; i++) {
        if (snpfrag[slist[i]].parent == 0 && hap[slist[i]] == '1') hap[slist[i]] = '0';
        else if (snpfrag[slist[i]].parent == 0 && hap[slist[i]] == '0') hap[slist[i]] = '1';
    }

    for (j = 0; j < component->frags; j++) {
        f = component->flist[j];
        update_fragscore(Flist, f, hap);
    }
    return newscore;
}

float find_cutvalue(struct SNPfrags* snpfrag, struct fragment* Flist, char* hap, int* slist, int N) {
    int j = 0, f = 0, i = 0, n = 0, k = 0, node = 0, newnode = 0;
    char allele, allele1;
    float q2, q1, p1, p2, ew, match;
    float curr_cut = 0;

    for (i = 0; i < N; i++) {
        node = slist[i];
        // for all vertices linked to node -> update score of the neighbors, can be done using list of fragments for snpfrag[node]
        for (n = 0; n < snpfrag[node].frags; n++) {
            f = snpfrag[node].flist[n]; // index into Flist global 
            for (j = 0; j < Flist[f].blocks; j++) // first find the allele at 'node' in the fragment 'f' 
            {
                for (k = 0; k < Flist[f].list[j].len; k++) {
                    if (Flist[f].list[j].offset + k == node) {
                        allele = Flist[f].list[j].hap[k];
                        q1 = Flist[f].list[j].pv[k];
                    }
                }
            }
            for (j = 0; j < Flist[f].blocks; j++) {
                for (k = 0; k < Flist[f].list[j].len; k++) {
                    allele1 = Flist[f].list[j].hap[k];
                    q2 = Flist[f].list[j].pv[k];
                    newnode = Flist[f].list[j].offset + k;
                    if (newnode == node) continue;
                    p1 = q1 * q2 + (1 - q1)*(1 - q2);
                    p2 = q1 * (1 - q2) + q2 * (1 - q1);
                    ew = log10(p1 / p2);
                    if (hap[node] == hap[Flist[f].list[j].offset + k] && allele == allele1) match = 1;
                    else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele != allele1) match = 1;
                    else if (hap[node] == hap[Flist[f].list[j].offset + k] && allele != allele1) match = -1;
                    else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele == allele1) match = -1;
                    if (snpfrag[node].parent != snpfrag[newnode].parent) curr_cut += match * 0.5 * ew;
                    // multiply by 1/2 since we double count pairs of edges 
                }
            }
        }
    }
    return curr_cut;
}

// function not working, going in loop ?? 10/25/14
// function to see if moving individual SNP from one end to other end of cut improves cut score...
// reason it is not working is that cut-score + node-score are not on same scale

float improve_cutvalue(struct SNPfrags* snpfrag, struct fragment* Flist, char* hap, int* slist, int N, float* curr_cut, int* c1, int* c2) {
    int j = 0, f = 0, i = 0, n = 0, k = 0, node = 0, newnode = 0, moved = 0;
    char allele, allele1;
    float q2, q1, p1, p2, ew, match;

    while (1) {
        moved = 0;
        for (i = 0; i < N; i++) {
            node = slist[i];
            snpfrag[node].score = 0;
            for (n = 0; n < snpfrag[node].frags; n++) {
                f = snpfrag[node].flist[n]; // index into Flist global 
                for (j = 0; j < Flist[f].blocks; j++) // first find the allele at 'node' in the fragment 'f' 
                {
                    for (k = 0; k < Flist[f].list[j].len; k++) {
                        if (Flist[f].list[j].offset + k == node) {
                            allele = Flist[f].list[j].hap[k];
                            q1 = Flist[f].list[j].pv[k];
                        }
                    }
                }
                for (j = 0; j < Flist[f].blocks; j++) {
                    for (k = 0; k < Flist[f].list[j].len; k++) {
                        allele1 = Flist[f].list[j].hap[k];
                        q2 = Flist[f].list[j].pv[k];
                        newnode = Flist[f].list[j].offset + k;
                        if (newnode == node) continue;
                        p1 = q1 * q2 + (1 - q1)*(1 - q2);
                        p2 = q1 * (1 - q2) + q2 * (1 - q1);
                        ew = log10(p1 / p2);
                        if (hap[node] == hap[Flist[f].list[j].offset + k] && allele == allele1) match = 1;
                        else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele != allele1) match = 1;
                        else if (hap[node] == hap[Flist[f].list[j].offset + k] && allele != allele1) match = -1;
                        else if (hap[node] != hap[Flist[f].list[j].offset + k] && allele == allele1) match = -1;
                        if (snpfrag[newnode].parent == 0) snpfrag[node].score += match * ew;
                        else if (snpfrag[newnode].parent == 1) snpfrag[node].score -= match * ew;
                    }
                }
                if (snpfrag[node].parent == 0 && snpfrag[node].score < 0 && *c1 > 1) {
                    snpfrag[node].parent = 1;
                    *curr_cut += snpfrag[node].score;
                    moved++;
                    (*c1)--;
                    (*c2)++;
                } else if (snpfrag[node].parent == 1 && snpfrag[node].score > 0 && *c2 > 1) {
                    snpfrag[node].parent = 0;
                    *curr_cut -= snpfrag[node].score;
                    moved++;
                    (*c2)--;
                    (*c1)++;
                }
            }
        }
        if (moved == 0) break;
    }
    return 0;
}
