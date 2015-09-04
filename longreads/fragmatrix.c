#include "fragmatrix.h"
#include "printhaplotypes.c"
#include "find_starting_haplotypes.c"
#include "MECscore.c"
#define MAXBUF 10000

//////////////////////////////////////// edge list is only used in the two functions below add_edges (Feb 4 2013) //////////////////////////////

void label_node(struct SNPfrags* snpfrag, int node, int comp) // DFS search routine for connected component 
{
    int i = 0;
    if (snpfrag[node].component == -1) {
        //	fprintf(stdout," called %d node edges %d %d \n",node,snpfrag[node].edges,comp);
        snpfrag[node].component = comp;
        snpfrag[comp].csize++;
        for (i = 0; i < snpfrag[node].edges; i++) label_node(snpfrag, snpfrag[node].elist[i].snp, comp);
    }
}

void add_edges_fosmids(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components) {
    int i = 0, j = 0, t = 0, k = 0, iter = 0, maxdeg = 0, avgdeg = 0, mdelta = 0, t1 = 0, t2 = 0, l = 0;
    char allele;
    int csnps = 0;
    for (i = 0; i < snps; i++) snpfrag[i].edges = 0;

    int varlist[4096];
    char allelelist[4096];
    int vars = 0;
    for (i = 0; i < fragments; i++) {
        // generate list of variants and alleles 
        vars = 0;
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                varlist[vars] = Flist[i].list[j].offset + k;
                allelelist[vars] = Flist[i].list[j].hap[k];
                vars++;
                if (vars >= 4096) break;
            }
        }
        // add edge between adjacent pair of variants in each fragment 
        for (j = 0; j < vars - 1; j++) {
            t1 = varlist[j];
            t2 = varlist[j + 1];
            snpfrag[t1].elist[snpfrag[t1].edges].snp = t2;
            snpfrag[t1].elist[snpfrag[t1].edges].frag = i;
            snpfrag[t1].elist[snpfrag[t1].edges].p[0] = allelelist[j];
            snpfrag[t1].elist[snpfrag[t1].edges].p[1] = allelelist[k];
            snpfrag[t1].edges++;
            snpfrag[t2].elist[snpfrag[t2].edges].snp = t1;
            snpfrag[t2].elist[snpfrag[t2].edges].frag = i;
            snpfrag[t2].elist[snpfrag[t2].edges].p[1] = allelelist[j];
            snpfrag[t2].elist[snpfrag[t2].edges].p[0] = allelelist[k];
            snpfrag[t2].edges++;
        }

    }
    // elist contains duplicates (due to multiple fragments), telist does not, feb 5 2013 
    // sort all edges lists once for all by snp number, this can be done faster using QSORT, see later code...
    for (i = 0; i < snps; i++) qsort(snpfrag[i].elist, snpfrag[i].edges, sizeof (struct edge), edge_compare);
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].edges > maxdeg) maxdeg = snpfrag[i].edges;
        avgdeg += snpfrag[i].frags;
        if (snpfrag[i].edges == 0) continue;
        csnps++;
        if (snpfrag[i].component != -1) continue; // already labeled with component
        snpfrag[i].component = i;
        for (j = 0; j < snpfrag[i].edges; j++) label_node(snpfrag, snpfrag[i].elist[j].snp, i);
    }
    for (i = 0; i < fragments; i++) Flist[i].component = snpfrag[Flist[i].list[0].offset].component; // each fragment has a component fixed 

    *components = 0;
    int nodes_in_graph = 0;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1) {
            (*components)++;
            nodes_in_graph += snpfrag[i].csize;
        }
        //else if (snpfrag[i].component ==i || snpfrag[i].edges ==0) singletons++;
    }
    fprintf(stdout, "Number of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
    //fprintf(stderr, "Number of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
}

// for each fragment: add all pairwise edges between all variants in it, complexity = O(k^2) for 'k' length fragment

void add_edges(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components) {
    int i = 0, j = 0, t = 0, k = 0, iter = 0, maxdeg = 0, avgdeg = 0, mdelta = 0;
    int csnps = 0;
    for (i = 0; i < snps; i++) snpfrag[i].edges = 0;
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                for (t = 0; t < Flist[i].blocks; t++) {
                    for (iter = 0; iter < Flist[i].list[t].len; iter++) {
                        if (Flist[i].list[j].offset + k == Flist[i].list[t].offset + iter) continue;
                        if (Flist[i].list[j].offset + k - Flist[i].list[t].offset + iter > mdelta) mdelta = Flist[i].list[j].offset + k - Flist[i].list[t].offset + iter;
                        snpfrag[Flist[i].list[t].offset + iter].elist[snpfrag[Flist[i].list[t].offset + iter].edges].snp = Flist[i].list[j].offset + k;
                        snpfrag[Flist[i].list[j].offset + k].elist[snpfrag[Flist[i].list[j].offset + k].edges].frag = i;
                        snpfrag[Flist[i].list[t].offset + iter].elist[snpfrag[Flist[i].list[t].offset + iter].edges].frag = i;
                        snpfrag[Flist[i].list[t].offset + iter].elist[snpfrag[Flist[i].list[t].offset + iter].edges].p[0] = Flist[i].list[t].hap[iter];
                        snpfrag[Flist[i].list[t].offset + iter].elist[snpfrag[Flist[i].list[t].offset + iter].edges].p[1] = Flist[i].list[j].hap[k];
                        snpfrag[Flist[i].list[t].offset + iter].edges++;
                    }
                }
            }
        }
    }
    // elist contains duplicates (due to multiple fragments), telist does not, feb 5 2013 
    // sort all edges lists once for all by snp number, this can be done faster using QSORT, see later code...
    for (i = 0; i < snps; i++) qsort(snpfrag[i].elist, snpfrag[i].edges, sizeof (struct edge), edge_compare);
    for (i = 0; i < snps; i++) {
        //fprintf(stdout," snp %d edges %d || ",i,snpfrag[i].edges); for (j=0;j<snpfrag[i].edges;j++) fprintf(stdout,"%d ",snpfrag[i].elist[j]); fprintf(stdout,"\n"); getchar();
        if (snpfrag[i].edges > maxdeg) maxdeg = snpfrag[i].edges;
        avgdeg += snpfrag[i].frags;
        // edit here june 7 2012
        if (snpfrag[i].edges == 0) continue;
        csnps++;
        if (snpfrag[i].component != -1) continue; // already labeled with component
        snpfrag[i].component = i;
        for (j = 0; j < snpfrag[i].edges; j++) label_node(snpfrag, snpfrag[i].elist[j].snp, i);
    }
    for (i = 0; i < fragments; i++) Flist[i].component = snpfrag[Flist[i].list[0].offset].component; // each fragment has a component fixed 

    *components = 0;
    int nodes_in_graph = 0;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1) {
            (*components)++;
            nodes_in_graph += snpfrag[i].csize;
        }
        //else if (snpfrag[i].component ==i || snpfrag[i].edges ==0) singletons++;
    }
    fprintf(stdout, "\nno of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
    fprintf(stderr, "\nno of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
}
//////////////////////////////////////// edge list is only used in the two functions below //////////////////////////////

// populate the connected component data structure only for non-trivial connected components, at least 2 variants

void generate_clist_structure(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int components, struct BLOCK* clist) {
    // bcomp maps the SNP to the component number in clist since components << snps  
    // note that the previous invariant about the first SNP (ordered by position) being the root of the component is no longer true !! 
    // should we still require it ?? 
    int i = 0, j = 0, component = 0;
    for (i = 0; i < snps; i++) snpfrag[i].bcomp = -1;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1) // root node of component
        {
            snpfrag[i].bcomp = component;
            clist[component].slist = calloc(sizeof (int), snpfrag[i].csize);
            clist[component].phased = 0;
            component++;
        }
    }
    //fprintf(stderr,"non-trivial components in graph %d \n",components);
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component < 0) continue; // to allow for initialization to -1 in other code feb 15 2013
        if (snpfrag[i].csize <= 1 && snpfrag[i].component == i) continue; // ignore singletons that are not connected to other variants

        if (snpfrag[i].component != i) snpfrag[i].bcomp = snpfrag[snpfrag[i].component].bcomp;
        if (snpfrag[i].bcomp < 0) continue;
        component = snpfrag[i].bcomp;
        if (clist[component].phased == 0) clist[component].offset = i;
        clist[component].slist[clist[component].phased] = i;
        clist[component].phased++;
        clist[component].lastvar = i;
    }
    for (i = 0; i < components; i++) clist[i].length = clist[i].lastvar - clist[i].offset + 1;
    for (i = 0; i < components; i++) clist[i].frags = 0;
    for (i = 0; i < fragments; i++) {
        if (snpfrag[Flist[i].list[0].offset].bcomp < 0)continue; // ignore fragments that cover singleton vertices
        clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++;
    }
    for (i = 0; i < components; i++) clist[i].flist = calloc(sizeof (int), clist[i].frags);
    for (i = 0; i < components; i++) clist[i].frags = 0;
    for (i = 0; i < fragments; i++) {
        if (snpfrag[Flist[i].list[0].offset].bcomp < 0)continue;
        clist[snpfrag[Flist[i].list[0].offset].bcomp].flist[clist[snpfrag[Flist[i].list[0].offset].bcomp].frags] = i;
        clist[snpfrag[Flist[i].list[0].offset].bcomp].frags++;
    }
    for (i = 0; i < components; i++) {
        fprintf(stdout, "comp %d first %d last %d phased %d fragments %d \n", i, clist[i].offset, clist[i].lastvar, clist[i].phased, clist[i].frags);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this function updates the data structure snpfrag which links the heterozyous SNPs and the haplotype fragments
// // it generates a list of fragments (flist) that affect each SNP 

void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components) {
    int i = 0, j = 0, t = 0, k = 0, iter = 0, calls = 0; //maxdeg=0,avgdeg=0;

    // find the first fragment whose endpoint lies at snp 'i' or beyond
    for (i = 0; i < snps; i++) {
        snpfrag[i].frags = 0;
        snpfrag[i].ff = -1;
        snpfrag[i].split = 0;
        snpfrag[i].prune_status = 0;
        snpfrag[i].post_notsw = 0;
    }
    for (i = 0; i < fragments; i++) {
        j = Flist[i].list[0].offset;
        k = Flist[i].list[Flist[i].blocks - 1].len + Flist[i].list[Flist[i].blocks - 1].offset;
        // commented the line below since it slows program for long mate-pairs june 7 2012
        for (t=j;t<k;t++) { if (snpfrag[t].ff == -1) snpfrag[t].ff = i;  } 
    } //for (i=0;i<snps;i++) { fprintf(stdout,"SNP %d firstfrag %d start snp %d \n",i,snpfrag[i].ff,i); } 

    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            //if (Flist[i].list[j].offset+k < 0 || Flist[i].list[j].offset+k >= snps) fprintf(stdout,"%d %d %s\n",Flist[i].list[j].offset,snps,Flist[i].id);
            for (k = 0; k < Flist[i].list[j].len; k++) snpfrag[Flist[i].list[j].offset + k].frags++;
        }
    }
    for (i = 0; i < snps; i++) {
        snpfrag[i].flist = (int*) malloc(sizeof (int)*snpfrag[i].frags);
        snpfrag[i].alist = (char*) malloc(snpfrag[i].frags);
    }

    for (i = 0; i < snps; i++) {
        snpfrag[i].component = -1;
        snpfrag[i].csize = 1;
        snpfrag[i].frags = 0;
        snpfrag[i].edges = 0;
        snpfrag[i].best_mec = 10000;
    }

    for (i = 0; i < fragments; i++) {
        calls = 0;
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                snpfrag[Flist[i].list[j].offset + k].flist[snpfrag[Flist[i].list[j].offset + k].frags] = i;
                snpfrag[Flist[i].list[j].offset + k].alist[snpfrag[Flist[i].list[j].offset + k].frags] = Flist[i].list[j].hap[k];
                snpfrag[Flist[i].list[j].offset + k].frags++;
                calls += Flist[i].list[j].len;
            }
        }

        if (FOSMIDS == 1) // long reads 
        {
            // 2 edges for every node in fragment ( 00000----0---[1]----10011 ) adjacent left and right
            for (j = 0; j < Flist[i].blocks; j++) {
                for (k = 0; k < Flist[i].list[j].len; k++) snpfrag[Flist[i].list[j].offset + k].edges += 2;
            }
            j = 0;
            snpfrag[Flist[i].list[j].offset].edges -= 1;
            j = Flist[i].blocks - 1;
            snpfrag[Flist[i].list[j].offset + Flist[i].list[j].len - 1].edges -= 1;
            // single edge for first and last node in fragment 
        } else {
            for (j = 0; j < Flist[i].blocks; j++) {
                for (k = 0; k < Flist[i].list[j].len; k++) snpfrag[Flist[i].list[j].offset + k].edges += calls - 1;
            }
        }
    }
}
