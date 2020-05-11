#include "khash.h"
KHASH_SET_INIT_INT(32)
#include "variantgraph.h"
/*
functions to build the read-haplotype graph (nodes are SNPfrag objects), build edge list, find connected components
*/

extern int VERBOSE;
extern int LONG_READS;
//////////////////////////////////////// edge list is only used in the two functions below add_edges (Feb 4 2013) //////////////////////////////


void print_variant(struct SNPfrags* snpfrag,int i,FILE* OUTFILE)
{ 
 fprintf(OUTFILE,"VAR %d %d %s %s %s %d\n",i,snpfrag[i].position,snpfrag[i].allele0,snpfrag[i].allele1,snpfrag[i].genotypes,snpfrag[i].frags);
}


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


// non-recursive version of label_node
void label_node_alt(struct SNPfrags* snpfrag, int init_node, int comp, khash_t(32) *label_node_hash) {
	int ret;
	int *nodes  = NULL;
    int m_nodes = 16;
    int n_nodes = 0;

    // init
	kh_clear(32, label_node_hash);
    nodes    = (int*)calloc(m_nodes, sizeof(int));
    nodes[0] = init_node;
    n_nodes  = 1;
	kh_put(32, label_node_hash, init_node, &ret);

    while (0 < n_nodes) {
        // get the next node
        int node = nodes[n_nodes-1];
        n_nodes--;

        if (snpfrag[node].component != -1) continue;

        // process
        snpfrag[node].component = comp;
        snpfrag[comp].csize++;
        int i;
        for (i = snpfrag[node].edges - 1; 0 <= i; i--) { // reverse order for DFS
			int cur_node = snpfrag[node].elist[i].snp;
			// check the hash
			if (kh_get(32, label_node_hash, cur_node) != kh_end(label_node_hash)) continue;

            // make room

            while (m_nodes <= n_nodes) {
                m_nodes <<= 1;
				if ((1 << 29) <= m_nodes) {
					fprintf(stderr, "Too many nodes allocated: %d\n", m_nodes);
					exit(1);
				}
                nodes = (int*)realloc(nodes, sizeof(int)*m_nodes);
            }
            nodes[n_nodes] = cur_node;
            n_nodes++;
			kh_put(32, label_node_hash, cur_node, &ret);
        }
    }

    // destroy
    free(nodes);
}

void label_node(struct SNPfrags* snpfrag, int node, int comp, khash_t(32) *label_node_hash) // DFS search routine for connected component
{
	/*
    int i = 0, ret;
	if (kh_get(32, label_node_hash, node) != kh_end(label_node_hash)) return;
	kh_put(32, label_node_hash, node, &ret);
	if (ret == 0) {
		fprintf(stderr, "kh_put returned non-empty\n");
		exit (1);
	}
    if (snpfrag[node].component == -1) {
        //  fprintf(stdout," called %d node edges %d %d \n",node,snpfrag[node].edges,comp);
        snpfrag[node].component = comp;
        snpfrag[comp].csize++;
        for (i = 0; i < snpfrag[node].edges; i++) label_node(snpfrag, snpfrag[node].elist[i].snp, comp, label_node_hash);
    }
	*/
    label_node_alt(snpfrag, node, comp, label_node_hash);
}

int edge_compare(const void *a, const void *b) {
    const struct edge *ia = (const struct edge*) a;
    const struct edge *ib = (const struct edge*) b;
    return ia->snp - ib->snp;
}

void add_edges_longreads(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components) {
    int i = 0, j = 0, k = 0, maxdeg = 0, avgdeg = 0, t1 = 0, t2 = 0;
    //char allele;
    int csnps = 0;
    int max_vars = 65536;
	khash_t(32) *label_node_hash = kh_init(32);
    for (i = 0; i < snps; i++) snpfrag[i].elist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));
    for (i = 0; i < snps; i++) snpfrag[i].edges = 0;

    int varlist[max_vars];
    char allelelist[max_vars];

    for (i = 0; i < max_vars; i++){
        varlist[i] = 0;
        allelelist[i] = 0;
    }
    int vars = 0;
    for (i = 0; i < fragments; i++) {
        // generate list of variants and alleles
        vars = 0;
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
		if (snpfrag[Flist[i].list[j].offset + k].phase == '0') continue; // filter for homozygous variants 
                varlist[vars] = Flist[i].list[j].offset + k;
                allelelist[vars] = Flist[i].list[j].hap[k];
                vars++;
                if (vars >= 65536) break;
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
	//fprintf(stderr, "Iterating through SNPs\n0");
	//int every = (snps < 10000) ? snps : (snps / 10000.0);
    for (i = 0; i < snps; i++) {
		//if (0 == (i % every)) fprintf(stderr, "\r%.2lf%% complete", i * 100.0 / snps);
        if (snpfrag[i].edges > maxdeg) maxdeg = snpfrag[i].edges;
        avgdeg += snpfrag[i].frags;
        if (snpfrag[i].edges == 0) continue;
        csnps++;
        if (snpfrag[i].component != -1) continue; // already labeled with component
        snpfrag[i].component = i;
        for (j = 0; j < snpfrag[i].edges; j++) {
			kh_clear(32, label_node_hash);
			label_node(snpfrag, snpfrag[i].elist[j].snp, i, label_node_hash);
		}
    }
	///fprintf(stderr, "\r%.2lf%% complete\n", i * 100.0 / snps);
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
	kh_destroy(32, label_node_hash);
    fprintf(stdout, "Number of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
    //fprintf(stderr, "Number of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
}

// for each fragment: add all pairwise edges between all variants in it, complexity = O(k^2) for 'k' length fragment

void add_edges(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components) {
	khash_t(32) *label_node_hash = kh_init(32);
    int i = 0, j = 0, t = 0, k = 0, iter = 0, maxdeg = 0, avgdeg = 0, mdelta = 0;
    int csnps = 0;
    for (i = 0; i < snps; i++) snpfrag[i].elist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));
    for (i = 0; i < snps; i++) snpfrag[i].edges = 0;
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) {
                for (t = 0; t < Flist[i].blocks; t++) {
                    for (iter = 0; iter < Flist[i].list[t].len; iter++) {
                        if (Flist[i].list[j].offset + k == Flist[i].list[t].offset + iter) continue;
			if (snpfrag[Flist[i].list[j].offset + k].phase == '0' || snpfrag[Flist[i].list[t].offset + iter].phase == '0') continue;  // filter for homozygous
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
	//fprintf(stderr, "Iterating through SNPs\n0");
	//int every = snps / 10000.0;
    for (i = 0; i < snps; i++) {
		//if (snps > 0 && every > 0 && (0 == (i % every))){
		//	 fprintf(stderr, "\r%.2lf%% complete", i * 100.0 / snps);
		//}
        //fprintf(stdout," snp %d edges %d || ",i,snpfrag[i].edges); for (j=0;j<snpfrag[i].edges;j++) fprintf(stdout,"%d ",snpfrag[i].elist[j]); fprintf(stdout,"\n"); getchar();
        if (snpfrag[i].edges > maxdeg) maxdeg = snpfrag[i].edges;
        avgdeg += snpfrag[i].frags;
        // edit here june 7 2012
        if (snpfrag[i].edges == 0) continue;
        csnps++;
        if (snpfrag[i].component != -1) continue; // already labeled with component
        snpfrag[i].component = i;
        for (j = 0; j < snpfrag[i].edges; j++) {
			kh_clear(32, label_node_hash);
			label_node(snpfrag, snpfrag[i].elist[j].snp, i, label_node_hash);
		}
    }
	//if (snps > 0){
	//	fprintf(stderr, "\r%.2lf%% complete\n", i * 100.0 / snps);
	//}
    /*
    fprintf_time(stderr,"FRAGMENTS=%d",fragments);
    for (i = 0; i < fragments; i++){
        fprintf_time(stderr,"i=%d",i);
        Flist[i].component = snpfrag[Flist[i].list[0].offset].component; // each fragment has a component fixed
    }*/

    *components = 0;
    int nodes_in_graph = 0;
    for (i = 0; i < snps; i++) {
        if (snpfrag[i].component == i && snpfrag[i].csize > 1) {
            (*components)++;
            nodes_in_graph += snpfrag[i].csize;
        }
        //else if (snpfrag[i].component ==i || snpfrag[i].edges ==0) singletons++;
    }
	kh_destroy(32, label_node_hash);
    //fprintf(stdout, "\nno of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
    fprintf_time(stderr, "no of non-trivial connected components %d max-Degree %d connected variants %d coverage-per-variant %f \n", *components, maxdeg, nodes_in_graph, (double) avgdeg / (double) csnps);
}
//////////////////////////////////////// edge list is only used in the two functions below //////////////////////////////

// this function updates the data structure snpfrag which links the heterozyous SNPs and the haplotype fragments
// // it generates a list of fragments (flist) that affect each SNP

void init_variant(struct fragment* Flist,int fragments,struct SNPfrags* snpfrag, int snps) 
{ 
    int i=0,j=0,k=0;
    for (i=0;i<snps;i++) snpfrag[i].frags =0;
    // calculate number of fragments covering each variant 
    for (i = 0; i < fragments; i++) {
        for (j = 0; j < Flist[i].blocks; j++) {
            for (k = 0; k < Flist[i].list[j].len; k++) snpfrag[Flist[i].list[j].offset + k].frags++;
        }
    }

    for (i = 0; i < snps; i++) {
        snpfrag[i].flist = (int*) malloc(snpfrag[i].frags*sizeof (int));
        //snpfrag[i].alist = (char*) malloc(snpfrag[i].frags*sizeof(char));
        snpfrag[i].jlist = (int*) malloc(snpfrag[i].frags *sizeof(int));
        snpfrag[i].klist = (int*) malloc(snpfrag[i].frags *sizeof(int));
    }

    for (i = 0; i < snps; i++) {
        snpfrag[i].component = -1; // this will automatically be < 0 for homozygous variants
        snpfrag[i].csize = 1;
        snpfrag[i].frags = 0; // reset to 0 since it is used and updated again (see below)
        snpfrag[i].edges = 0; // this will be zero for all homozygous variants
        snpfrag[i].tedges = 0; 
	snpfrag[i].parent = -1; 
	snpfrag[i].score = 0.0;
	snpfrag[i].post_notsw = 0;
	snpfrag[i].post_hap = 0;
	snpfrag[i].pruned_discrete_heuristic = 0;
	
    }
}

void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps) {
    int f = 0,h = 0 , j = 0, s=0, k = 0, calls = 0; //maxdeg=0,avgdeg=0;
    int prev = 0,curr=0;

    init_variant(Flist,fragments,snpfrag,snps);  // calculate number of fragments covering the variant and allocate space for lists

    // find the first fragment whose endpoint lies at snp 'i' or beyond
    for (f = 0; f < fragments; f++) {
        calls = 0;
        for (j = 0; j < Flist[f].blocks; j++) {
            for (k = 0; k < Flist[f].list[j].len; k++) {
                s = Flist[f].list[j].offset + k;           // index in snp list
                h = snpfrag[s].frags; // index into snpfrag.flist

                // save the f,j,k indices that we found this fragment spot at
                // so that from the index in snpfrag.flist alone, we can rapidly
                // get to the spot in the  global fragment* Flist.
                snpfrag[s].flist[h] = f;
                snpfrag[s].jlist[h] = j;
                snpfrag[s].klist[h] = k;
 //               snpfrag[s].alist[h] = Flist[f].list[j].hap[k];

                snpfrag[s].frags++;
                if (snpfrag[s].phase == '1') calls +=1;  // only non-homozygous variants 
            }
        }

        // calculate the number of edges per variant 
	    prev = -1; curr = -1;
            for (j = 0; j < Flist[f].blocks; j++) {
                for (k = 0; k < Flist[f].list[j].len; k++) 
		{
			if (snpfrag[Flist[f].list[j].offset+k].phase == '0') continue; 
			curr = Flist[f].list[j].offset+k;
			if (LONG_READS ==0) snpfrag[Flist[f].list[j].offset + k].edges += calls - 1;
			else // long reads, therefore only maximum of 2 edges for every node in fragment
			{
				if (prev >= 0 && curr >= 0) { snpfrag[prev].edges += 1; snpfrag[curr].edges += 1; } 
			}
			prev = curr; 
		}
            }
    }
    //for (i=0;i<snps;i++) { fprintf(stdout,"SNP %d firstfrag %d start snp %d \n",i,snpfrag[i].frags,i); }
}

