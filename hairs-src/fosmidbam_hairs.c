#include<stdint.h>
int MIN_CLUSTER_DISTANCE = 10000;
#include "print_clusters.c"
#include "sam.h"
int BLOCK_SIZE = 100; // 100 bp blocks
int BF = 10;
#define BS1 20

// code for parsing fosmid pool bam files

// ./extractFOSMID --bam fosmid-data/SRR799544_GGACTCCTAAGGAGTA.rmdup.bam.chr1 --VCF fosmid-data/NA12878.1kg.hets.chr1.liftoverhg19.vcf --fosmids 1 --ref fosmid-data/chr1.fa > na

// block of reads in a window of size 'x' bp

struct BLOCK_READ {
    int index;
    int cluster;
    int firstread;
    int lastread;
    short reads;
    int start, end;
    float GC; // GC percentage of window // GC content of the full fragment is important |  both GC-rich fragments and AT-rich fragments are underrepresented 
    float mappability; // mappability of short reads in this window 
    short subcounts[BS1];

    double bscore;
    int previous;
    int previous_cluster;
    int reads_window;
    double score_window;
    double mean, variance;
    double bscore1;
    int previous1;
    double score_window1;
    // bed file with 100bp windows -> read count for unique and non-uniquely mapping reads 
};


// block based analysis gives us probability that a block has 1 or more reads -> background poisson distribution for reads...

// should be called on reads from single chromosome....

int cluster_reads(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, REFLIST* reflist) {
    int reads = 0, i = 0, j = 0, k = 0, first = 0, last = 0, empty = 0;
    double ws = 0, ws_corrected = 0;
    double blockscore = 0, prior_size = 0, mean = 0, variance = 0;

    double mean_block_size = 40000;
    double std_block_size = 8000;
    double pconst = log(std_block_size) + 0.5 * log(2 * 3.14159);
    double block_open_penalty = log(1e-6);
    //double exp_mean = 1.0/500000;
    double gap_length_penalty = 0;

    int blocks = readlist[e - 1]->position / BLOCK_SIZE + 1;
    struct BLOCK_READ* RL = calloc(sizeof (struct BLOCK_READ), blocks);
    for (i = 0; i < blocks; i++) {
        RL[i].reads = 0;
        RL[i].firstread = -1;
        for (j = 0; j < 10; j++) RL[i].subcounts[j] = 0;
    }
    int creads = 0;
    int empty_blocks = 0, block = 0, nonempty_blocks = 0;


    for (i = s; i < e; i++) {
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;

        block = (readlist[i]->position / BLOCK_SIZE);
        RL[block].reads++;
        RL[block].lastread = i;
        RL[block].end = readlist[i]->position;
        if (RL[block].firstread < 0) {
            RL[block].firstread = i;
            RL[block].start = readlist[i]->position;
        }
        readlist[i]->blockid = block;
        k = (readlist[i]->position % BLOCK_SIZE) / (BLOCK_SIZE / BS1);
        RL[block].subcounts[k]++;
        creads++;
    }
    for (i = 0; i < blocks; i += BF) {
        if (i + BF >= blocks) continue; // extra check at boundaey
        empty = 0;
        for (j = 0; j < BF; j++) empty += RL[i + j].reads;
        if (empty == 0) empty_blocks++;
        else nonempty_blocks++;
        for (j = 0; j < BF; j++) {
            if (RL[i + j].reads < 4) continue;
            empty = 0;
            for (k = 0; k < BS1; k++) if (RL[i + j].subcounts[k] > 0) empty++;
            //for (k=0;k<BS1;k++) fprintf(stdout,"%d:",RL[i+j].subcounts[k]);
            //fprintf(stdout," block %d reads %d %d first %d last %d \n",i+j,RL[i+j].reads,empty,RL[i+j].firstread,RL[i+j].lastread);
            RL[i + j].reads = empty;
        }
    }
    double read_density = (1.0 - (float) empty_blocks / (empty_blocks + nonempty_blocks)) / (BF * BLOCK_SIZE);
    double score0, lambda;
    //read_density = 0.00009;//*BLOCK_SIZE;
    int eb = 0;

    fprintf(stderr, "blocks %d %d %0.8f\n", blocks, empty_blocks, read_density);
    // use bigger block size of 1000bp to calculate background density =  1/10000 bases 


    fprintf(stderr, "clustering reads using dynamic programming algorithm reads %d...%d %d %d %f \n", readlist[s]->position, readlist[e - 1]->position, e - s + 1, creads, read_density);
    fprintf(stdout, "\nclustering reads using dynamic programming algorithm reads %d %d %f \n\n", e - s + 1, creads, read_density);

    // initialize all values	
    for (i = 0; i < blocks; i++) {
        RL[i].previous = -1;
        RL[i].bscore = -10000000;
        RL[i].GC = 0;
        RL[i].previous1 = -1;
        RL[i].bscore1 = -10000000;
        if (reflist->current < 0) continue;
        for (k = 0; k < BLOCK_SIZE; k++) {
            if (reflist->sequences[reflist->current][i * BLOCK_SIZE + k] == 'G' || reflist->sequences[reflist->current][i * BLOCK_SIZE + k] == 'C') RL[i].GC += 1;
        }
        RL[i].GC /= BLOCK_SIZE;
    }

    RL[0].bscore = RL[0].reads * log(read_density) - read_density*BLOCK_SIZE;
    for (i = 1; i < blocks; i++) {
        score0 = RL[i].reads * log(read_density) - read_density*BLOCK_SIZE;
        RL[i].bscore = RL[i - 1].bscore + score0;
        RL[i].previous = i - 1;
        RL[i].bscore1 = RL[i - 1].bscore1 + score0;
        RL[i].previous1 = i - 1;
        RL[i].previous_cluster = RL[i - 1].previous_cluster;
        RL[i].score_window = score0;

        j = i - 1;
        reads = RL[i].reads;
        eb = 0;
        ws = 0;
        ws_corrected = 0;
        variance = RL[i].reads * RL[i].reads;
        while ((i - j) < (int) (200000 / BLOCK_SIZE) && j >= 0 && eb * BLOCK_SIZE < 10000) {
            if (RL[j].reads > 0) eb = 0;
            else eb++;
            reads += RL[j].reads;
            variance += RL[j].reads * RL[j].reads;
            if (RL[j].GC < 0.999) ws_corrected += BLOCK_SIZE;
            //if (RL[j].GC < 0.67) ws_corrected += BLOCK_SIZE;
            ws += BLOCK_SIZE;

            if (reads >= 5 && ws_corrected >= 1000) {
                mean = (double) reads / (i - j);
                lambda = (double) reads / ws_corrected;
                blockscore = reads * (log(lambda) - 1);
                //gap_length_penalty = log(1.0-exp((RL[j].previous_cluster-j-1)*exp_mean*BLOCK_SIZE)); 
                //fprintf(stderr,"penalty %d %f \n",RL[j].previous_cluster-j,gap_length_penalty);
                if (j == 0) score0 = blockscore + block_open_penalty + prior_size + gap_length_penalty;
                else score0 = RL[j - 1].bscore + blockscore + block_open_penalty + prior_size + gap_length_penalty;

                if (score0 > RL[i].bscore) {
                    RL[i].bscore = score0;
                    RL[i].previous = (j - 1);
                    RL[i].reads_window = reads;
                    RL[i].score_window = blockscore;
                    RL[i].previous_cluster = i;
                    RL[i].variance = variance / (i - j) - mean*mean;
                    RL[i].mean = mean;
                }

                prior_size = -1 * (ws - mean_block_size)*(ws - mean_block_size) / (2 * std_block_size * std_block_size) - pconst;
                if (score0 + prior_size > RL[i].bscore1) {
                    RL[i].bscore1 = score0 + prior_size;
                    RL[i].previous1 = (j - 1);
                    RL[i].score_window1 = blockscore;
                }
            }
            j--;
        }
        j = RL[i].previous;
        //if (RL[i].reads > 0) fprintf(stdout,"score1 for %d..%d ws %d reads %d score %f | %d:%d:%d %f \n",i,j,ws,reads,RL[i].bscore,RL[i].start,RL[i].end,RL[i].reads,RL[i].bscore);
        if (i % 500000 == 0) fprintf(stderr, "done for block %d %d\n", i, blocks);
    }

    int flag = 0;
    int block0_start = -1;
    i = blocks - 1;
    while (i >= 0) {
        j = RL[i].previous;
        ws = (i - j) * BLOCK_SIZE;
        if (j < i - 1) {
            if (flag * BLOCK_SIZE < 3000) fprintf(stdout, "potential merge... ");
            if (reflist->current >= 0 && flag * BLOCK_SIZE < 3000) {
                for (k = 0; k < flag; k++) fprintf(stdout, "%d:%0.2f | ", (k + block0_start) * BLOCK_SIZE, RL[block0_start - k].GC);
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "empty stretch %d %d-%d\n\n", flag*BLOCK_SIZE, (block0_start - flag) * BLOCK_SIZE, block0_start * BLOCK_SIZE);
            fprintf(stdout, "\n==block %d...%d length %0.0f score %0.2f %0.2f 1 | ", j*BLOCK_SIZE, (i) * BLOCK_SIZE, ws, RL[i].bscore, RL[i].score_window);
            //fprintf(stdout,"%d %0.2f | ",BLOCK_SIZE*(i-RL[i].previous1),RL[i].bscore1);
            fprintf(stdout, "%d....%d reads %d density %0.2f var %0.4f:%0.4f %0.1f ==\n", j, i, RL[i].reads_window, (float) ws / RL[i].reads_window, RL[i].mean, RL[i].variance, RL[i].mean * RL[i].mean / (RL[i].variance - RL[i].mean));

            first = RL[j + 1].firstread;
            k = 1;
            while (first < 0) {
                first = RL[j + 1 + k].firstread;
                k++;
            }
            last = RL[i].lastread;
            k = 1;
            while (last < 0) {
                last = RL[i - k].lastread;
                k++;
            }
            generate_single_fragment(readlist, first, last, ws, (float) ws / RL[i].reads_window, flist, varlist);
            // bug here last can be -1 sometimes 
            print_reads_window(readlist, first, last, flist, varlist, 1);

            flag = 0;
            block0_start = -1;
        } else {
            if (block0_start < 0) block0_start = j;
            // empty block... 
            flag++;
        }
        i = j;
    }
    return 1;

    free(RL);
}

void process_chunk(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, REFLIST* reflist) {
    find_matepair(readlist, s, e);
    fprintf(stderr, "e %d s %d \n", e, s);

    //int cl = init_clusters(readlist,s,e); // cluster using a maximum intra-cluster distance value

    //estimate_readdistance_distribution(readlist,s,e,cluster); // estimate distances between start positions of adjacent reads within same cluster 
    // cluster size distribution from data | probability of a read being a singleton read
    cluster_reads(readlist, s, e, flist, varlist, reflist);

    fprintf(stdout, "\n\n");
    //print_clusters(readlist,s,e,flist,varlist);  // print clusters
}

// want to process all reads for a single chromosome since calculation of read-density requires data... 
// readlist is too big,not cleaned up after every chromosome... done 02/11/2014 
// extract haplotype informative reads from sorted bam file //
// need to discard reads that are marked as duplicates using flag //

int parse_bamfile_fosmid(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist) {
    fprintf(stderr, "reading sorted bamfile %s for fosmid pool\n", bamfile);
    int reads = 0;
    int MAX_READS = 5000000; // 10 million for now
    //struct alignedread* read = (struct alignedread*)malloc(sizeof(struct alignedread));
    struct alignedread** readlist = calloc(MAX_READS, sizeof (struct alignedread*));
    for (reads = 0; reads < MAX_READS; reads++) readlist[reads] = calloc(1, sizeof (struct alignedread));
    struct alignedread* read_pt;

    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*10000);
    FRAGMENT* flist = (FRAGMENT*) malloc(sizeof (FRAGMENT) * MAX_READS / 5);
    int fragments = 0;

    int chrom = 0;
    int r = 0, i = 0;
    int prevchrom = -1;
    int prevtid = -1; //int prevposition = -1; // position of previous read in sorted bam file
    int lastread = 0;

    samfile_t *fp;
    if ((fp = samopen(bamfile, "rb", 0)) == 0) {
        fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
        return -1;
    }
    bam1_t *b = bam_init1();

    while (samread(fp, b) >= 0) {
        //readlist[r] = calloc(1,sizeof(struct alignedread));
        fetch_func(b, fp, readlist[r]);
        if ((readlist[r]->flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) // unmapped reads, PCR/optical dups are ignored
        {
            free_readmemory(readlist[r]);
            continue;
        }
        // find the chromosome in reflist that matches read->chrom if the previous chromosome is different from current chromosome
        // if too many reads, break off when distance between adjacent reads is large 
        if (readlist[r]->tid != prevtid || r >= MAX_READS - 1) {
            if (r >= MAX_READS - 1) fprintf(stderr, "limit on max reads %d exceeded.. need to clean up buffer \n", MAX_READS);
            if (prevtid >= 0) {
                fprintf(stderr, "reads in buffer %d \n", r);
                process_chunk(readlist, lastread, r, flist, varlist, reflist);
                // free up reads in list and move up index 
                for (i = lastread; i < r; i++) free_readmemory(readlist[i]);
                read_pt = readlist[0];
                readlist[0] = readlist[r];
                readlist[r] = read_pt;
                r = 0;
                for (i = 0; i < fragments; i++) {
                    free(flist[i].alist);
                    free(flist[i].id);
                }
                fprintf(stderr, "free memory for reads from chrom %d cleaning up of fragment list %d\n", prevtid, fragments);
                fragments = 0;
            }

            chrom = getindex(ht, readlist[r]->chrom);
            get_chrom_name(readlist[r], ht, reflist); // reflist->current has chromosome index, -1 if no reflist 
            lastread = r;
        } else chrom = prevchrom;

        fragment.variants = 0; //fragment.id = readlist[r]->readid;
        if (chrom >= 0) extract_variants_read(readlist[r], ht, chromvars, varlist, 1, &fragment, chrom, reflist);
        if (fragment.variants > 0) {
            add_fragment(flist, &fragment, readlist[r], fragments);
            readlist[r]->findex = fragments++;
        } else readlist[r]->findex = -1;

        reads += 1;
        if (reads % 2000000 == 0) fprintf(stderr, "processed %d reads, useful fragments \n", reads);
        prevchrom = chrom;
        prevtid = readlist[r]->tid;
        //if (readlist[r]->IS == 0) prevposition = readlist[r]->position; 
        //else if (readlist[r]->IS > 0) prevposition = readlist[r]->position + readlist[r]->IS; // outer end of fragment r1....r2
        r++;
    }
    process_chunk(readlist, lastread, r, flist, varlist, reflist);
    for (i = lastread; i < r; i++) free_readmemory(readlist[i]);
    for (reads = 0; reads < MAX_READS; reads++) free(readlist[reads]);
    free(readlist);
    bam_destroy1(b);
    return 0;
}

