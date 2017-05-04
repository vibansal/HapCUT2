#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<time.h>
#include<zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "readfasta.h"
#define MAX_BUF_SIZE 4096

// initialize the reflist datastructure

REFLIST* init_reflist(char* fastafile, REFLIST* reflist) {
    int i = 0;
    reflist = (REFLIST*) malloc(sizeof (REFLIST));
    //if this statement is used, we have to return reflist pointer for memory to remain valid after function call is complete, otherwise there will be a segfault
    reflist->ns = 0;
    reflist->names = NULL;
    reflist->lengths = NULL;
    reflist->sequences = NULL;
    reflist->current = -1;
    if (read_fastaheader(fastafile, reflist) > 0) {
        reflist->sequences = calloc(reflist->ns, sizeof (char*));
        for (i = 0; i < reflist->ns; i++) {
            reflist->sequences[i] = NULL;
            //reflist->sequences[i] = calloc(reflist->lengths[i]+1,sizeof(char)); // do not allocate memory
            if (i < 4) fprintf(stderr, "chrom %s length %d \n", reflist->names[i], reflist->lengths[i]);
        }
        if (reflist->ns > 4) fprintf(stderr, ".....\n.....\nchrom %s length %d\n\n", reflist->names[reflist->ns - 1], reflist->lengths[reflist->ns - 1]);
    } else return NULL;
    return reflist;
}


// this should actually be read_fastaindex

int read_fastaheader(char* fastafile, REFLIST* reflist) // assumed to be reference.fasta.fai
{
    char buffer[MAX_BUF_SIZE];
    char fastaindexfile[1024];
    sprintf(fastaindexfile, "%s.fai", fastafile);
    FILE* fp;
    fp = fopen(fastafile, "r"); // check if fasta file is present
    if (fp == NULL) {
        fprintf(stderr, "fasta file not found, please provide a valid reference fasta file\n");
        return -1;
    } else fclose(fp);

    fp = fopen(fastaindexfile, "r");
    if (fp == NULL) {
        fprintf(stderr, "fasta index file not found, please index the reference fasta file\n");
        return -1;
    } else fprintf(stderr, "reading fasta index file %s ... ", fastaindexfile);

    reflist->ns = 0;
    while (fgets(buffer, MAX_BUF_SIZE, fp) != NULL) {
        //fprintf(stderr,"string read %s",buffer);
        reflist->ns++;
    }
    fclose(fp);

    int i = 0, j = 0, c = 0, c1 = 0, k = 0; // char length[64];
    reflist->names = (char**) malloc(sizeof (char*)*reflist->ns);
    reflist->lengths = (int*) malloc(sizeof (int)*reflist->ns);
    reflist->offsets = (uint64_t*) malloc(sizeof (uint64_t) * reflist->ns);
    reflist->sequences = (unsigned char**) malloc(sizeof (unsigned char*)*reflist->ns);
    for (i = 0; i < reflist->ns; i++) {
        reflist->sequences[i] = NULL;
        reflist->names[i] = (char*) malloc(4096);
    }
    i = 0;
    fp = fopen(fastaindexfile, "r");
    while (fgets(buffer, MAX_BUF_SIZE, fp) != NULL) {
        c = 0;
        j = c;
        c1 = -1;
        while (buffer[c] != '\t') {
            // add additional check to allow for spaces in reference name
            if (buffer[c] == ' ' && c1 < 0) c1 = c;
            c++;
        }
        if (c1 < 0) c1 = c;
        for (k = j; k < c1; k++) reflist->names[i][k - j] = buffer[k];
        reflist->names[i][k - j] = '\0';
        //fprintf(stderr,"name ..%s.. c %d c1 %d \n",reflist->names[i],c,c1);

        while (buffer[c] == ' ' || buffer[c] == '\t') c++;
        j = c;
        while (buffer[c] != ' ' && buffer[c] != '\t') c++;
        reflist->lengths[i] = 0;
        for (k = j; k < c; k++) {
            reflist->lengths[i] *= 10;
            reflist->lengths[i] += (int) buffer[k] - 48;
        }

        while (buffer[c] == ' ' || buffer[c] == '\t') c++;
        j = c; // move to next character that is not space or tab
        while (buffer[c] != ' ' && buffer[c] != '\t') c++;
        reflist->offsets[i] = 0;
        for (k = j; k < c; k++) {
            reflist->offsets[i] *= 10;
            reflist->offsets[i] += (uint64_t) buffer[k] - 48;
        }
        //fprintf(stdout,"reflist %d %d off %ld buffer %s\n",i,reflist->lengths[i],reflist->offsets[i],buffer);
        i++;
    }
    fclose(fp);
    fprintf(stderr, "fasta file %s has %d chromosomes/contigs\n\n", fastafile, reflist->ns);
    return 1;
}

// new function that reads the chromosome 'chrom' directly from fasta file using seek

int read_chromosome(REFLIST* reflist, int chrom, FILE* fp) {
    fseek(fp, reflist->offsets[chrom], SEEK_SET);
    fprintf(stderr, "reading chromosome %s offset %ld ", reflist->names[chrom], reflist->offsets[chrom]);
    if (reflist->lengths[chrom] < 1) {
        fprintf(stderr, "size of chromosome is too small,error length %d\n", reflist->lengths[chrom]);
        return -1;
    }
    reflist->sequences[chrom] = (unsigned char*) malloc(reflist->lengths[chrom] + 1);
    int j = 0, i = 0, k = 0, bases = 0;
    char c = fgetc(fp);
    //fprintf(stderr,"first char %c \n",c);
    while (bases < reflist->lengths[chrom]) {
        if (c != '\n' && c != '\t' && c != ' ') {
            reflist->sequences[chrom][j] = toupper(c);
            j++;
            // 06/22/13 what if 'c' is not 'A|C|T|G|a|c|t|g' or 'N'
            // R = A/G | Y = C/T | S = G/c | W = A/T | K = G/T | M = A/C | B = C/G/T | D = A/G/T | H = A/C/T | V = A/C/G | N
            bases++;
        }
        c = fgetc(fp);
    }
    //fprintf(stderr,"bases %d %d \n",bases,reflist->lengths[chrom]);
    //for (i=0;i<20;i++) fprintf(stderr,"%c",reflist->sequences[chrom][i]); fprintf(stderr," ");
    // now mark the targeted bases in lower case
    if (reflist->intervals <= 0 || reflist->first_interval_chrom[chrom] < 0) {
        fprintf(stderr, "\n");
        return 1;
    }
    i = reflist->first_interval_chrom[chrom];
    while (i < reflist->intervals && reflist->intervallist[i].chrom == chrom) {
        k += reflist->intervallist[i].end - reflist->intervallist[i].start;
        /*
        for (j=reflist->intervallist[i].start;j<=reflist->intervallist[i].end;j++)
        {
                if (isupper(reflist->sequences[chrom][j]) > 0)
                {
                        reflist->sequences[chrom][j] = tolower(reflist->sequences[chrom][j]); k++;
                }
        }
         */
        i++;
    }
    fprintf(stderr, "# targeted bases on chrom is %d/%d \n", k, reflist->lengths[chrom]);
    return 1;
}

int read_next_chromosome(REFLIST* reflist, int chrom, FILE* fp) {
    fprintf(stderr, "reading next chromosome %s ", reflist->names[chrom]);
    reflist->sequences[chrom] = (unsigned char*) malloc(reflist->lengths[chrom] + 1);
    int j = 0, i = 0, k = 0;
    // first character is assumed to be '>'
    char c = fgetc(fp);
    while (c != '\n') c = fgetc(fp);
    while (c != EOF && c != '>') {
        c = fgetc(fp);
        if (c != '\n' && c != '\t' && c != ' ') {
            reflist->sequences[chrom][j] = toupper(c);
            j++;
        }
    }
    // now mark the targeted bases in lower case
    if (reflist->intervals <= 0) {
        fprintf(stderr, "\n");
        return 1;
    }
    i = reflist->first_interval_chrom[chrom];
    if (i < 0) {
        fprintf(stderr, "\n");
        return 1;
    }
    while (i < reflist->intervals && reflist->intervallist[i].chrom == chrom) {
        k += reflist->intervallist[i].end - reflist->intervallist[i].start;
        i++;
    }
    fprintf(stderr, "# targeted bases on chrom is %d/%d \n", k, reflist->lengths[chrom]);
    return 1;
}

int read_fasta(char* seqfile, REFLIST* reflist) {
    clock_t t;
    kseq_t *seq;
    gzFile fp = gzopen(seqfile, "r");
    seq = kseq_init(fp);
    if (fp == NULL) {
        fprintf(stderr, "file %s not found \n", seqfile);
        return -1;
    }
    fprintf(stderr, "reading reference sequence file %s with %d sequences\n", seqfile, reflist->ns);
    t = clock();
    int i=0, j=0;
    while (kseq_read(seq) >= 0) {
        memcpy(reflist->sequences[i], seq->seq.s, seq->seq.l);
        for (j = 0; j < seq->seq.l; j++) {
            reflist->sequences[i][j] = toupper(reflist->sequences[i][j]);
        }
        i++;
    }
    gzclose(fp); fp=NULL;
    for (i = 0; i < reflist->ns; i++) {
        reflist->sequences[i][reflist->lengths[i]] = '\0';
        if (i < 10) {
            fprintf(stderr, "%s %d ", reflist->names[i], reflist->lengths[i]);
            for (j = 0; j < 30; j++) fprintf(stderr, "%c", reflist->sequences[i][j]);
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "read reference sequence file in %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    return 1;
}

int compare_intervals(const void* a, const void* b) {
    const INTERVAL* ia = (const INTERVAL*) a;
    const INTERVAL* ib = (const INTERVAL*) b;
    if (ia->chrom == ib->chrom) {
        if (ia->start == ib->start) return ia->end - ib->end;
        else return ia->start - ib->start;
    } else return ia->chrom - ib->chrom;
}

// modified on august 16 2011 to ignore lines in bed file that have no match of chromosome name to reference fasta file

int read_bedfile(char* bedfile, REFLIST* reflist) {
    char buffer[MAX_BUF_SIZE];
    int s = 0, c = 0, i = 0, e = 0, j = 0, isann = 0;
    char chrom[1024];
    int start;
    int end;
    char prevchrom[1024];
    char annotation[1024];
    strcpy(prevchrom, "-");
    int index = 0;
    reflist->intervals = 0;
    reflist->intervallist = NULL;
    if (strcmp(bedfile, "None") == 0) return -1; // no bedfile
    FILE* fp = fopen(bedfile, "r");
    if (fp == NULL) {
        fprintf(stderr, "bed file %s not found, please check the filename\n", bedfile);
        return -1;
    }
    while (fgets(buffer, MAX_BUF_SIZE, fp) != NULL) reflist->intervals++;
    if (reflist->intervals > 0) {
        reflist->intervallist = (INTERVAL*) malloc(sizeof (INTERVAL) * reflist->intervals);
        reflist->first_interval_chrom = (int*) malloc(sizeof (int)*reflist->ns);
        for (i = 0; i < reflist->ns; i++) reflist->first_interval_chrom[i] = -1;
        fclose(fp);
    } else {
        fprintf(stderr, "no intervals in bed file \n");
        fclose(fp);
        return -1;
    }

    fp = fopen(bedfile, "r");
    j = 0;
    while (fgets(buffer, MAX_BUF_SIZE, fp) != NULL) {
        c = 0;
        while (buffer[c] == ' ' || buffer[c] == '\t') c++;
        s = c;
        while (buffer[c] != ' ' && buffer[c] != '\t') c++;
        e = c;
        for (i = s; i < e; i++) chrom[i - s] = buffer[i];
        chrom[i - s] = '\0';

        while (buffer[c] == ' ' || buffer[c] == '\t') c++;
        s = c;
        while (buffer[c] != ' ' && buffer[c] != '\t') c++;
        e = c;
        start = 0;
        for (i = s; i < e; i++) start = start * 10 + (int) buffer[i] - 48;

        while (buffer[c] == ' ' || buffer[c] == '\t') c++;
        s = c;
        while (buffer[c] != ' ' && buffer[c] != '\t' && buffer[c] != '\n' && isdigit(buffer[c])) c++;
        e = c;
        end = 0;
        for (i = s; i < e; i++) end = end * 10 + (int) buffer[i] - 48;

        isann = 0;
        if (buffer[c] != '\n') {
            while (buffer[c] == ' ' || buffer[c] == '\t') c++;
            s = c;
            while (buffer[c] != ' ' && buffer[c] != '\t' && buffer[c] != '\n') c++;
            e = c;
            if (e - s > 1) // annotation
            {
                for (i = s; i < e; i++) annotation[i - s] = buffer[i];
                annotation[i - s] = '\0';
                isann = 1;
            }
        }


        if (strcmp(chrom, prevchrom) != 0) {
            index = -1;
            for (i = 0; i < reflist->ns; i++) {
                if (strcmp(chrom, reflist->names[i]) == 0) {
                    index = i;
                    break;
                }
            }
            if (index == -1) fprintf(stderr, "no match for target name, ignoring this chromosome in bedfile %s", buffer);
        }

        // start is 1-offset while sequences stored as 0-offset
        if (index != -1) {
            reflist->intervallist[j].chrom = index;
            reflist->intervallist[j].start = start;
            reflist->intervallist[j].end = end;
            if (isann == 2) {
                reflist->intervallist[j].annotation = (char*) malloc(strlen(annotation) + 1);
                strcpy(reflist->intervallist[j].annotation, annotation);
                fprintf(stdout, "|| %d %d-%d %s \n", reflist->intervallist[j].chrom, reflist->intervallist[j].start, reflist->intervallist[j].end, reflist->intervallist[j].annotation);
            } else reflist->intervallist[j].annotation = NULL;
            j++;
            //for (j=start;j<=end;j++) reflist->sequences[index][j-1] = tolower(reflist->sequences[index][j-1]);
        }
        strcpy(prevchrom, chrom);
        //fprintf(stderr,"target|%s| |%d| |%d|\n",chrom,start,end);
    }
    fclose(fp);
    reflist->intervals = j;
    fprintf(stderr, "read %d intervals from bedfile %s\n\n", reflist->intervals, bedfile);

    // sort the list of intervals by chromosome and (start,end) pairs
    qsort(reflist->intervallist, reflist->intervals, sizeof (INTERVAL), compare_intervals);
    j = 0;
    reflist->first_interval_chrom[0] = 0;
    for (i = 0; i < reflist->intervals; i++) {
        if (reflist->intervallist[i].chrom != j) {
            reflist->first_interval_chrom[reflist->intervallist[i].chrom] = i;
            j = reflist->intervallist[i].chrom;
        }
        //fprintf(stdout,"%s %d-%d \n",reflist->names[reflist->intervallist[i].chrom],reflist->intervallist[i].start,reflist->intervallist[i].end);
    }
    //for (i=0;i<reflist->ns;i++) fprintf(stdout,"chrom %s first index %d \n",reflist->names[i],reflist->first_interval_chrom[i]);
    return 1;
}
