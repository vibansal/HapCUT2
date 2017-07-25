#include "bamread.h"

char INT_CIGAROP[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', 'E', 'X'};

extern int DATA_TYPE;

int QVoffset = 33;

int fetch_func(const bam1_t *b, void *data, struct alignedread* read) {
    samfile_t *fp = (samfile_t*) data;
    uint32_t *cigar = bam1_cigar(b);
    const bam1_core_t *c = &b->core;
    int i, op, ol;
    read->cigs = 0;
    read->alignedbases = 0;
    read->clipped = 0;
    read->span = 0;
    read->gapped = 0;
    read->cflag = 0;
    read->readlength = b->core.l_qseq;
    read->sequence = (char*) malloc(b->core.l_qseq + 1);
    read->quality = (char*) malloc(b->core.l_qseq + 1);

    uint8_t* sequence = bam1_seq(b);
    uint8_t* quality = bam1_qual(b);
    for (i = 0; i < b->core.l_qseq; i++) read->sequence[i] = bam_nt16_rev_table[bam1_seqi(sequence, i)];
    read->sequence[i] = '\0';
    if (quality[0] == 255) // quality string is missing, 01/29/2014, quality is set to minimum quality value specified using --minq
    {
        for (i = 0; i < b->core.l_qseq; i++) read->quality[i] = (char) (MINQ + 33);
        read->quality[i] = '\0';
    } else {
        for (i = 0; i < b->core.l_qseq; i++) read->quality[i] = (char) (quality[i] + 33);
        read->quality[i] = '\0';
    }
    //fprintf(stderr,"quality |%d| \n",quality[1]);

    read->flag = c->flag;
    read->mquality = c->qual;
    read->position = c->pos + 1;
    read->mateposition = c->mpos + 1;
    read->IS = c->isize;
    read->strand = '+';
    if ((read->flag & 16) == 16) read->strand = '-'; // fixed sept 29 2011

    read->cigarlist = (int*) malloc(sizeof (int)*c->n_cigar);
    read->cigs = c->n_cigar;
    for (i = 0; i < c->n_cigar; ++i) {
        read->cigarlist[i] = cigar[i];
        op = cigar[i]&0xf;
        ol = cigar[i] >> 4;
        if (op == BAM_CMATCH) {
            read->alignedbases += ol;
            read->span += ol;
        } else if (op == BAM_CDEL) {
            read->gapped += 1;
            read->span += ol;
        } else if (op == BAM_CINS) {
            read->alignedbases += ol;
            read->gapped += 1;
        } else if (op == BAM_CREF_SKIP) read->span += ol;
        else if (op == BAM_CSOFT_CLIP) read->clipped += ol;
        else if (op == BAM_CHARD_CLIP) {
        } else read->cflag = 1;
    }
    //      fprintf(stderr," read IS %d \n",c->isize);
    //uint8_t* barcode = NULL;
    char* barcode = NULL;
    read->barcode = NULL;
	if ( !(read->flag & 4) && DATA_TYPE == 2) // 10X reads
	{
		barcode = (char *) bam_aux_get(b,"BX");
		if (barcode && barcode != NULL){
            read->barcode = (char*)malloc(strlen(barcode)); //else read->barcode = NULL;
        	for (op=1;op<strlen(barcode);op++){
                read->barcode[op-1] = barcode[op];
            }
            read->barcode[op-1]= '\0';
        }else read->barcode = NULL;
	}

    //if (read->mquality >= 60) read->mquality = 60; // cap it at 60 april 18 2012
    read->readid = (char*) malloc(c->l_qname + 1);
    unsigned char* qs = b->data;
    for (i = 0; i < c->l_qname; i++) read->readid[i] = qs[i];
    read->readid[i] = '\0';

    if (c->tid >= 0) read->chrom = fp->header->target_name[c->tid];
    else read->chrom = NULL;
    if (c->mtid >= 0) read->matechrom = fp->header->target_name[c->mtid];
    else read->matechrom = NULL;
    read->tid = c->tid;
    read->mtid = c->mtid;
    //fprintf(stdout,"%s %s %d %d\n",read->chrom,read->matechrom,read->IS,c->mtid);
    // for MAQ bam files, mtid is not set resulting in lack of paired-end reads, may 1 2012
    return 0;
}

void free_readmemory(struct alignedread* read) {
    free(read->readid);
    free(read->sequence);
    free(read->quality);
    free(read->barcode);
    if (read->cigs > 0) free(read->cigarlist);
}
