#include "hapfragments.h"

// sort such that mate pairs are together and reverse sorted by starting position of second read in a mate-piar

int compare_fragments(const void *a, const void *b) {
    FRAGMENT* f1 = (FRAGMENT*) a;
    FRAGMENT* f2 = (FRAGMENT*) b;
    if (f1->matepos == f2->matepos) return strcmp(f1->id, f2->id);
    else return f2->matepos - f1->matepos;
}

int compare_alleles(const void *a, const void *b) {
    if (((allele*) a)->varid == ((allele*) b)->varid) return ((allele*) a)->allele - ((allele*) b)->allele;
    else return ((allele*) a)->varid - ((allele*) b)->varid;
}

int print_fragment(FRAGMENT* fragment, VARIANT* varlist, FILE* outfile) {
    if (PRINT_FRAGMENTS == 0) return 0;
    int i = 0;
    /*
       fprintf(stdout,"HAIR %s %d \t",fragment->id,fragment->variants);
       j = fragment->alist[i].varid;
       for (i=0;i<fragment->variants;i++) fprintf(stdout,"%d %s %d %c/%c %c %c \t",j,varlist[j].chrom,varlist[j].position,varlist[j].allele1,varlist[j].allele2,fragment->alist[i].allele,fragment->alist[i].qv);
       fprintf(stdout,"\n");
     */
    // varid is the index of the variant in the list of variants (this should actually  be the index in the VCF file)
    // fragment is printed using 1-based coordinate system instead of 0-based since this is encoded in HapCUT

    fragment->blocks = 1;
    for (i = 0; i < fragment->variants - 1; i++) {
        if (fragment->alist[i + 1].varid - fragment->alist[i].varid != 1) fragment->blocks++;
    }
    fprintf(outfile, "%d %s", fragment->blocks, fragment->id);

    //new format prints col 3 as data type (0 for normal, 1 for HiC) and col 4 as mate 2 index
    if (DATA_TYPE == 2){
        if(fragment->barcode != NULL)
            fprintf(outfile, " 2 %s -1", fragment->barcode);
        else
            fprintf(outfile, " 2 NULL -1");
    }else if (NEW_FORMAT)
        fprintf(outfile, " %d -1 -1", DATA_TYPE);

    //for (i=0;i<fragment->variants;i++) fprintf(stdout,"%c",fragment->alist[i].qv);
    // varid is printed with offset of 1 rather than 0 since that is encoded in the Hapcut program
    fprintf(outfile, " %d %c", fragment->alist[0].varid + 1, fragment->alist[0].allele);
    for (i = 1; i < fragment->variants; i++) {
        if (fragment->alist[i].varid - fragment->alist[i - 1].varid == 1) fprintf(outfile, "%c", fragment->alist[i].allele);
        else fprintf(outfile, " %d %c", fragment->alist[i].varid + 1, fragment->alist[i].allele);
    }
    fprintf(outfile, " ");
    for (i = 0; i < fragment->variants; i++) fprintf(outfile, "%c", fragment->alist[i].qv);
    fprintf(outfile, "\n");

    return 0;
}

// make sure they are in the correct order, i+1 could be < i

int print_matepair(FRAGMENT* f1, FRAGMENT* f2, VARIANT* varlist, FILE* outfile) {
    if (PRINT_FRAGMENTS == 0) return 0;
    int i = 0;
    f1->blocks = 1;
    for (i = 0; i < f1->variants - 1; i++) {
        if (f1->alist[i + 1].varid - f1->alist[i].varid != 1) f1->blocks++;
    }
    if (f2->alist[0].varid - f1->alist[f1->variants - 1].varid != 1) f1->blocks++;
    for (i = 0; i < f2->variants - 1; i++) {
        if (f2->alist[i + 1].varid - f2->alist[i].varid != 1) f1->blocks++;
    }

    fprintf(outfile, "%d %s_MP", f1->blocks, f1->id);
    //fprintf(outfile,"%d %s_%s_MP",f1->blocks,varlist[f1->alist[0].varid].chrom,f1->id);
    //	for (i=0;i<f1->variants;i++) fprintf(outfile,"%c",f1->alist[i].qv);
    //	for (i=0;i<f2->variants;i++) fprintf(outfile,"%c",f2->alist[i].qv);

    //new format prints col 3 as data type (0 for normal, 1 for HiC) and col 4 as mate 2 index
    if (DATA_TYPE == 2){
        if(f1->barcode != NULL)
            fprintf(outfile, " 2 %s -1", f1->barcode);
        else
            fprintf(outfile, " 2 NULL -1");
    }
    else if (NEW_FORMAT)
        fprintf(outfile, " %d %d %d", DATA_TYPE, f2->alist[0].varid+1, f1->absIS);

    // varid is printed with offset of 1 rather than 0 since that is encoded in the Hapcut program
    fprintf(outfile, " %d %c", f1->alist[0].varid + 1, f1->alist[0].allele);

    for (i = 1; i < f1->variants; i++) {
        if (f1->alist[i].varid - f1->alist[i - 1].varid == 1) fprintf(outfile, "%c", f1->alist[i].allele);
        else fprintf(outfile, " %d %c", f1->alist[i].varid + 1, f1->alist[i].allele);
    }
    if (f2->alist[0].varid - f1->alist[f1->variants - 1].varid == 1) fprintf(outfile, "%c", f2->alist[0].allele);
    else fprintf(outfile, " %d %c", f2->alist[0].varid + 1, f2->alist[0].allele);
    for (i = 1; i < f2->variants; i++) {
        if (f2->alist[i].varid - f2->alist[i - 1].varid == 1) fprintf(outfile, "%c", f2->alist[i].allele);
        else fprintf(outfile, " %d %c", f2->alist[i].varid + 1, f2->alist[i].allele);
    }
    fprintf(outfile, " ");
    for (i = 0; i < f1->variants; i++) fprintf(outfile, "%c", f1->alist[i].qv);
    for (i = 0; i < f2->variants; i++) fprintf(outfile, "%c", f2->alist[i].qv);
    fprintf(outfile, "\n");

    //	fprintf(outfile," type:");
    //	for (i=0;i<f1->variants;i++) fprintf(outfile,"%d:%d,",varlist[f1->alist[i].varid].type,varlist[f1->alist[i].varid].position);
    //	for (i=0;i<f2->variants;i++) fprintf(outfile,"%d:%d,",varlist[f2->alist[i].varid].type,varlist[f2->alist[i].varid].position);
    /*
       fprintf(outfile,"mated frag %s matepos %d vars %d \t",f1->id,f1->matepos,f1->variants);
       for (j=0;j<f1->variants;j++)
       {
       k = f1->alist[j].varid;
       fprintf(outfile,"%d %s %d %c/%c %c %c \t",k,varlist[k].chrom,varlist[k].position,varlist[k].allele1,varlist[k].allele2,f1->alist[j].allele,f1->alist[j].qv);
       }
       fprintf(outfile,"\n");
       fprintf(outfile,"mated frag %s matepos %d vars %d \t",f2->id,f2->matepos,f2->variants);
       for (j=0;j<f2->variants;j++)
       {
       k = f2->alist[j].varid;
       fprintf(stdout,"%d %s %d %c/%c %c %c \t",k,varlist[k].chrom,varlist[k].position,varlist[k].allele1,varlist[k].allele2,f2->alist[j].allele,f2->alist[j].qv);
       }
       fprintf(stdout,"\n");
     */
    return 0;
}


// sort the fragment list by 'mate-position or position of 2nd read' so that reads that are from the same DNA fragment are together
// also takes care of overlapping paired-end reads to avoid duplicates in fragments

void clean_fragmentlist(FRAGMENT* flist, int* fragments, VARIANT* varlist, int currchrom, int currpos, int prevchrom) {
    int i = 0, j = 0, k = 0, first = 0, sl = 0, bl = 0;
    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*1000);
    if (*fragments > 1) qsort(flist, *fragments, sizeof (FRAGMENT), compare_fragments);
    // sort such that mate pairs are together and reverse sorted by starting position of second read in a mate-piar
    //for (i=0;i<*fragments;i++) fprintf(stdout,"frag %s %d vars %d \n",flist[i].id,flist[i].alist[0].varid,flist[i].variants);
    if (currchrom == prevchrom) // need to ignore the top of the fragment list
    {
        first = 0;
        while (flist[first].matepos >= currpos && first < *fragments) first++;
    }
    //	fprintf(stdout,"cleaning the fragment list: current chrom %d %d first %d fragments %d\n",currchrom,currpos,first,*fragments);

    if (*fragments > 1) // bug fixed jan 13 2012, when there is only one fragment, we don't need to check if it is part of mate-pair
    {
        // serious bug fixed here: mate-pairs being examined twice April 5 2012
        // check this code for corrrectness: mate-pairs will be adjacent to each other.
        i = first;
        while (i < (*fragments) - 1) {
            if (strcmp(flist[i].id, flist[i + 1].id) == 0) // mate pair with both ends having at least one variant
            {
                //fprintf(stdout,"mate-pair %s %s %s\n",flist[i].id);
                if (flist[i].alist[flist[i].variants - 1].varid < flist[i + 1].alist[0].varid) print_matepair(&flist[i], &flist[i + 1], varlist, fragment_file);
                else if (flist[i + 1].alist[flist[i + 1].variants - 1].varid < flist[i].alist[0].varid) print_matepair(&flist[i + 1], &flist[i], varlist, fragment_file);
                else if (flist[i].variants + flist[i + 1].variants > 2) {
                    j = 0;
                    k = 0;
                    fragment.variants = 0;
                    while (j < flist[i].variants || k < flist[i + 1].variants) {
                        if (j >= flist[i].variants) {
                            fragment.alist[fragment.variants].varid = flist[i + 1].alist[k].varid;
                            fragment.alist[fragment.variants].allele = flist[i + 1].alist[k].allele;
                            fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.variants++;
                            k++;
                            continue;
                        }
                        if (k >= flist[i + 1].variants) {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            fragment.variants++;
                            j++;
                            continue;
                        }

                        if (flist[i].alist[j].varid < flist[i + 1].alist[k].varid) {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            fragment.variants++;
                            j++;
                        } else if (flist[i].alist[j].varid > flist[i + 1].alist[k].varid) {
                            fragment.alist[fragment.variants].varid = flist[i + 1].alist[k].varid;
                            fragment.alist[fragment.variants].allele = flist[i + 1].alist[k].allele;
                            fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.variants++;
                            k++;
                        } else if (flist[i].alist[j].allele == flist[i + 1].alist[k].allele) // consistent
                        {
                            fragment.alist[fragment.variants].varid = flist[i].alist[j].varid;
                            fragment.alist[fragment.variants].allele = flist[i].alist[j].allele;
                            fragment.alist[fragment.variants].qv = flist[i].alist[j].qv;
                            if (flist[i + 1].alist[k].qv > flist[i].alist[j].qv) fragment.alist[fragment.variants].qv = flist[i + 1].alist[k].qv;
                            fragment.variants++;
                            j++;
                            k++;
                        } else {
                            j++;
                            k++;
                        }
                    }
                    if (fragment.variants >= 2) {
                        sl = strlen(flist[i].id);
                        fragment.id = (char*) malloc(sl + 1);
                        for (j = 0; j < sl; j++) fragment.id[j] = flist[i].id[j];
                        fragment.id[j] = '\0';
                        if (flist[i].barcode == NULL){
                            fragment.barcode = NULL;
                        }else{
                            bl = strlen(flist[i].barcode);
                            fragment.barcode = (char*) malloc(bl + 1);
                            for (j = 0; j < bl; j++) fragment.barcode[j] = flist[i].barcode[j];
                            fragment.barcode[j] = '\0';
                        }

                        //for (j=0;j<flist[i].variants;j++) fprintf(stdout,"%d ",flist[i].alist[j].varid); fprintf(stdout,"| ");
                        //for (j=0;j<flist[i+1].variants;j++) fprintf(stdout,"%d ",flist[i+1].alist[j].varid);
                        //fprintf(stdout,"order of variants not correct %s \t",flist[i].id);
                        print_fragment(&fragment, varlist, fragment_file);
                        free(fragment.id);
                    }

                }
                else if (flist[i].variants+flist[i+1].variants ==2 && SINGLEREADS ==1)print_fragment(&flist[i],varlist,fragment_file); // added 05/31/2017 for OPE

                //else if (flist[i].variants ==1 && flist[i+1].variants >1) print_fragment(&flist[i+1],varlist);
                //else if (flist[i].variants > 1 && flist[i+1].variants ==1) print_fragment(&flist[i],varlist);
                // april 27 2012 these PE reads were being ignored until now
                i += 2;
                // what about overlapping paired-end reads.... reads..... ???? jan 13 2012,
            } else if (flist[i].variants >= 2 || SINGLEREADS == 1) {
                print_fragment(&flist[i], varlist, fragment_file);
                i++;
            } else i++;

        }
        // last read examined if it is not paired
        if (i < *fragments) {
            if (flist[i].variants >= 2 || SINGLEREADS == 1) print_fragment(&flist[i], varlist, fragment_file);
        }
    } else // only one fragment in fraglist single end
    {
        if (flist[first].variants >= 2 || SINGLEREADS == 1) print_fragment(&flist[first], varlist, fragment_file);
    }

    // free the fragments starting from first....
    if (*fragments > 0)// check added jan 13 2012
    {
        for (i = first; i<*fragments; i++) {
            free(flist[i].id);
            free(flist[i].alist);
        }
    }
    (*fragments) = first;
    free(fragment.alist);
}
