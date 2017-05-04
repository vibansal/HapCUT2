
int get_chrom_name(struct alignedread* read, HASHTABLE* ht, REFLIST* reflist) {
    int i = read->tid;
    if (reflist->ns > 0) {
        reflist->current = i;
        if (i >= reflist->ns || i < 0 || strcmp(reflist->names[i], read->chrom) != 0) {
            reflist->current = -1;
            for (i = 0; i < reflist->ns; i++) {
                if (strcmp(reflist->names[i], read->chrom) == 0) {
                    reflist->current = i;
                    break;
                }
            }
        }
    }
    return 1;
}


// find mate in readlist and also identify PCR duplicates

void find_matepair(struct alignedread** readlist, int s, int e) {
    int i = 0, j = 0;
    int prevstart = 0, prevIS = 0, duplicates = 0;
    for (i = s; i < e; i++) {
        if (readlist[i]->IS < 0) continue;
        else if (readlist[i]->IS >= 0) // insert size is positive
        {
            if (readlist[i]->position == prevstart && readlist[i]->IS == prevIS) {
                readlist[i]->flag |= 1024; //
                duplicates++;
            }
            prevstart = readlist[i]->position;
            prevIS = readlist[i]->IS;

        }
        readlist[i]->mateindex = -1;
        if (readlist[i]->IS > 0 && readlist[i]->IS < 2000) {
            // search forward in list for mate
            for (j = i + 1; j < e && readlist[j]->position <= readlist[i]->mateposition; j++) {
                if (readlist[j]->position == readlist[i]->mateposition && readlist[j]->mateposition == readlist[i]->position && strcmp(readlist[i]->readid, readlist[j]->readid) == 0) {
                    readlist[i]->mateindex = j;
                    readlist[j]->mateindex = i;
                    break;
                }
            }

        }
    }
}
/*
int estimate_readdistance_distribution(struct alignedread** readlist, int s, int e) {
    int i = 0, j = 0, k = 0;
    int ndistances = 0;
    int MAX_SIZE = 10000;
    int intra_read_dist[MAX_SIZE];
    for (i = 0; i < MAX_SIZE; i++) intra_read_dist[i] = 0;
    float intra_read_pdf[MAX_SIZE];
    for (i = 0; i < MAX_SIZE; i++) intra_read_pdf[i] = 0.0;

    for (i = s; i < e - 1; i++) {
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;
        j = i + 1;
        while (readlist[j]->IS < 0 || ((((readlist[j]->flag & 1024) == 1024)) && j < e)) j++;

        // This appears to be a bug?
        // !!!!!!!!!!!!!!!!!!!!!!!!!!
        // readlist[i]->cluster == readlist[i]->cluster always true, should be  readlist[i]->cluster == readlist[j]->cluster?
        if (j < e && readlist[i]->cluster == readlist[i]->cluster && readlist[j]->position - readlist[i]->position < MAX_SIZE && readlist[i]->mquality > 0 && readlist[j]->mquality > 0) {
            intra_read_dist[readlist[j]->position - readlist[i]->position] += 1.0;
            ndistances++;
        }
        // !!!!!!!!!!!!!!!!!!!!!!!!!
    }
    for (i = 0; i < MAX_SIZE; i++) {
        if (intra_read_dist[i] >= 10) intra_read_pdf[i] = log10((float) intra_read_dist[i] / ndistances);
        else if (i >= 5) {
            j = i - 5;
            while (j + 10 >= MAX_SIZE) j--;
            for (k = j; k < j + 10; k++) intra_read_pdf[i] += intra_read_dist[k];
            if (intra_read_dist[i] >= 10) intra_read_pdf[i] = log10(intra_read_pdf[i]) - log10(11) - log10(ndistances);
            else {
                j = i - 50;
                while (j + 100 >= MAX_SIZE) j--;
                for (k = j; k < j + 50; k++) intra_read_pdf[i] += intra_read_dist[k];
                if (intra_read_pdf[i] < 1) intra_read_pdf[i] = 0.001;
                intra_read_pdf[i] = log10(intra_read_pdf[i]) - log10(101) - log10(ndistances);
            }
            // use a window of 10-100 around 'i' to get an average value
        }
        //fprintf(stdout,"IRD %d %d %0.2f \n",i,intra_read_dist[i],intra_read_pdf[i]);
    }
    return 0;
}
*/

int generate_single_fragment(struct alignedread** readlist, int s, int e, int length, double read_density, FRAGMENT* flist, VARIANT* varlist) {
    int j = 0, i = 0, k = 0;
    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*4096);
    for (k = s; k < e; k++) {
        i = k;
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;
        if (readlist[i]->findex >= 0) {
            for (j = 0; j < flist[readlist[i]->findex].variants; j++) {
                fragment.alist[fragment.variants].varid = flist[readlist[i]->findex].alist[j].varid;
                fragment.alist[fragment.variants].allele = flist[readlist[i]->findex].alist[j].allele;
                fragment.alist[fragment.variants].qv = flist[readlist[i]->findex].alist[j].qv;
                fragment.variants++;
            }
        }
        i = readlist[i]->mateindex;
        if (i >= 0 && readlist[i]->findex >= 0) {
            for (j = 0; j < flist[readlist[i]->findex].variants; j++) {
                fragment.alist[fragment.variants].varid = flist[readlist[i]->findex].alist[j].varid;
                fragment.alist[fragment.variants].allele = flist[readlist[i]->findex].alist[j].allele;
                fragment.alist[fragment.variants].qv = flist[readlist[i]->findex].alist[j].qv;
                fragment.variants++;
            }
        }
    }

    int unique_variants = 1;
    int hets = 0;
    int counts[4];
    int qv = 0;
    qsort(fragment.alist, fragment.variants, sizeof (allele), compare_alleles);
    for (i = 0; i < fragment.variants; i++) {
        j = fragment.alist[i].varid;
        if (i > 0 && fragment.alist[i].varid != fragment.alist[i - 1].varid) unique_variants++;
        if (i > 0 && j == fragment.alist[i - 1].varid && fragment.alist[i].allele != fragment.alist[i - 1].allele) hets++;
    }
    if (hets >= 2 || hets * 3 >= unique_variants || unique_variants < 2) // fragment only has single variant or has 2 or more heterzygous variants
    {
        free(fragment.alist);
        return 0;
    }

    FRAGMENT fp;
    fp.variants = 0;
    fp.alist = (allele*) malloc(sizeof (allele) * unique_variants);

    counts[0] = counts[1] = counts[2] = counts[3] = 0;
    counts[(int) fragment.alist[0].allele - 48]++;
    counts[(int) fragment.alist[0].allele - 48 + 2] += (int) fragment.alist[0].qv - QVoffset;

    j = 0;
    for (i = 1; i <= fragment.variants; i++) {
        if (i == fragment.variants || fragment.alist[i].varid != fragment.alist[i - 1].varid) {
            // print consensus base
            if (counts[0] > counts[1] && counts[1] <= 1) {
                fp.alist[j].varid = fragment.alist[i - 1].varid;
                fp.alist[j].allele = '0';
                qv = (QVoffset + counts[2] - counts[3]);
                if (counts[2] - counts[3] >= 60) qv = 60 + QVoffset;
                fp.alist[j].qv = (char) (qv);
                if (qv - QVoffset >= MINQ) j++;
            } else if (counts[1] > counts[0] && counts[0] <= 1) {
                fp.alist[j].varid = fragment.alist[i - 1].varid;
                fp.alist[j].allele = '1';
                qv = (QVoffset + counts[3] - counts[2]);
                if (counts[3] - counts[2] >= 60) qv = 60 + QVoffset;
                fp.alist[j].qv = (char) (qv);
                if (qv - QVoffset >= MINQ) j++;
            }
            counts[0] = counts[1] = counts[2] = counts[3] = 0;
        }
        if (i < fragment.variants) {
            counts[(int) fragment.alist[i].allele - 48]++;
            counts[(int) fragment.alist[i].allele - 48 + 2] += (int) fragment.alist[i].qv - QVoffset;
        }
    }
    /*
     */
    fprintf(stdout, "fragment %d %d \n", unique_variants, j);
    fp.id = (char*) malloc(1024);
    //if (GROUPNAME != NULL) sprintf(fp.id,"%s:%s:%d_%d_%d_%0.1f",GROUPNAME,varlist[fp.alist[0].varid].chrom,readlist[s].position,readlist[e-1].position,length,read_density);
    //else sprintf(fp.id,"%s:%d_%d_%d_%0.1f",varlist[fp.alist[0].varid].chrom,readlist[s].position,readlist[e-1].position,length,read_density);
    sprintf(fp.id, "%s:%d_%d_%d_%0.1f", varlist[fp.alist[0].varid].chrom, readlist[s]->position, readlist[e - 1]->position, length, read_density);

    fp.variants = j;
    if (j >= 2) {
        fprintf(stdout, "FRAGMENT ");
        print_fragment(&fp, varlist, stdout);
        //fprintf(stderr,"fragfile %s \n",fragment_file);
        //if (fragment_file != stdout)
        print_fragment(&fp, varlist, fragment_file);
    }
    free(fp.alist);
    free(fp.id);

    for (i = 0; i < fragment.variants; i++) {
        j = fragment.alist[i].varid;
        if (i == 0 || j != fragment.alist[i - 1].varid) fprintf(stdout, "\n %d:%d %s/%s %c:%d | ", j, varlist[j].position, varlist[j].allele1, varlist[j].allele2, fragment.alist[i].allele, fragment.alist[i].qv - 33);
        else if (fragment.alist[i].allele != fragment.alist[i - 1].allele) fprintf(stdout, "%c:%d:HET | ", fragment.alist[i].allele, fragment.alist[i].qv - 33);
        else fprintf(stdout, "%c:%d | ", fragment.alist[i].allele, fragment.alist[i].qv - 33);
    }
    fprintf(stdout, "\n");
    free(fragment.alist);
    return 1;
}

int print_read(struct alignedread** readlist, int i, int prevpos, FRAGMENT* flist, VARIANT* varlist) {
    int j = 0, v = 0;
    fprintf(stdout, "dist %5d %5d %5d ", readlist[i]->position - prevpos, readlist[i]->cluster, readlist[i]->blockid);
    fprintf(stdout, "%s %s %d-%d %d %d ", readlist[i]->readid, readlist[i]->chrom, readlist[i]->position, readlist[i]->position + readlist[i]->IS, readlist[i]->IS, readlist[i]->mquality);
    for (j = 0; j < readlist[i]->cigs; j++) fprintf(stdout, "%d%c", readlist[i]->cigarlist[j] >> 4, INT_CIGAROP[readlist[i]->cigarlist[j]&0xf]);
    if (readlist[i]->findex >= 0) {
        //fprintf(stdout," vars %d ",flist[readlist[i]->findex].variants);
        for (j = 0; j < flist[readlist[i]->findex].variants; j++) {
            v = flist[readlist[i]->findex].alist[j].varid;
            fprintf(stdout, " %d:%c:%d:%s/%s:%c%c%c ", v + 1, flist[readlist[i]->findex].alist[j].allele, varlist[v].position, varlist[v].RA, varlist[v].AA, varlist[v].genotype[0], varlist[v].genotype[1], varlist[v].genotype[2]);
        }
    }

    if (readlist[i]->mateindex < 0 || readlist[readlist[i]->mateindex]->findex < 0) {
        fprintf(stdout, " \n");
        return 1;
    }

    i = readlist[i]->mateindex;
    //fprintf(stdout," || mate:%d:%d:",readlist[i]->position,readlist[i]->mquality);
    //for (j=0;j<readlist[i]->cigs;j++) fprintf(stdout,"%d%c",readlist[i]->cigarlist[j]>>4,INT_CIGAROP[readlist[i]->cigarlist[j]&0xf]);
    if (readlist[i]->findex >= 0) {
        //fprintf(stdout," vars %d ",flist[readlist[i]->findex].variants);
        for (j = 0; j < flist[readlist[i]->findex].variants; j++) {
            v = flist[readlist[i]->findex].alist[j].varid;
            fprintf(stdout, " %d:%c:%d:%s/%s:%c%c%c ", v + 1, flist[readlist[i]->findex].alist[j].allele, varlist[v].position, varlist[v].RA, varlist[v].AA, varlist[v].genotype[0], varlist[v].genotype[1], varlist[v].genotype[2]);
            fprintf(stdout, " %d:%c:%d ", v + 1, flist[readlist[i]->findex].alist[j].allele, varlist[v].position);
        }
    }
    fprintf(stdout, " \n");
    return 1;
}

// print reads within a range (s-e)

void print_reads_window(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, int if_variant_read) {
    int i = 0, prevpos = -1;
    FRAGMENT fragment;
    fragment.variants = 0;
    fragment.alist = (allele*) malloc(sizeof (allele)*1024);
    fprintf(stdout, "cluster chr%s:%d-%d length %d reads %d\n", readlist[s]->chrom, readlist[s]->position, readlist[e - 1]->position, readlist[e - 1]->position - readlist[s]->position, e - s);
    for (i = s; i < e; i++) {
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;
        if (readlist[i]->IS > 1000) {
            fprintf(stdout, "P-DEL ");
            print_read(readlist, i, prevpos, flist, varlist);
        } else if (readlist[i]->position - prevpos > 1000 && prevpos > 0) print_read(readlist, i, prevpos, flist, varlist);
        //else if (if_variant_read ==0) print_read(readlist,i,prevpos,flist,varlist);
        //else if (readlist[i]->findex >= 0 || (readlist[i]->mateindex > 0 && readlist[readlist[i]->mateindex].findex >=0)) print_read(readlist,i,prevpos,flist,varlist);
        prevpos = readlist[i]->position;
    }
    free(fragment.alist);
}

int print_clusters(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist) {
    // function called after clusters have been generated using init_clusters function
    int cluster_start = 0, cluster_end = 0, reads = 0, blocksize = 0, i = 0, j = 0, prevpos = -1;
    double mean = 0, std = 0, value = 0, obs = 0;

    cluster_start = s;
    for (i = s; i < e - 1; i++) {
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;
        j = i + 1;
        while (readlist[j]->IS < 0 || (((readlist[j]->flag & 1024) == 1024) && j < e)) j++;
        if (j < e && readlist[i]->cluster == readlist[j]->cluster) {
            cluster_end = j;
            print_read(readlist, i, prevpos, flist, varlist);
            if (readlist[j]->mquality > 0) {
                value = (readlist[j]->position - readlist[i]->position);
                mean += value;
                std += value*value;
                reads++;
                obs++;
            }
            prevpos = readlist[i]->position;
        } else if (j < e && readlist[i]->cluster != readlist[j]->cluster) {
            // print cluster information
            print_read(readlist, i, prevpos, flist, varlist);
            blocksize = readlist[cluster_end]->position - readlist[cluster_start]->position;
            fprintf(stdout, "cluster chr%s:%d-%d length %d reads %d ", readlist[cluster_start]->chrom, readlist[cluster_start]->position, readlist[cluster_end]->position, blocksize, reads);
            if (reads >= 5) {
                mean /= obs;
                std /= obs;
                std -= mean*mean;
                if (std > 0) std = sqrt(std);
                fprintf(stdout, "mean %0.2f %0.2f std %0.2f ", mean, mean * (obs), std);
            }
            if (reads >= 5) fprintf(stdout, "good \n\n");
            else if (reads == 0) fprintf(stdout, "repeat \n\n");
            else if (reads == 1) fprintf(stdout, "singleton \n\n");
            else fprintf(stdout, "weak \n\n");
            if (readlist[j]->mquality > 0) reads = 1;
            else reads = 0;
            cluster_start = j;
            cluster_end = j;
            prevpos = readlist[i]->position;
            mean = 0;
            std = 0;
            obs = 0;
        }
    }
    return 1;
}

int init_clusters(struct alignedread** readlist, int s, int e) {
    int i = 0, prevpos = -1;
    int prevtid = -1;
    int cluster_start = 0, cl = 0;

    for (i = s; i < e; i++) readlist[i]->cluster = -1; // initialized to -1
    for (i = s; i < e; i++) {
        if (readlist[i]->IS < 0 || ((readlist[i]->flag & 1024) == 1024)) continue;
        if (prevtid != readlist[i]->tid) // new chromosome
        {
            cluster_start = i;
            cl++;
            prevtid = readlist[i]->tid;
            prevpos = readlist[i]->position;
        } else if (readlist[i]->position - prevpos >= MIN_CLUSTER_DISTANCE && prevpos > 0) {
            cluster_start = i;
            cl++;
        }
        prevpos = readlist[i]->position;
        readlist[i]->cluster = cluster_start;
    }
    fprintf(stderr, "clusters %d \n", cl);
    return cl;
}
