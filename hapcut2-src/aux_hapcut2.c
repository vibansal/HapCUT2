
/*
auxilliary functions for hapcut2.c
*/

int detect_long_reads(struct fragment* Flist,int fragments)
{
    int i=0;
    int long_reads=0;
    float mean_snps_per_read = 0;
    if (AUTODETECT_LONGREADS){
        for (i = 0; i < fragments; i++){
            mean_snps_per_read += Flist[i].calls;
        }
        mean_snps_per_read /= fragments;
        if (mean_snps_per_read >= 3){
            long_reads = 1;
        }else{
            long_reads = 0;
        }
    }
    fprintf(stderr,"mean number of variants per read is %0.2f \n",mean_snps_per_read);
    return long_reads;
}

int read_input_files(char* fragmentfile,char* fragmentfile2,char* variantfile,DATA* data)
{
    // READ FRAGMENT MATRIX, allows for second fragment file as well. 
    data->full_fragments = get_num_fragments(fragmentfile);
    int offset = data->full_fragments;
    if (strcmp(fragmentfile2,"None") !=0) data->full_fragments += get_num_fragments(fragmentfile2);
    data->full_Flist = (struct fragment*) malloc(sizeof (struct fragment)* data->full_fragments);
    int flag = read_fragment_matrix(fragmentfile, data->full_Flist,offset,0);
    if (strcmp(fragmentfile2,"None") !=0) read_fragment_matrix(fragmentfile2, data->full_Flist,data->full_fragments-offset,offset);

    int new_fragments = 0;
    struct fragment* new_Flist;
    int i=0,j=0;

    if (MAX_IS != -1){
        // we are going to filter out some insert sizes, some memory being lost here, need to FIX
        new_fragments = 0;
        new_Flist = (struct fragment*) malloc(sizeof (struct fragment)* data->full_fragments);
        for(i = 0; i < data->full_fragments; i++){
            if (data->Flist[i].isize < MAX_IS) new_Flist[new_fragments++] = data->Flist[i];
        }
        data->full_Flist = new_Flist;
        data->full_fragments = new_fragments;
    }
    if (flag < 0) {
        fprintf_time(stderr, "unable to read fragment matrix file %s \n", fragmentfile);
        return 0;
    }
    data->snps = count_variants_vcf(variantfile);
    if (data->snps < 0) {
        fprintf_time(stderr, "unable to read variant file %s \n", variantfile);
        return 0;
    }
    else // read VCF file
    {
	   data->snpfrag = (struct SNPfrags*) malloc(sizeof (struct SNPfrags)*data->snps);
    	   read_vcffile(variantfile, data->snpfrag, data->snps);
    }
    // check consistency of fragments with variant list (if maximum # listed in fragments is > # of variants)
    for (i=0;i<data->full_fragments;i++)
    {
	for (j=0;j<data->full_Flist[i].blocks;j++) 
	{
		if (data->full_Flist[i].list[j].offset + data->full_Flist[i].list[j].len > data->snps) 
		{
			fprintf(stderr,"ERROR, maximum variant index in fragment exceeds number of variants %d \n",data->snps);
			return -1;
		}
	}
    }
    return 1;	
}

void free_memory(struct SNPfrags* snpfrag,int snps,struct BLOCK* clist,int components)
{
    int i=0;
    for (i = 0; i < snps; i++) free(snpfrag[i].elist);
    for (i = 0; i < snps; i++) free(snpfrag[i].telist);
    int component = 0;
    for (i = 0; i < snps; i++) {
        free(snpfrag[i].flist);
        //free(snpfrag[i].alist);
        free(snpfrag[i].jlist);
        free(snpfrag[i].klist);

        if (snpfrag[i].component == i && snpfrag[i].csize > 1) // root node of component
        {
            free(clist[component].slist);
            component++;
        }
    }
    for (i = 0; i < components; i++) free(clist[i].flist);
    free(snpfrag);
    free(clist);
}

void print_output_files(DATA* data,char* variantfile, char* outputfile)
{

    fprintf_time(stderr, "OUTPUTTING PRUNED HAPLOTYPE ASSEMBLY TO FILE %s\n", outputfile);
    print_contigs(data->clist,data->components,data->HAP1,data->Flist,data->fragments,data->snpfrag,outputfile);

    float N50length = calculate_N50(data->clist,data->components,data->snpfrag,data->HAP1);
    fprintf_time(stderr,"N50 haplotype length is %0.2f kilobases \n",N50length);

    char outvcffile[4096];  sprintf(outvcffile,"%s.phased.VCF",outputfile);
    if (OUTPUT_VCF ==1) {
    	fprintf_time(stderr, "OUTPUTTING PHASED VCF TO FILE %s\n", outvcffile);
	output_vcf(variantfile,data->snpfrag,data->snps,data->HAP1,data->Flist,data->fragments,outvcffile,0);
    }

    char assignfile[4096];  sprintf(assignfile,"%s.tags",outputfile);
    if (OUTPUT_HAPLOTAGS ==1) {
	fprintf_time(stderr,"OUTPUTTING read-haplotype assignments to file %s \n",assignfile);
	// even singleton fragments can be potentially assigned.. should we pass full_fragment_list
	// segfault if we use reduced fragment list... format issue
	fragment_assignments(data->full_Flist,data->full_fragments,data->snpfrag,data->HAP1,assignfile); // added 03/10/2018 to output read-haplotype assignments
    }
}

void optimization_using_maxcut(DATA* data,int max_global_iter,int maxcut_iter)
{
    fprintf_time(stderr, "starting Max-Likelihood-Cut based haplotype assembly algorithm\n");
    int iter=0; int global_iter=0; int k=0; int i=0;
    int* slist = (int*) malloc(sizeof (int)*data->snps);
    int converged_count=0;

    // compute the component-wise score for 'initHAP' haplotype
    float miscalls = 0;
    float bestscore = 0;
    for (k = 0; k < data->components; k++) {
        data->clist[k].SCORE = 0;
        data->clist[k].bestSCORE = 0;
        for (i = 0; i < data->clist[k].frags; i++) {
            update_fragscore1(data->Flist, data->clist[k].flist[i],data->HAP1);
            data->clist[k].SCORE += data->Flist[data->clist[k].flist[i]].currscore;
        }
        data->clist[k].bestSCORE = data->clist[k].SCORE;
        bestscore += data->clist[k].bestSCORE;
        miscalls += data->clist[k].SCORE;
    }

    HTRANS_MAXBINS = 0;
    if (HIC) init_HiC(data->Flist,data->fragments,HTRANS_DATA_INFILE);
    float HIC_LL_SCORE = -80;
    float OLD_HIC_LL_SCORE = -80;

    OLD_HIC_LL_SCORE = bestscore;
    for (global_iter = 0; global_iter < max_global_iter; global_iter++){ // single iteration except for HiC 
        if (VERBOSE)
            fprintf_time(stdout, "HIC ITER %d\n", global_iter);
        for (k = 0; k < data->components; k++){
            data->clist[k].iters_since_improvement = 0;
        }
        for (i=0; i<data->snps; i++){
            data->snpfrag[i].post_hap = 0;
        }
        // RUN THE MAX_CUT ALGORITHM ITERATIVELY TO IMPROVE LIKELIHOOD
        for (iter = 0; iter < maxcut_iter; iter++) {
            if (VERBOSE)
                fprintf_time(stdout, "PHASING ITER %d\n", iter);
            converged_count = 0;
            for (k = 0; k < data->components; k++){
                if(VERBOSE && iter == 0)
                    fprintf_time(stdout, "component %d length %d phased %d %d...%d\n", k, data->clist[k].length, data->clist[k].phased,data->clist[k].offset,data->clist[k].lastvar);
                if (data->clist[k].SCORE > 0)
                    converged_count += evaluate_cut_component(data->Flist,data->snpfrag,data->clist, k, slist, data->HAP1);
                else converged_count++;
            }

            if (converged_count == data->components) {
                //fprintf(stdout, "Haplotype assembly terminated early because no improvement seen in blocks after %d iterations\n", CONVERGE);
                break;
            }
        }

        // H-TRANS ESTIMATION FOR HIC
        if (max_global_iter > 1){

            // Possibly break if we're done improving
            HIC_LL_SCORE = 0;
            for (k = 0; k < data->components; k++){
                HIC_LL_SCORE += data->clist[k].bestSCORE;
            }
            if (HIC_LL_SCORE >= OLD_HIC_LL_SCORE){
                break;
            }
            OLD_HIC_LL_SCORE = HIC_LL_SCORE;

            likelihood_pruning(data->snps,data->Flist,data->snpfrag,data->HAP1, 0); // prune for only very high confidence SNPs
            // estimate the h-trans probabilities for the next round
            estimate_htrans_probs(data->Flist,data->fragments,data->HAP1,data->snpfrag,HTRANS_DATA_OUTFILE);
        }
    }
    free(slist);
}

void post_processing(DATA* data,int split_blocks) // block splitting and pruning individual variants
{
    fprintf_time(stderr, "starting to post-process phased haplotypes to further improve accuracy\n");
    // BLOCK SPLITTING
    int split_count, new_components;
    int k=0;
    new_components = data->components;

    if (split_blocks ==1){
        split_count = 0;
        for (k=0; k<data->components; k++){
            split_count += split_block(data->HAP1, data->clist,k, data->Flist, data->snpfrag, &new_components); // attempt to split block
        }
        if (split_count > 0){
            // regenerate clist if necessary
            free(data->clist);
            data->clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*new_components);
            generate_contigs(data->Flist, data->fragments, data->snpfrag, data->snps, new_components, data->clist);
        }
        data->components = new_components;
    }else if(ERROR_ANALYSIS_MODE && !HIC){
        for (k=0; k<data->components; k++){
            // run split_block but don't actually split, just get posterior probabilities
            split_block(data->HAP1, data->clist, k, data->Flist, data->snpfrag, &new_components);
        }
    }
    if (!SKIP_PRUNE){
        likelihood_pruning(data->snps,data->Flist, data->snpfrag,data->HAP1, CALL_HOMOZYGOUS);
    }
}

