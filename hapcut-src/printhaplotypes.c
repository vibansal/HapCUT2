
// THIS FUNCTION PRINTS THE CURRENT HAPLOTYPE ASSEMBLY in a new file block by block 

void print_hapcut_options() {
    fprintf(stdout, "\nHAPCUT: haplotype assembly using cut computations, last updated 07/8/2015 \n\n");
    fprintf(stdout, "USAGE : ./HAPCUT [options]* -v <vcf_file> -f <fragment_file>\n\n");
    fprintf(stderr, "=============== REQUIRED ARGUMENTS ======================================== \n\n");
    fprintf(stderr, "-f, --fragments   <FILENAME>: file with haplotype-informative reads generated using the extracthairs program, or simulated with happysim.py\n\n");
    fprintf(stderr, "-v,  --vcf        <FILENAME>: variant file in VCF format (use same file as used for running the extracthairs program).\n");
    fprintf(stderr, "=============== OPTIONAL ARGUMENTS ======================================== \n\n");
    fprintf(stderr, "-o,  --output     <FILENAME>: print detailed haplotype output file (defaults to hapcut2.out, assuming VCF was provided). The program will write best current haplotypes to this file after every 10 iterations\n");
    fprintf(stderr, "-s   --simple     <FILENAME>: print simple haplotype output file for benchmarking (this option allows you to skip VCF file (-v))\n");
    fprintf(stderr, "-m,  --maxiter    <int>     : maximum number of global iterations for HAPCUT (default 100)\n");
    fprintf(stderr, "-q,  --qvoffset   <33/48/64>: quality value offset for base quality scores (default 33) (use same value as for extracthairs)\n");
    fprintf(stderr, "-mc, --maxcutiter <int>     : maximum number of iterations to find max cut for each haplotype block in a given iteration, default value is 100. If this is set to -1, the number of iterations = N/10 where N is the number of nodes in the block\n");
    fprintf(stderr, "-l,  --longreads            : set this flag for long read data if program uses too much memory (default False)\n");
    fprintf(stderr, "-fo, --fosmid               : set this flag for phasing fosmid pooled sequencing data (default False)\n");
    fprintf(stderr, "-sw   --switch               : use switch error rate (for fosmid data) instead of MEC (default False)\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "NOTES:\n");
    fprintf(stderr, " 1. Always use the same VCF file for running extracthairs and HapCUT\"\n");
    fprintf(stderr, " 2. For phasing fosmid data or synthetic long read data, use the options \"-fo -s\"\n");
    fprintf(stderr, " 3. For phasing Hi-C data, use a large value for the maxIS option in extracthairs and use option \" -mc 100\" \n ");
    fprintf(stderr, "\n");
}

int print_hapfile(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fname, int score, char* outfile) {
    // print a new file containing one block phasing and the corresponding fragments 
    if (is_none(outfile)){
        return 0;
    }
    int i = 0, t = 0, k = 0, span = 0;
    char c, c1, c2;
    FILE* fp;

    //char fn[200]; sprintf(fn,"%s-%d.phase",fname,score);

    fp = fopen(outfile, "w");
    float delta = 0;

    for (i = 0; i < blocks; i++) {
        span = snpfrag[clist[i].lastvar].position - snpfrag[clist[i].offset].position;
        fprintf(fp, "BLOCK: offset: %d len: %d phased: %d ", clist[i].offset + 1, clist[i].length, clist[i].phased);
        fprintf(fp, "SPAN: %d MECscore %2.2f fragments %d\n", span, clist[i].bestMEC, clist[i].frags);
        for (k = 0; k < clist[i].phased; k++) {
            t = clist[i].slist[k];
            //fprintf(fp,"frags %d | ",snpfrag[t].frags); 
            //if (clist[i].haplotype[k] =='-') fprintf(fp,"%s_%s_%d_%s_%s_%s\t%10c\t%10c\n",snpfrag[t].id,snpfrag[t].chromosome,snpfrag[t].position,snpfrag[t].allele0,snpfrag[t].allele1,snpfrag[t].genotypes,'-','-'); 
            // changed code here to use the component
            //if (snpfrag[clist[i].offset+k].component == snpfrag[clist[i].offset].component)	fprintf(fp,"%d\t%c\t%c\t%s\t%d\t%s\t%s\t%s\n",t+1,'-','-',snpfrag[t].chromosome,snpfrag[t].position,snpfrag[t].allele0,snpfrag[t].allele1,snpfrag[t].genotypes); 
            //if (snpfrag[clist[i].offset+k].component == snpfrag[clist[i].offset].component)	
            {
                if (h1[t] == '0') c = '1';
                else if (h1[t] == '1') c = '0';
                delta = snpfrag[t].L00 - snpfrag[t].L01;
                if (snpfrag[t].L11 - snpfrag[t].L01 > delta) delta = snpfrag[t].L11 - snpfrag[t].L01;
                // print this line to keep consistency with old format 

                if (snpfrag[t].genotypes[0] == '2' || snpfrag[t].genotypes[2] == '2') {
                    //fprintf(stderr,"tri-allelic variant %d:%d %s \n",t,snpfrag[t].position,snpfrag[t].genotypes);
                    if (h1[t] == '0') {
                        c1 = snpfrag[t].genotypes[0];
                        c2 = snpfrag[t].genotypes[2];
                    } else if (h1[t] == '1') {
                        c2 = snpfrag[t].genotypes[0];
                        c1 = snpfrag[t].genotypes[2];
                    }
                    fprintf(fp, "%d\t%c\t%c\t", t + 1, c1, c2); // two alleles that are phased in VCF like format
                } else {
                    fprintf(fp, "%d\t%c\t%c\t", t + 1, h1[t], c);
                }
                fprintf(fp, "%s\t%d\t%s\t%s\t%s\t%d,%d:%0.1f,%0.1f,%0.1f:%0.1f:%0.1f", snpfrag[t].chromosome, snpfrag[t].position, snpfrag[t].allele0, snpfrag[t].allele1, snpfrag[t].genotypes, snpfrag[t].R0, snpfrag[t].R1, snpfrag[t].L00, snpfrag[t].L01, snpfrag[t].L11, delta, snpfrag[t].rMEC);
                if (delta >= 3 && snpfrag[t].rMEC >= 2) fprintf(fp, ":FV");

                // print genotype read counts and likelihoods
                if (snpfrag[t].G00 < 0 || snpfrag[t].G01 < 0 || snpfrag[t].G11 < 0) fprintf(fp, "\t%d,%d:%0.1f,%0.1f,%0.1f\n", snpfrag[t].A0, snpfrag[t].A1, snpfrag[t].G00, snpfrag[t].G01, snpfrag[t].G11);
                else fprintf(fp, "\n");
                //fprintf(fp,"%s_%s_%d_%s_%s_%s\t%10c\t%10c\tDELTA:%f\n",snpfrag[t].id,snpfrag[t].chromosome,snpfrag[t].position,snpfrag[t].allele0,snpfrag[t].allele1,snpfrag[t].genotypes,h1[t],c,snpfrag[t].deltaLL); 
                //fprintf(fp,"%s\t%c\t%c \n",snpfrag[t].id,h1[t],c); 
            }
        }
        if (i < blocks - 1) fprintf(fp, "******** \n");
    }
    fclose(fp);

    return 0;
}

int print_simplefile(char* h1, char* h2, char* simplefile){
    
    if (is_none(simplefile)){
        return 0;
    }
    // print haplotypes 1 and 2 to files for simple benchmarking
    // credit to dasblinkenlight on stackoverflow for full buffering soln: http://stackoverflow.com/questions/11558540/c-why-is-a-fprintfstdout-so-slow
    char buf[10000];
    setvbuf(stdout, buf, _IOFBF, sizeof(buf));
    FILE* fp;
    fp = fopen(simplefile, "w");
    fprintf(fp,"h1\n");
    fprintf(fp, h1);
    fprintf(fp,"\nh2\n");
    fprintf(fp, h2); 
    fclose(fp);
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// important NOTE: to get from SNP i to its component in 'clist', we have variable snpfrag[i].bcomp, feb 1 2012

void print_haplotypes_vcf(struct BLOCK* clist, int blocks, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, char* outfile) {
    // print the haplotypes in VCF like format 
    FILE* fp;
    fp = fopen(outfile, "w");
    int i = 0, k = 0, span = 0, component;
    char c1, c2;
    for (i = 0; i < snps; i++) {
        fprintf(fp, "%s\t%d\t%s\t%s\t%s\t", snpfrag[i].chromosome, snpfrag[i].position, snpfrag[i].id, snpfrag[i].allele0, snpfrag[i].allele1);
        fprintf(fp, "%s\t", snpfrag[i].genotypes);
        component = snpfrag[i].component;
        if (component < 0) fprintf(fp, "./.\n");
        else if (snpfrag[component].csize < 2) fprintf(fp, "./.\n");
        else {
            k = snpfrag[i].bcomp;
            if (clist[k].offset == i) // main node of component
                //if (component ==i) // main node of component changed feb 14 2013
            {
                if (h1[i] == '0') {
                    c1 = '0';
                    c2 = '1';
                } else {
                    c1 = '1';
                    c2 = '0';
                }
                span = snpfrag[clist[k].lastvar].position - snpfrag[clist[k].offset].position;
                fprintf(fp, "SP:%c|%c:%d ", c1, c2, clist[k].offset + 1);
                fprintf(fp, "BLOCKlength: %d phased %d SPAN: %d MECscore %2.2f fragments %d\n", clist[k].length, clist[k].phased, span, clist[k].bestMEC, clist[k].frags);
            } else {
                if (h1[i] == '0') {
                    c1 = '0';
                    c2 = '1';
                } else {
                    c1 = '1';
                    c2 = '0';
                }
                fprintf(fp, "SP:%c|%c:%d\n", c1, c2, clist[k].offset + 1);
            }
        }
    }
}

