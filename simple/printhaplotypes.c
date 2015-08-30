
#include "common.h"


// THIS FUNCTION PRINTS THE CURRENT HAPLOTYPE ASSEMBLY in a new file block by block 

void print_hapcut_options() {
    fprintf(stdout, "\nHAPCUT: haplotype assembly using cut computations, last updated 09/11/13 \n\n");
    fprintf(stdout, "USAGE : ./HAPCUT --fragments fragment_file --VCF variantcalls.vcf --output haplotype_output_file > hapcut.log\n\n");
    fprintf(stderr, "=============== PROGRAM OPTIONS ======================================== \n\n");
    fprintf(stderr, "--fragments <FILENAME>: file with haplotype-informative reads generated using the extracthairs program\n");
    //fprintf(stderr,"--variants : variant file in hapCUT format (same file as used for extracthairs)\n");
    fprintf(stderr, "--VCF <FILENAME>: variant file in VCF format ((same file as used for running the extracthairs program)\n");
    fprintf(stderr, "--output <FILENAME> : file to which phased haplotype segments/blocks will be output | the program will write best current haplotypes to this file after every 10 iterations\n");
    fprintf(stderr, "--maxiter <int> : maximum number of global iterations for HAPCUT, default is 100\n");
    fprintf(stderr, "--qvoffset <33/48/64> : quality value offset for base quality scores, default is 33 (use same value as for extracthairs)\n");
    fprintf(stderr, "--maxcutiter <int> : maximum number of iterations to find max cut for each haplotype block in a given iteration, default value is 100 | if this is set to -1, the number of iterations = N/10 where N is the number of nodes in the block\n");
    fprintf(stderr, "--longreads <0/1> : set to 1 for phasing long read data if program uses too much memory, default is 0\n");
    //fprintf(stderr,"--fosmids <0/1> : set to 1 for phasing fosmid pooled sequencing data, default is 0\n");
    fprintf(stderr, "--sf <0/1> : scoring function for comparing reads against haplotypes: default is 0 (MEC score), set to 1 (switch error rate) for fosmid data\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "NOTES:\n");
    fprintf(stderr, " 1. Make sure to use the same VCF file for running extracthairs and HapCUT\"\n");
    fprintf(stderr, " 2. For phasing fosmid data or synthetic long read data, use the options \"--fosmids 1 --sf 1\"\n");
    fprintf(stderr, " 3. For phasing Hi-C data, use a large value for the maxIS option in extracthairs and set the paramter maxcutiter = 100 \n ");
    fprintf(stderr, "\n");

}

void print_hapfile(struct component* clist, int components, int snps, char* h1, struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, char* fname, char* outfile, char* pruned) {
    // print a new file containing one block phasing and the corresponding fragments 
    int i = 0, j=0, t = 0, numblocks=0, firstSNP = 1;
    char c, c1, c2;
    FILE* fp;
    fp = fopen(outfile, "w");
    
    for (i = 0; i < components; i++) {
        firstSNP = 1;
        for (j = 0; j < clist[i].size; j++){
            t = clist[i].slist[j];
            if (h1[t] == '-') continue;
            
            if (firstSNP || pruned[t] == 4){
                if (numblocks > 0) fprintf(fp, "******** \n");
                fprintf(fp, "BLOCK: offset: %d\n", t+1);
                numblocks++;
                firstSNP = 0;
            }

            if (h1[t] == '0') c = '1';
            else if (h1[t] == '1') c = '0';


            // print this line to keep consistency with old format
            // if SNP was pruned then print '-'s
            if (pruned[t] == 1)
                fprintf(fp, "%d\t-\t-", t + 1);
            else if (pruned[t] == 2)
                fprintf(fp, "%d\t0\t0", t + 1);
            else if (pruned[t] == 3)
                fprintf(fp, "%d\t1\t1", t + 1);
            else{
                if (snpfrag[t].genotypes[0] == '2' || snpfrag[t].genotypes[2] == '2') {
                    //fprintf(stderr,"tri-allelic variant %d:%d %s \n",t,snpfrag[t].position,snpfrag[t].genotypes);
                    if (h1[t] == '0') {
                        c1 = snpfrag[t].genotypes[0];
                        c2 = snpfrag[t].genotypes[2];
                    } else if (h1[t] == '1') {
                        c2 = snpfrag[t].genotypes[0];
                        c1 = snpfrag[t].genotypes[2];
                    }
                    fprintf(fp, "%d\t%c\t%c", t + 1, c1, c2); // two alleles that are phased in VCF like format
                } else {
                    fprintf(fp, "%d\t%c\t%c", t + 1, h1[t], c);
                }
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
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

