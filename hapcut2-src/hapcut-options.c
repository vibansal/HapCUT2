void print_hapcut_options() {
    fprintf(stdout, "\nHapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies\n\n");
    fprintf(stdout, "USAGE : ./HAPCUT2 --fragments fragment_file --VCF variantcalls.vcf --output haplotype_output_file\n\n");
    fprintf(stderr, "Basic Options:\n");
    fprintf(stderr, "--fragments, --f <FILENAME>:        file with haplotype-informative reads generated using the extracthairs program\n");
    fprintf(stderr, "--VCF <FILENAME>:                   variant file in VCF format (use EXACT SAME file that was used for the extracthairs program)\n");
    fprintf(stderr, "--output, --o <OUTFILE> :          file to which phased haplotype segments/blocks will be output \n");
    fprintf(stderr, "--outvcf <0/1> :                    output phased variants to VCF file (<OUTFILE>.phased.vcf), default: 0 \n");
    fprintf(stderr, "--converge, --c <int>:              cut off iterations (global or maxcut) after this many iterations with no improvement. default: 5\n");
    fprintf(stderr, "--verbose, --v <0/1>:               verbose mode: print extra information to stdout and stderr. default: 0\n");

    fprintf(stderr, "\nRead Technology Options:\n");
    fprintf(stderr, "--hic <0/1> :                       increases accuracy on Hi-C data; models h-trans errors directly from the data. default: 0\n");
    fprintf(stderr, "--hic_htrans_file, --hf <FILENAME>  optional tab-delimited input file where second column specifies h-trans error probabilities for insert size bins 0-50Kb, 50Kb-100Kb, etc.\n");
    fprintf(stderr, "--qv_offset, --qo <33/48/64> :      quality value offset for base quality scores, default: 33 (use same value as for extracthairs)\n");
    fprintf(stderr, "--long_reads, --lr <0/1> :          reduces memory when phasing long read data with many SNPs per read. default: automatic.\n");

    fprintf(stderr, "\nHaplotype Post-Processing Options:\n");
    fprintf(stderr, "--threshold, --t <float>:           PHRED SCALED threshold for pruning low-confidence SNPs (range 0-100, larger values prune more.). default: 6.98\n");
    fprintf(stderr, "--skip_prune, --sp <0/1>:           skip default likelihood pruning step (prune SNPs after the fact using column 11 of the output). default: 0\n");
    //fprintf(stderr, "--split_blocks, --sb <0/1>:         split blocks using simple likelihood score to reduce switch errors. default: 0\n");
    //fprintf(stderr, "--split_threshold, --st <float>:    PHRED SCALED threshold for splitting blocks (range 0-100, larger values split more). default: 0\n");
    fprintf(stderr, "--call_homozygous, --ch <0/1>:      call positions as homozygous if they appear to be false heterozygotes. default: 0\n");
    fprintf(stderr, "--discrete_pruning, --dp <0/1>:     use discrete heuristic to prune SNPs. default: 0\n");
    fprintf(stderr, "--error_analysis_mode, --ea <0/1>:  compute switch confidence scores and print to haplotype file but don't split blocks or prune. default: 0\n");

    fprintf(stderr, "\nAdvanced Options:\n");
    fprintf(stderr, "--new_format, --nf <0/1>:           use new Hi-C fragment matrix file format (but don't do h-trans error modeling). default: 0\n");
    fprintf(stderr, "--max_iter, --mi <int> :            maximum number of global iterations. Preferable to tweak --converge option instead. default: 10000\n");
    fprintf(stderr, "--maxcut_iter, --mc <int> :         maximum number of max-likelihood-cut iterations. Preferable to tweak --converge option instead. default: 10000\n");
    fprintf(stderr, "--htrans_read_lowbound, --hrl <int> with --hic on, h-trans probability estimation will require this many matepairs per window. default: 500\n");
    fprintf(stderr, "--htrans_max_window, --hmw <int>    with --hic on, the insert-size window for h-trans probability estimation will not expand larger than this many basepairs. default: 4000000\n");

    fprintf(stderr, "\n\nHi-C-specific Notes:\n");
    fprintf(stderr, "  (1) When running extractHAIRS, must use --hic 1 option to create a fragment matrix in the new Hi-C format.\n");
    fprintf(stderr, "  (2) When running HapCUT2, use --hic 1 if h-trans probabilities are unknown. Use --hic_htrans_file if they are known\n");
    fprintf(stderr, "  (3) Using --hic_htrans_file is faster than --hic and may yield better results at low read coverage (<30x).\n");
    fprintf(stderr, "\n");

}
