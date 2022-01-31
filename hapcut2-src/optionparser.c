/* specification of the command line arguments for HapCUT2 and I/O of these
 */

// Printing related
int VERBOSE = 0;
int PRINT_FRAGMENT_SCORES = 0; // output the MEC/switch error score of erroneous reads/fragments to a file for evaluation
int OUTPUT_HAPLOTAGS =0; // read-haplotype assignments
int OUTPUT_VCF =1;     // output phased VCF or not 

// Quality-score related parameters
int QVoffset = 33; // for fragment quality values
int MINQ = 6; // additional base quality filter in hapcut added april 18 2012
int MAXQ = 40; // maximum quality score, if quality score in input fragment is > MAXQ, it will be capped

// Number of iterations
int MAXITER = 10000;     // maximum number of global iterations
int MAXCUT_ITER = 10000; // maximum number of iterations for max-cut algorithm, if this is proportional to 'N' -> complexity is 'N^2', added march 13 2013
int CONVERGE = 5; // stop iterations on a given block if exceed this many iterations since improvement

// Post-processing related variables
float THRESHOLD = 0.8;
float SPLIT_THRESHOLD = 0.8;
float HOMOZYGOUS_PRIOR = -80; // in log form. assumed to be really unlikely
int DISCRETE_PRUNING =0;    // not used 
int CALL_HOMOZYGOUS = 0;  
int SPLIT_BLOCKS = 0;
int ERROR_ANALYSIS_MODE = 0;  
int SKIP_PRUNE = 0;
int GENOTYPING = 0;  // if set to 1, do genotyping of variants
int SNVS_BEFORE_INDELS = 0; // ??? 

int AUTODETECT_LONGREADS = 1;
int LONG_READS = 0; // if this variable is 1, the data contains long read data
int PLOIDY =2; // diploid

// HiC-related global variables
int HIC = 0;
int MAX_HIC_EM_ITER = 1;
int NEW_FRAGFILE_FORMAT = 0;
int HTRANS_BINSIZE = 5000;
int HTRANS_MAXBINS = 10000; // this value will be overwritten at startup
int HTRANS_READ_LOWBOUND = 500;
int HTRANS_MAX_WINDOW = 4000000; // maximum window size for h-trans estimation
char HTRANS_DATA_INFILE[10000];
char HTRANS_DATA_OUTFILE[10000];
int MAX_IS = -1;

void print_hapcut_options();

int parse_arguments(int argc,char* argv[],char* fragfile,char* fragfile2,char* VCFfile,char* hapfile)
{
    // input arguments are initial fragment file, variant file with variant information and alleles for each variant
    // number of iterations total, when to output the solution, file to output solution .....
    int i = 0;
    int flag = 0;

    if (argc % 2 != 1){
        fprintf(stderr, "\nERROR: Invalid number of arguments specified.\n");
        exit(1);
    }

    for (i = 1; i < argc; i += 2) {
        if (argc < 6) break;

        // BASIC OPTIONS
        if (strcmp(argv[i], "--fragments") == 0 || strcmp(argv[i], "--f") == 0) {
            strcpy(fragfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--fragments2") == 0 ) {
            strcpy(fragfile2, argv[i + 1]);
        } else if (strcmp(argv[i], "--VCF") == 0 || strcmp(argv[i], "--vcf") == 0) {
            strcpy(VCFfile, argv[i + 1]);
            flag++;
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "--out") == 0|| strcmp(argv[i], "--o") == 0) {
            strcpy(hapfile, argv[i + 1]);
            flag++;
        }else if ((strcmp(argv[i], "--converge") == 0) || (strcmp(argv[i], "--c") == 0)) {
            CONVERGE = atoi(argv[i + 1]);
        }else if ((strcmp(argv[i], "--rh") == 0) || (strcmp(argv[i], "--tags") == 0)) { // read-haplotype assignments
            check_input_0_or_1(argv[i + 1]);
            OUTPUT_HAPLOTAGS = atoi(argv[i + 1]);
        }else if ((strcmp(argv[i], "--outvcf") == 0) ) { // output VCF or not, default is 1
            OUTPUT_VCF = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "--v") == 0){
            check_input_0_or_1(argv[i + 1]);
            VERBOSE = atoi(argv[i + 1]);
        }
        // READ-TECHNOLOGY OPTIONS
        else if (strcmp(argv[i], "--HiC") == 0 || strcmp(argv[i], "--hic") == 0){
            check_input_0_or_1(argv[i + 1]);
            if (atoi(argv[i + 1])){
                MAX_HIC_EM_ITER = 100; //atoi(argv[i + 1]);
                NEW_FRAGFILE_FORMAT = 1;
                HIC = 1;
            }
        }else if (strcmp(argv[i], "--long_reads") == 0 || strcmp(argv[i], "--lr") == 0){
            check_input_0_or_1(argv[i + 1]);
            LONG_READS = atoi(argv[i + 1]);
            AUTODETECT_LONGREADS = 0;
        }else if (strcmp(argv[i], "--QV_offset") == 0 || strcmp(argv[i], "--qv_offset") == 0 || strcmp(argv[i], "--qo") == 0){
            QVoffset = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--hic_htrans_file") == 0 || strcmp(argv[i], "--hf") == 0){
            NEW_FRAGFILE_FORMAT = 1;
            strcpy(HTRANS_DATA_INFILE, argv[i + 1]);
            HIC = 1;
        }
        // HAPLOTYPE POST-PROCESSING OPTIONS
        else if (strcmp(argv[i], "--threshold") == 0 || strcmp(argv[i], "--t") == 0){
            THRESHOLD = 1.0 - unphred(atof(argv[i + 1]));
        //}else if (strcmp(argv[i], "--split_blocks") == 0 || strcmp(argv[i], "--sb") == 0){
            //SPLIT_BLOCKS = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--split_threshold") == 0 || strcmp(argv[i], "--st") == 0){
            SPLIT_THRESHOLD = 1.0 - unphred(atof(argv[i + 1]));
        }else if (strcmp(argv[i], "--call_homozygous") == 0 || strcmp(argv[i], "--ch") == 0){
            check_input_0_or_1(argv[i + 1]);
            CALL_HOMOZYGOUS = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--error_analysis_mode") == 0 || strcmp(argv[i], "--ea") == 0){
            check_input_0_or_1(argv[i + 1]);
            ERROR_ANALYSIS_MODE = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--SNVs_before_indels") == 0 || strcmp(argv[i], "--si") == 0){
            check_input_0_or_1(argv[i + 1]);
            SNVS_BEFORE_INDELS = atoi(argv[i + 1]);
        }
        // ADVANCED OPTIONS
        else if (strcmp(argv[i], "--nf") == 0 || strcmp(argv[i], "--new_format") == 0){
            check_input_0_or_1(argv[i + 1]);
            NEW_FRAGFILE_FORMAT = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--max_iter") == 0 || strcmp(argv[i], "--mi") == 0){
            MAXITER = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--maxcut_iter") == 0 || strcmp(argv[i], "--mc") == 0) {
            MAXCUT_ITER = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--htrans_read_lowbound") == 0 || strcmp(argv[i], "--hrl") == 0){
            HTRANS_READ_LOWBOUND = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--htrans_max_window") == 0 || strcmp(argv[i], "--hmw") == 0){
            HTRANS_MAX_WINDOW = atoi(argv[i + 1]);
        }
        // HIDDEN OPTIONS
        else if (strcmp(argv[i], "--htrans_data_outfile") == 0 || strcmp(argv[i], "--ohf") == 0){
            strcpy(HTRANS_DATA_OUTFILE, argv[i + 1]);
        }else if (strcmp(argv[i], "--printscores") == 0 || strcmp(argv[i], "--scores") == 0){
            PRINT_FRAGMENT_SCORES = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--mbq") == 0){
            MINQ = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--maxbq") == 0){
            MAXQ = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--skip_prune") == 0 || strcmp(argv[i], "--sp") == 0){
            check_input_0_or_1(argv[i + 1]);
            SKIP_PRUNE = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--genotyping") == 0 || strcmp(argv[i], "--geno") == 0){
            check_input_0_or_1(argv[i + 1]);
            GENOTYPING = atoi(argv[i + 1]);
        }else if (strcmp(argv[i], "--max_IS") == 0 || strcmp(argv[i], "--mi") == 0){
            MAX_IS = atoi(argv[i + 1]);
        }else{
            fprintf(stderr, "\nERROR: Invalid Option \"%s\" specified.\n",argv[i]);
            exit(1);
        }
    }
    if (ERROR_ANALYSIS_MODE && HIC){
        fprintf_time(stderr,"WARNING: Switch error quality scores are not intended for use with Hi-C data. Scores will be left blank.\n");
    }

    if (flag != 3) // three essential arguments are not supplied
    {
        print_hapcut_options();
        return 0;
    }

	fprintf(stderr, "\n\n");
    fprintf_time(stderr, "input fragment file: %s\n", fragfile);
    fprintf_time(stderr, "input variantfile (VCF format):%s\n", VCFfile);
    fprintf_time(stderr, "haplotypes will be output to file: %s\n", hapfile);
    fprintf_time(stderr, "solution convergence cutoff: %d\n", CONVERGE);
    //fprintf_time(stderr, "QVoffset: %d\n\n", QVoffset);
    return 0;
}


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
    //fprintf(stderr, "--error_analysis_mode, --ea <0/1>:  compute switch confidence scores and print to haplotype file but don't split blocks or prune. default: 0\n");

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
