
/* global parameters for the segmentation */

int BLOCK_SIZE = 200; // blocks
int BF = 10; // used for making bigger blocks for background density estimation..
int BS1 = 20; // actual number of sub-blocks per BLOCK | < BS1_max 

int COMPARE_PHASE = 0; // if VCF file is phased, compare haplotype for each fragment to phased haplotypes
int VERBOSE=1;
int FIRSTPASS = 1;

int MIN_READS_PER_FRAGMENT = 8; // minimum number of reads to consider a segment as a fragment 
int MAX_FRAG_LENGTH = 250000; // max length of fragment...
int MIN_FRAG_LENGTH = 1000; 
double MIN_READ_DENSITY= 0.0005; // for a fragment, what are the units, per base pair, 5 reads per 10 kb
double max_empty_stretch = 50000;

double read_density_global=-1;// background read density -> corresponds to either singleton reads (noise) or incorrectly mapped reads 
double block_penalty_global = -7.72, length_penalty_global = -4.61; //

int MIN_CLUSTER_DISTANCE = 10000;
int MIN_SEPARATION = 3000; // max distance between two fragments to consider them as potentially part of same fragment

#define GCBINS 20
double MIN_MAPPABILITY = 0.9; 
double MIN_GC = 0.25, MAX_GC = 0.7; 
int GC_CORRECTION = 1;// use GC correction or not
int SINGLE_READS = 0;

void init_parameters(int data_type)
{
	if (FOSMIDS ==1)  // fosmid 40 kb data
	{
		//BLOCK_SIZE = 100; BS1 = 10;
		MIN_READS_PER_FRAGMENT = 8; // minimum number of reads to consider a segment as a fragment, changed to non-empty-blocks 10/24/16
		MAX_FRAG_LENGTH = 200000; // max length of fragment...
		MIN_FRAG_LENGTH = 1000; 
		MIN_READ_DENSITY= 0.0005; // for a fragment, what are the units, per base pair, 5 reads per 10 kb
		max_empty_stretch = 30000;
		block_penalty_global = -7.72, length_penalty_global = -4.61; //
		MIN_SEPARATION = 3000; // max distance between two fragments to consider them as potentially part of same fragment
		MIN_MAPPABILITY = 0.9; // making this lose increases background read-depth !!
		MIN_GC = 0.25; MAX_GC = 0.7; GC_CORRECTION = 1; 
	}
	else if (FOSMIDS ==2) // single cell MDA data 
	{
		BLOCK_SIZE = 500; BS1 = 1; 
		MIN_READS_PER_FRAGMENT = 10; // minimum number of reads to consider a segment as a fragment 
		MAX_FRAG_LENGTH = 500000; // max length of fragment...
		MIN_FRAG_LENGTH = 5000; 
		MIN_READ_DENSITY= 0.0002; // for a fragment, what are the units, per base pair, 5 reads per 10 kb
		max_empty_stretch = 50000; // disconnect fragments
		block_penalty_global = -7, length_penalty_global = -3; //
		MIN_SEPARATION = 30000; // max distance between two fragments to consider them as potentially part of same fragment
		MIN_MQ = 30; 
		MIN_MAPPABILITY = 0.5; 
		MIN_GC = 0.2, MAX_GC = 0.8; GC_CORRECTION =0;
		
	}
}

// information about reads in a block of size 'x' bp
struct BLOCK_READ
{
	int index; int cluster; 
	int firstread; int lastread; short counts; float reads; 
	int start, end; 
	float GC; // GC percentage of window // GC content of the full fragment is important |  both GC-rich fragments and AT-rich fragments are underrepresented 
	float mappability; // mappability of short reads in this window 
	uint8_t* subcounts; // count of reads in smaller blocks of 10 base pairs each 
	short adjusted_size;
	double score_partial; double W; double W_log;
	//int poslist[60];
	
	double bscore[5];  int previous[5]; double pmean[5]; double tscore; 
	int previous_cluster; 
	int reads_window; double score_window;
	double mean,variance;
	char valid,pruned;  // if block has low or high GC content, it is marked as valid=0 
	float GC_corr; float GC_corr_log;
	//int nv,pv; // next valid and previous valid block.... 
	// bed file with 100bp windows -> read count for unique and non-uniquely mapping reads 
	double bgscore;
	double curr_score; char ignore; uint8_t peak;
};

struct FOSMID // information for the genomic interval/segment corresponding to a HMW DNA fragment, could also store the HapCUT haplotype fragment as part of it...
{
        int firstblock,lastblock; int firstread,lastread; int start,end,start1,end1; // (start1,end1) is for overlapping fragments
        double mean,reads_window,score_window,ws;
	int* varlist;  int hets,unique_variants;
	int cluster; int distance_next; char ignore; 
	int fv,lv; // BL[]
	double delta; // difference in likelihood 
	char middle; // overlapping region of two fragments 
        //FRAGMENT fragment; fragment.variants =0; fragment.alist = (allele*)malloc(sizeof(allele)*4096);

};

// global variables, priors, parameters
struct global_params
{
	struct BLOCK_READ* RL; int blocks;
	int* BL; int validblocks; 
	double* length_PDF; int lbins; // fragment length PDF 
	float counts[20]; float GCbins[20]; 
	double read_density; // background read density outside fragments
	double *log_sizes; //[BLOCK_SIZE]; 
	double frag_penalty; // penalty for a single fragment log(k/B) where k is # of fragments and B is # of bins 
	int BLOCK_SIZE; 
};

// initialize the data structure and count # of reads with starting point in each block
void init_blocklist(struct BLOCK_READ* RL,int blocks,struct alignedread** readlist,int s,int e)
{
        int creads = 0; int i=0,j=0,k=0;int block=0;
        int offset =0;   // set this to BS1
        if (SINGLE_READS ==1) offset= BS1;
        for (i=0;i<blocks;i++)
        {
                RL[i].reads=0; RL[i].counts = 0; RL[i].firstread = -1; RL[i].lastread = -1; // init variables
                for (j=0;j<BS1+offset;j++) RL[i].subcounts[j] =0;
        }

        // count # of reads in each block, starting read and last read for each block   
        for (i=s;i<e;i++)
        {
                if (readlist[i]->IS < 0 ||  ((readlist[i]->flag & 1024) ==1024)) continue; // only consider first read in pair for PE reads
                if (readlist[i]->mquality < MIN_MQ) continue;  // don't use low MQ reads

                block = (readlist[i]->position/BLOCK_SIZE);
                RL[block].counts++;
                RL[block].lastread = i; RL[block].end = readlist[i]->position;
                if (RL[block].firstread < 0) { RL[block].firstread = i; RL[block].start = readlist[i]->position; }
                readlist[i]->blockid = block;
                k = (readlist[i]->position%BLOCK_SIZE)/(BLOCK_SIZE/BS1);
                if ( (readlist[i]->flag & 1) == 1 ) // PE read  
		{
			RL[block].subcounts[k]=1; // paired-end read... 
		}
                else // single end read currently not using strand information
                {
			//if ((readlist[i]->flag & 16) == 16) k = ((readlist[i]->position-200)%BLOCK_SIZE)/(BLOCK_SIZE/BS1); // add 200 to start position of read 
                        if ((readlist[i]->flag & 16) == 16) RL[block].subcounts[k+offset]=1;     else RL[block].subcounts[k]=1; 
                }
                creads++;
        }
        for (i=0;i<blocks;i++)
        {
                RL[i].reads = RL[i].counts;  RL[i].mappability = 1; RL[i].GC = 0.5; RL[i].valid = '1'; RL[i].adjusted_size = BLOCK_SIZE;
                for (j=0;j<5;j++) {RL[i].bscore[j] = -100000000; RL[i].previous[j] = -1; }
        }
}


