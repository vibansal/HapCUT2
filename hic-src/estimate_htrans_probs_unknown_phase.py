# -*- coding: utf-8 -*-

import argparse
import pysam
import sys

min_mapq = 20

# parse command line arguments to benchmark.py
def parse_args():

    parser = argparse.ArgumentParser(description='use a single-chromosome Hi-C bamfile to estimate h-trans probability based on interactions with other chromosomes')

    parser.add_argument('-c', '--chromosome_length', nargs='?', type = int, help='length of chromosome', default=247249719)
    parser.add_argument('-g', '--genome_length', nargs='?', type = int, help='length of genome', default=2897310462)
    parser.add_argument('-s', '--bin_size', nargs='?', type = int, help='size in base pairs of insert size bins', default=50000)
    parser.add_argument('-b', '--bamfile', nargs='?', type = str, help='', default=None)
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='', default='estimated_htrans_probs_unknown_phase.txt')

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def main():
    
    # get parameters
    args = parse_args()    
    C_len = args.chromosome_length
    G_len = args.genome_length
    binsize = args.bin_size
    infile = args.bamfile
    outfile = args.outfile
    
    Lg_minus_Lc = G_len-C_len
    
    map_to_other_chrom = 0
         
    max_is = -1

    with pysam.AlignmentFile(infile,"rb") as file:
        for a in file:
            if abs(a.isize) > max_is:
                max_is = abs(a.isize)

    # bins of size binsize
    # bin i concerns insert size bin starting at i*binsize    
    numbins = int(max_is/binsize) + 1
    bins = [0]*numbins
    
    # names of chromosomes to consider
    csuffixes = list(range(1,23))
    csuffixes.append("X")
    valid_chroms = set(['chr{}'.format(x) for x in csuffixes])

    with pysam.AlignmentFile(infile,"rb") as file:
       
        # go through bam file
        # count number of mates within current chrom of a given insert size,
        # and number that map to other chromosomes
        for a in file:
            if a.reference_id == -1 or a.next_reference_id == -1:
                continue
            if (a.is_unmapped or a.mate_is_unmapped
                or a.is_duplicate
                or (a.next_reference_name not in valid_chroms)):
                continue
            
            # if we are on the same chrom, only count the mate pair once (using first mate)
            if a.reference_id == a.next_reference_id and a.reference_start > a.next_reference_start:
                continue

            mapq = a.mapping_quality
            mate_mapq = a.get_tag('MQ')
            isize = abs(a.isize)        
            
            if mapq < min_mapq or mate_mapq < min_mapq:
                continue
                
            if a.reference_name == a.next_reference_name:
                
                bin_ix = int(isize / binsize)
                bins[bin_ix] += 1
            else:
                map_to_other_chrom += 1 # this mate maps to another chrom
       
            mapq = a.mapping_quality
            mate_mapq = a.get_tag('MQ')
            
        # Rc is the estimated number of interactions with other chromosomes per base pair
        Rc = map_to_other_chrom / (2 * Lg_minus_Lc)
        # for filling in missing data points with previous bin
        last_prob = 0
        
        # estimate probabilities and write to output file
        with open(outfile,'w') as o:
            for i, b in enumerate(bins):
                binstart = i*binsize
                binend   = (i+1)*binsize
                binstr   = "{}-{}".format(binstart, binend)
                
                if b > 0:
                    p_htrans = Rc * binsize / b
                else:
                    p_htrans = 1
                
                if p_htrans <= 0.5 and p_htrans >= 0:
                    last_prob = p_htrans
                else:
                    p_htrans = last_prob
                
                print("{}\t{}".format(binstr, p_htrans), file=o)
    
if __name__ == '__main__':
    main()
