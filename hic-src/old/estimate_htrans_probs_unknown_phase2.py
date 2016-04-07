# -*- coding: utf-8 -*-

import argparse
import pysam
import sys
from collections import defaultdict

min_mapq = 20
emax = 1000 # estimate max: max insert size for estimating number of htrans per bin

# parse command line arguments to benchmark.py
def parse_args():

    parser = argparse.ArgumentParser(description='use a single-chromosome Hi-C bamfile to estimate h-trans probability based on interactions with other chromosomes')

    parser.add_argument('-s', '--bin_size', nargs='?', type = int, help='size in base pairs of insert size bins', default=50000)
    parser.add_argument('-f', '--frag_mat', nargs='?', type = str, help='', default=None)
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='', default='estimated_htrans_probs_unknown_phase.txt')

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

class fragment():
    id          = None
    datatype    = None
    mate2_ix    = None
    insert_size = None
    alleles     = None # list of (index, call, quality)

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

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            f = fragment()
            el = line.strip().split()

            num_blks      = int(el[0])
            f.id          = el[1]
            f.datatype    = el[2]
            f.mate2_ix    = int(el[3])-1
            f.insert_size = int(el[4])

            alist  = el[5:(5+2*num_blks)]              # extract base call part of line
            alist  = zip(*[iter(alist)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            alist  = [(int(a)-1, b) for a,b in alist]  # convert index to 0-based integer
            alist2 = []

            for ix, blk in alist:
                curr_ix = ix
                for a in blk:
                    alist2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            f.alleles = [(a,b,c) for ((a,b),c) in zip(alist2,qlist)]

            flist.append(f)

    number_in_ebin   = 0
    consistency_dict = defaultdict(int)

    for f in flist:
        if f.insert_size > emax:
            continue

        # only sample fragments with exactly 1 snp per mate
        alist = []
        for a in f.alleles:
            if call1 == '-' or call2 == '-'
            continue
            alist.append(a)

        if len(alist != 2):
            continue

        number_in_ebin += 1

        for (ix1,call1,qual1), (ix2,call2,qual2) in itertools.combinations(f.alleles,2):


            # consistency_dict[(ix1,ix2,0)] has the count of consistent fragments
            consistency_dict[(ix1,ix2,0)] =
            consistency_dict[(ix1,ix2,1)] = 0

    num_htrans_ebin = number_in_ebin * p_htrans_ebin
    num_htrans = num_htrans_ebin * (binsize / number_in_ebin)


    # estimate probabilities and write to output file
    with open(outfile,'w') as o:
        for i, b in enumerate(bins):
            binstart = i*binsize
            binend   = (i+1)*binsize
            binstr   = "{}-{}".format(binstart, binend)

            if b > 0:
                p_htrans = num_htrans / b
            else:
                p_htrans = 1

            #if p_htrans <= 0.5 and p_htrans >= 0:
            #    last_prob = p_htrans
            #else:
            #    p_htrans = last_prob

            print("{}\t{}".format(binstr, p_htrans), file=o)

if __name__ == '__main__':
    main()
