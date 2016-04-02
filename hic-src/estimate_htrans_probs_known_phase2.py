# -*- coding: utf-8 -*-

import argparse
import sys

min_mapq = 20

# parse command line arguments to benchmark.py
def parse_args():

    parser = argparse.ArgumentParser(description='use a vcf with known phase to estimate Hi-C h-trans probability based on interactions with other chromosomes')

    parser.add_argument('-f', '--frag_matrix', nargs='?', type = str, help='hapcut fragment matrix of Hi-C data. Should be new format (col 3 is data type, col 4 is mate 2 position)')
    parser.add_argument('-v', '--vcf', nargs='?', type = str, help='vcf file (corresponding to the fragment matrix) with phase data')
    parser.add_argument('-s', '--bin_size', nargs='?', type = int, help='size in base pairs of insert size bins', default=50000)
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='', default='estimated_htrans_probs_known_phase.txt')

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# given a VCF file, simply count the number of SNPs present.
def count_SNPs(vcf_file):
    count = 0
    with open(vcf_file,'r') as infile:

        for line in infile:
            if line[:1] == '#':
                continue
            if len(line.strip().split('\t')) < 5:
                continue
            count += 1

    return count

def parse_vcf_phase(vcf_file):

    num_snps = count_SNPs(vcf_file)

    hap1 = ['-']*num_snps
    hap2 = ['-']*num_snps

    with open(vcf_file, 'r') as vcf:
        i = 0
        for line in vcf:
            if line[:1] == '#':
                continue

            elements = line.strip().split('\t')
            if len(elements) < 10:
                continue

            phase_data = elements[9]

            if phase_data[0:3] == '1|0' or phase_data[0:3] == '0|1':
                hap1[i] = phase_data[0:1]
                hap2[i] = phase_data[2:3]

            i += 1

    return hap1

class fragment():
    id          = None
    datatype    = None
    mate2_ix    = None
    insert_size = None
    alleles     = None # list of (index, call, quality)


def main():

    args        = parse_args()
    frag_matrix = args.frag_matrix
    vcf         = args.vcf
    bin_size    = args.bin_size
    outfile     = args.outfile

    hap1 = parse_vcf_phase(vcf)

    flist = []

    # get the maximum insert size in fragment file
    max_is = -1

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue
            el = line.strip().split()
            insert_size = int(el[4])

            if insert_size > max_is:
                max_is = insert_size

    # bins of size binsize
    # bin i concerns insert size bin starting at i*binsize
    numbins = int(max_is/bin_size)+1
    bins = [0]*numbins

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

    for f in flist:

        # we only care about reads with mates
        if f.mate2_ix == -1 or f.insert_size == -1:
            continue

        # bins of size binsize
        # bin i concerns insert size bin starting at i*binsize
        numbins = int(max_is/binsize) + 1
        IS_count = [0]*numbins

        # only sample fragments with exactly 1 snp per mate
        alist = []
        for a in f.alleles:
            if call1 == '-' or call2 == '-'
            continue
            alist.append(a)

        if len(alist != 2 or alist[1][0] != f.mate2_ix):
            continue

        ix = int(f.insert_size / bin_size)
        post_list[ix].append(p)

    #import pdb
    #pdb.set_trace()
    # estimate probabilities and write to output file
    with open(outfile,'w') as o:
        for i, b in enumerate(bins):
            binstart = i*bin_size
            binend   = (i+1)*bin_size
            binstr   = "{}-{}".format(binstart, binend)

            plist    = post_list[i]
            p_htrans = sum(plist) / len(plist) if len(plist) > 0 else -1

            print("{}\t{}".format(binstr, p_htrans), file=o)


if __name__ == '__main__':
    main()
