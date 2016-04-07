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
    if len(sys.argv) < 4:
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

def main():

    args        = parse_args()
    frag_matrix = args.frag_matrix
    vcf         = args.vcf
    bin_size    = args.bin_size
    outfile     = args.outfile

    hap1 = parse_vcf_phase(vcf)


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
    MLE_IS_count = [0]*numbins     # only two-snp mate pairs
    total_IS_count = [0]*numbins   # all mate pairs for an insert size
    MLE_sum = [0]*numbins          # running sums for the numerator of MLE


    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            mate2_ix    = int(el[3])-1
            insert_size = int(el[4])

            call_list  = el[5:(5+2*num_blks)]              # extract base call part of line
            call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            alist= [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)]

            # we only care about reads with mates
            if mate2_ix == -1 or insert_size == -1:
                continue
    
            # insert size bin
            b = int(insert_size / bin_size)      
            
            total_IS_count[b] += 1
    
            if len(alist) != 2 or alist[1][0] != mate2_ix:
                continue
            
            i1 = alist[0][0]
            i2 = alist[1][0]
            h1 = hap1[i1]
            h2 = hap1[i2]
            
            if h1 == '-' or h2 == '-':
                continue
    
            a1 = alist[0][1]
            a2 = alist[1][1]
            q1 = alist[0][2]
            q2 = alist[1][2]
            
            MLE_IS_count[b] += 1        
            
            if (a1 == a2) == (h1 == h2):        # alleles match haplotype
                MLE_sum[b] += (1-q1)*q2 + (1-q2)*q1
            else:                               # alleles don't match haplotype
                MLE_sum[b] += (1-q1)*(1-q2) + q1*q2
            
    # estimate probabilities and write to output file
    with open(outfile,'w') as o:
        for i, b in enumerate(bins):
            binstart = i*bin_size
            binend   = (i+1)*bin_size
            binstr   = "{}-{}".format(binstart, binend)

            p_htrans = MLE_sum[i] / MLE_IS_count[i] if MLE_IS_count[i] > 0 else 0
            e_interactions = p_htrans * total_IS_count[i]
            
            print("{}\t{}\t{}\t{}".format(binstr, p_htrans, e_interactions, total_IS_count[i]), file=o)


if __name__ == '__main__':
    main()
