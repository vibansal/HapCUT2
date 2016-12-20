# -*- coding: utf-8 -*-

import argparse
import sys
import fileIO

verbose = False
min_mapq = 20
min_reads = 500
bin_cap   = 2000000

# parse command line arguments to benchmark.py
def parse_args():

    parser = argparse.ArgumentParser(description='use a vcf with known phase to estimate Hi-C h-trans probability based on interactions with other chromosomes')

    parser.add_argument('-f', '--frag_matrix', nargs='?', type = str, help='hapcut fragment matrix of Hi-C data. Should be new format (col 3 is data type, col 4 is mate 2 position)')
    parser.add_argument('-v', '--vcf', nargs='?', type = str, help='vcf file (corresponding to the fragment matrix) with phase data')
    parser.add_argument('-hb', '--hapblock', type = str, help='haplotype block file to use for phase instead of VCF file')
    parser.add_argument('-s', '--bin_size', nargs='?', type = int, help='size in base pairs of insert size bins', default=50000)
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='', default='estimated_htrans_probs_known_phase.txt')

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def estimate_htrans_probs(frag_matrix, vcf, outfile, hapblock=None, bin_size=50000):

    if vcf != None and hapblock !=None:
        print("Specify VCF or hapblock but not both")
        exit(1)
    elif vcf != None:
        hap1 = fileIO.parse_vcf_phase(vcf,use_SNP_index=True)
    elif hapblock != None:
        hap1 = fileIO.parse_hapblock_file(hapblock,use_SNP_index=True)

    hap_dict = dict()
    for blocknum, block in enumerate(hap1):
        #print(blocknum)
        for pos, a1, a2 in block:
            hap_dict[pos] = (blocknum, a1)

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

            # we only care about reads with mates
            if mate2_ix == -1 or insert_size == -1:
                continue

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

            # insert size bin
            b = int(insert_size / bin_size)

            total_IS_count[b] += 1

            if len(alist) != 2 or alist[1][0] != mate2_ix:
                continue

            i1 = alist[0][0]
            i2 = alist[1][0]

            if i1 not in hap_dict or i2 not in hap_dict:
                continue

            blocknum1, h1 = hap_dict[i1]
            blocknum2, h2 = hap_dict[i2]
            a1 = alist[0][1]
            a2 = alist[1][1]
            q1 = alist[0][2]
            q2 = alist[1][2]

            if (blocknum1 != blocknum2
              or h1 == '-' or h2 == '-'
              or a1 == '-' or a2 == '-'):
                continue

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

            #p_htrans = MLE_sum[i] / MLE_IS_count[i] if MLE_IS_count[i] > 0 else 0

            j = i+1
            k = i-1
            new_bin_size = bin_size
            adj_MLE_sum   = MLE_sum[i]
            adj_MLE_count = MLE_IS_count[i]
            while True:

                if (j < numbins):
                    new_bin_size  += bin_size
                    adj_MLE_count += MLE_IS_count[j]
                    adj_MLE_sum   += MLE_sum[j]
                    j += 1
                if (k >= 0):
                    new_bin_size  += bin_size
                    adj_MLE_count += MLE_IS_count[k]
                    adj_MLE_sum   += MLE_sum[k]
                    k -= 1

                if adj_MLE_count > min_reads or new_bin_size > bin_cap or (k<0 and j>=numbins):
                    break

            p_htrans = adj_MLE_sum / adj_MLE_count if adj_MLE_count > 0 else 0

            e_interactions = p_htrans * total_IS_count[i]

            if verbose:
                print("{}\t{}\t{}\t{}\t{}\t{}".format(binstr, p_htrans, e_interactions, MLE_IS_count[i], adj_MLE_count, total_IS_count[i]), file=o)
            else:
                print("{}\t{}".format(binstr, p_htrans), file=o)

if __name__ == '__main__':
    args        = parse_args()
    frag_matrix = args.frag_matrix
    vcf         = args.vcf
    hapblock    = args.hapblock
    bin_size    = args.bin_size
    outfile     = args.outfile
    estimate_htrans_probs(frag_matrix, vcf, outfile, hapblock, bin_size)
