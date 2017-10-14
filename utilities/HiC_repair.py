# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:28:16 2015

@author: peter
"""

# imports
import sys
import argparse
#import os
import pysam
import re
import subprocess as sp

MD_REGEX = re.compile("([0-9]+)([A-Z]|\^[A-Z]+)") # from SPARTA

def run_process(cmd):
    # credit to roland smith for method of retrieving stdout and stderr: http://stackoverflow.com/questions/14059558/why-is-python-no-longer-waiting-for-os-system-to-finish
    prog = sp.Popen("{} 2>&1".format(cmd), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    out, err = prog.communicate()
    print(out.decode("ISO-8859-1"))


def parse_args():

    parser = argparse.ArgumentParser(description='HiC postprocesser')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-b1', '--bamfile1', nargs='?', type = str, help='Hi-C mate 1 bamfile aligned with BWA MEM', default='data/SRR1658767_1.sorted.bam')
    parser.add_argument('-b2', '--bamfile2', nargs='?', type = str, help='Hi-C mate 2 bamfile aligned with BWA MEM', default='data/SRR1658767_2.sorted.bam')
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='output name for repaired Hi-C bam file', default='output.sam')
    parser.add_argument('-m', '--min_mapq', nargs='?', type = int, help='minimum mapping quality of snippet to paste it onto a main alignment', default=20)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

# return the distance between two aligned segments
# ASSUMES ALIGNED SEGMENTS ARE ON THE SAME CHROMOSOME
def AS_dist(a1, a2, file1, file2):

    # if not on same chromosome distance is "Infinite"
    if file1.getrname(a1.rname) != file2.getrname(a2.rname):
        return float("Inf")

    if a1.reference_start < a2.reference_start:
        first = a1
        second = a2
    else:
        first = a2
        second = a1

    return second.reference_start - first.reference_end

def fill_in_tlen(a1, a2, file1, file2):
    d1 = abs(a1.reference_start - a2.reference_start)
    d2 = abs(a1.reference_start - a2.reference_end)
    d3 = abs(a1.reference_end - a2.reference_start)
    d4 = abs(a1.reference_end - a2.reference_end)

    return max(d1,d2,d3,d4)

# decide by our own criteria whether a2 is a supplementary alignment to a1
def is_supplementary(a1, a2, cutoff):

    l1 = 0
    l2 = 0
    for op, op_length in a1.cigartuples:
        l1 += op_length

    for op, op_length in a2.cigartuples:
        l2 += op_length

    L = max(l1,l2)

    s1 = [0]*L
    s2 = [0]*L

    curr_ix = 0
    for op, op_length in a1.cigartuples:
        if op == 0:
            s1[curr_ix:(curr_ix+op_length)] = [1]*op_length
        if op != 2 and op != 3 and op != 5:
            curr_ix += op_length

    for op, op_length in a2.cigartuples:
        if op == 0:
            s2[curr_ix:(curr_ix+op_length)] = [1]*op_length
        if op != 2 and op != 3 and op != 5:
            curr_ix += op_length

    overlap = sum([x and y for x,y in zip(s1,s2)])

    return (overlap < cutoff)


def join_MD(MD1str, MD2str):

#    # an X specifies the end of the MD tag
    MD1str += 'X' # so we catch it in the regex
    MD2str += 'X' # "
    MD1 = re.findall(MD_REGEX,MD1str)
    MD2 = re.findall(MD_REGEX,MD2str)

    firstpart = MD1[:-1]
    x1 = int(MD1[-1][0])
    x2 = int(MD2[0][0]) if len(MD2) > 0 else 0
    mid_numstr = str(x1+x2)

    mid_base = MD2[0][1] if len(MD2) > 0 else 'X'
    lastpart = MD2[1:]

    newMD = firstpart + [(mid_numstr, mid_base)] + lastpart
    # flatten to string | thanks to Joel Cornett for this elegant tuple list flattening solution: http://stackoverflow.com/questions/10632839/python-transform-list-of-tuples-in-to-1-flat-list-or-1-matrix
    newMDstr = ''.join(list(sum(newMD, ())))

    return newMDstr


def join_CIGAR(CIGAR1raw, CIGAR2raw, sep):

    # filter out clips
    CIGAR1 = []
    CIGAR2 = []

    for op, oplen in CIGAR1raw:
        if op != 4 and op != 5:
            CIGAR1.append((op, oplen))

    for op, oplen in CIGAR2raw:
        if op != 4 and op != 5:
            CIGAR2.append((op, oplen))

    if sep == 0:
        if CIGAR1[-1][0] == CIGAR2[0][0]:
            firstpart = CIGAR1[:-1]
            oplen1 = CIGAR1[-1][1]
            oplen2 = CIGAR2[0][1]
            mid_oplen = oplen1 + oplen2
            mid_op = CIGAR2[0][0]
            lastpart = CIGAR2[1:]

            return firstpart + [(mid_op, mid_oplen)] + lastpart
        else:
            return CIGAR1 + CIGAR2
    else:
        return CIGAR1 + [(3,sep)] + CIGAR2

def combine_aln(a1, a2, distance):

    # figure out which alignment is first in the genome
    if a1.pos < a2.pos:
        first = a1
        second = a2
    else:
        first = a2
        second = a1

    # new position
    new_pos = first.pos

    # new mapping quality

    new_mapq = max(first.mapping_quality, second.mapping_quality)

    # cigar strings are concatenated with (distance)N inbetween to specify gap
    new_cigar = []

    # amount of overlap sequences have
    amt_to_chomp = -distance if distance < 0 else 0

    qstart1 = first.query_alignment_start
    qend1   = first.query_alignment_end
    qstart2 = second.query_alignment_start
    qend2   = second.query_alignment_end

    # concatenate sequence strings
    new_seq = first.seq[qstart1:qend1] + second.seq[(qstart2+amt_to_chomp):qend2]

    # concatenate quality score strings
    new_qual = first.qual[qstart1:qend1] + second.qual[(qstart2+amt_to_chomp):qend2]

    new_MD = join_MD(first.get_tag("MD"), second.get_tag("MD"))
    new_cigar = join_CIGAR(first.cigartuples, second.cigartuples, distance)

    cigar_len = 0
    for op, op_len in new_cigar:
        if op in {0,1,7,8}:
            cigar_len += op_len

    assert(cigar_len == len(new_seq))
    assert(len(new_seq) == len(new_qual))

    # set all values
    a1.pos   = new_pos
    a1.mapq  = new_mapq
    a1.cigartuples = new_cigar
    a1.seq   = new_seq
    a1.qual  = new_qual
    a1.set_tag('XD',new_MD) # I haven't verified the new MDs yet so I don't want them stored in MD yet
    a1.set_tag('MD',None)
    a1.set_tag('AS',None)
    a1.set_tag('XS',None)

    return a1

# main program logic
def repair_chimeras(bamfile1, bamfile2, outfile, min_mapq):
    print('Performing Hi-C repair...')

    i = 0

    # open as bamfile objects
    file1 = pysam.AlignmentFile(bamfile1,"rb");
    file2 = pysam.AlignmentFile(bamfile2,"rb");
    out = pysam.AlignmentFile(outfile,"wh", template=file1);

    # look at first aligned segments
    f1 = next(file1)
    f1_next = next(file1)
    f2 = next(file2)
    f2_next = next(file2)

    mod_count = 0
    overlap_count = 0
    total = 0

    # main loop
    while (f1 or f2):

        if not(f1 and f2):
            qname = f1.qname if f1 else f2.qname
            print("ERROR: Both alignments not present at record {}, qname {}".format(total, qname))
            continue

        # build a list of alignments for current read.
        # f1_aln[0] has primary alignment for mate1, f1_aln[1:] are secondary alignments for mate 1
        # f2_aln[0] has primary alignment for mate2, f2_aln[1:] are secondary alignments for mate 2

        f1_aln = [f1]
        f2_aln = [f2]

        while f1_next and f1.qname == f1_next.qname:
            f1_aln.append(f1_next)
            f1_next = safenext(file1)

        while f2_next and f2.qname == f2_next.qname:
            f2_aln.append(f2_next)
            f2_next = safenext(file2)

        # check that we are actually looking at mates
        assert(f1.qname == f2.qname)

        total += 1

        modded_f1 = False
        modded_f2 = False

        if not f1.is_unmapped and not f2.is_unmapped:

           # consider pasting secondary alignments to the mate
            for i in range(0,len(f1_aln)):
                f1s = f1_aln[i]
                d = AS_dist(f1s, f2, file1, file2)

                if (file1.getrname(f1s.rname) == file2.getrname(f2.rname) and d < 2000 and is_supplementary(f1,f1s,10)
                and f1s.is_reverse != f2.is_reverse and min(f2.mapq, f1.mapq) >= min_mapq):
                    if d < 0:
                        overlap_count += 1
                        continue
                    if i == 0:
                        if len(f1_aln) > 1:
                            f1 = f1_aln[1]
                            f1_aln[0],f1_aln[1] = f1_aln[1],f1_aln[0]
                        else:
                            break

                    f2 = combine_aln(f2, f1s, d)
                    mod_count += 1
                    modded_f2 = True
                    break


            for i in range(0,len(f2_aln)):
                # don't want to mess with f2 if we've modded it
                if modded_f2 == True and i == 0:
                    continue

                f2s = f2_aln[i]
                d = AS_dist(f1, f2s, file1, file2)

                if (file2.getrname(f2s.rname) == file1.getrname(f1.rname) and d < 2000 and is_supplementary(f2,f2s,10)
                and f2s.is_reverse != f1.is_reverse and min(f2.mapq, f1.mapq) >= min_mapq):
                    if d < 0:
                        overlap_count += 1
                        continue
                    # sometimes the "primary" alignment is one of our
                    # need to reassign the secondary alignment as primary
                    if i == 0:
                        if len(f2_aln) > 1:
                            f2 = f2_aln[1]
                            f2_aln[0],f2_aln[1] = f2_aln[1],f2_aln[0]
                        else:
                            break

                    f1 = combine_aln(f1, f2s, d)
                    mod_count += 1
                    modded_f1 = True
                    break

        # set flag for mate's CIGAR
        f1.set_tag('MC', f2.cigarstring)
        f2.set_tag('MC', f1.cigarstring)

        # set the XX tag to mark repaired reads
        if modded_f1:
            f1.set_tag('XX', 1)
        else:
            f1.set_tag('XX', 0)
        if modded_f2:
            f2.set_tag('XX', 1)
        else:
            f2.set_tag('XX', 0)

        # whether reads are paired end
        f1.flag = f1.flag | 1
        f2.flag = f2.flag | 1
        # whether both reads mapped or not
        if not f1.is_unmapped and not f2.is_unmapped:
            f1.flag = f1.flag | 2
            f2.flag = f2.flag | 2

            if f1.is_reverse:
                f1.flag = f1.flag | 16
                f2.flag = f2.flag | 32
            if f2.is_reverse:
                f2.flag = f2.flag | 16
                f1.flag = f1.flag | 32
            # set flag that these are paired end
            if f1.pos < f2.pos:
                f1.flag = f1.flag | 64
                f2.flag = f2.flag | 128
            else:
                f1.flag = f1.flag | 128
                f2.flag = f2.flag | 64

        else:
            f1.flag = f1.flag & 4093
            f2.flag = f2.flag & 4093

        # set flag of mate's reversal status

        # make sure f1.flag and f2.flag are not secondary, not supplementary
        f1.flag = f1.flag & 1791
        f2.flag = f2.flag & 1791

        # fill in tlen
        if not f1.is_unmapped and not f2.is_unmapped and file1.getrname(f1.rname) == file2.getrname(f2.rname):
            tlen = fill_in_tlen(f1,f2,file1,file2)
            f1.template_length = tlen
            f2.template_length = tlen

        # print to output sam/bam
        if not f1.is_unmapped:
            out.write(f1)
        if not f2.is_unmapped:
            out.write(f2)

        # move on to the next set of mates
        i+=1
        f1 = f1_next
        f1_next = safenext(file1)
        f2 = f2_next
        f2_next = safenext(file2)

    file1.close()
    file2.close()
    out.close()

    print("{} mate-pairs processed.".format(total))
    print("{} repaired.".format(mod_count))

if __name__ == '__main__':
    args = parse_args()
    main(args.bamfile1, args.bamfile2, args.outfile, args.min_mapq)
