#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#Created on Jul 6 2017
#
#@author: Peter Edge, pedge@eng.ucsd.edu


from collections import defaultdict
import sys
import pickle
import pysam
import argparse
import os

barcode_tag = 'BX'

###############################################################################
###############################################################################
# these functions pulled from "getMolecules" script:
# https://github.com/RCollins13/10XWGS/blob/master/getMolecules.py

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Estimate original molecule sizes and coordinates from 10X linked-read WGS barcodes
"""
import argparse
from collections import defaultdict, Counter, namedtuple
from subprocess import call
import pysam

min_mapq = 20

def get_gemcode_regions(ibam, dist):
    """
    Estimates molecule coordinates from colinear reads with overlapping 10X barcodes

    Parameters
    ----------
    ibam : pysam.AlignmentFile
        Input 10X bam
    dist : int
        Partitioning distance (bp) for separating reads with overlapping
        barcodes into independent fragments

    Yields
    ------
    region : namedtuple
        chr, start, end, barcode, readcount
    """

    #Create namedtuples for storing read coordinate and molecule info
    coords = namedtuple('coords', ['chr', 'pos', 'end'])
    molecule = namedtuple('molecule', ['chr', 'start', 'end', 'barcode', 'readcount'])

    #Create defaultdict for storing gemcode tuples
    #Key is gemcodes, value is coords namedtuple
    gemcodes = defaultdict(list)

    #Set current_chr as reference contig of first read
    current_chr = None
    #print("bam file",ibam, file=sys.stderr); #current_chr = None

    #Iterate over reads in bamfile
    for read in ibam:

        if (read.mapq < min_mapq or read.is_unmapped or read.is_duplicate or read.is_secondary or
           read.is_qcfail):
            continue

        assert(read.reference_name  != None)
        assert(read.reference_start != None)
        assert(read.reference_end   != None)

        #Save 10X barcode as gem
        if not read.has_tag('BX'):
            continue

        gem = read.get_tag('BX')
	#print read.reference_name,read.reference_start,read.reference_end

        #If the read is from a new contig/chromosome, write out all molecules held in memory
        #then add read to emptied dictionary
        if read.reference_name != current_chr and current_chr is not None:
            for barcode in gemcodes:
                yield molecule(gemcodes[barcode][0].chr, min([pos for chr,
                               pos, end in gemcodes[barcode]]),
                               max([end for chr, pos, end in gemcodes[barcode]]),
                               barcode, len(gemcodes[barcode]))

            gemcodes = defaultdict(list)

            gemcodes[gem].append(coords(read.reference_name, read.reference_start, read.reference_end))  # added read.reference_end -- pedge

        #If barcode has been seen previously and new read is colinear but beyond
        #dist, yield old barcode as interval before adding new read to list
        elif gem in gemcodes and read.reference_start - gemcodes[gem][-1].pos > dist:

            yield molecule(gemcodes[gem][0].chr,
                           min([pos for chr, pos, end in gemcodes[gem]]),
                           max([end for chr, pos, end in gemcodes[gem]]),
                           gem, len(gemcodes[gem]))

            gemcodes[gem] = [coords(read.reference_name, read.reference_start, read.reference_end)]

        else:
            #Else just add read to preexisting dictionary
            gemcodes[gem].append(coords(read.reference_name, read.reference_start, read.reference_end))

        #Save read contig as current_chr
        current_chr = read.reference_name

    #Write out all remaining molecules at end of bam
    for barcode in gemcodes:
        yield molecule(gemcodes[barcode][0].chr, min([pos for chr,
                       pos, end in gemcodes[barcode]]),
                       max([end for chr, pos, end in gemcodes[barcode]]),
                       barcode, len(gemcodes[barcode]))

def get_molecules(bam,ref=None,dist=20000):

    #Open outfile
    bamf = pysam.AlignmentFile(bam,"rb");
    if ref == None:
        ibam = bamf
    else:
        ibam = bamf.fetch(reference=ref)

    #Get gemcode regions
    for bed in get_gemcode_regions(ibam, dist):
        #Turn molecule object into string
        yield (bed.chr, bed.start, bed.end, bed.barcode)

# end of getMolecules code
###############################################################################
###############################################################################


# a class representing a HapCUT2 format haplotype fragment
class fragment:

    def __init__(self, seq, name, barcode, dtype = 0):
        self.seq = seq                         # list of (snp index, genomic index, allele call, quality score) tuples
        self.name = name                       # fragment ID / name
        self.barcode = barcode
        self.dtype = dtype
        self.used = False

    def __str__(self):
        fragstr = ''
        num_pairs = 0
        prev_snp_ix = -2
        qual = ' '
        for snp_ix, genome_ix, allele, q_char in self.seq:

            diff = snp_ix - prev_snp_ix

            if diff == 1:
                fragstr += allele
            else:
                num_pairs += 1
                barcode = self.barcode if self.barcode != None else 'NULL'
                fragstr += ' {} {}'.format(snp_ix+1, allele)

            prev_snp_ix = snp_ix
            qual += q_char

        fragstr += qual

        barcode = self.barcode if self.barcode != None else "NULL"
        if self.dtype == 0:
            prefix = '{} {} 0 -1 -1'.format(num_pairs,self.name)
        elif self.dtype == 2:
            prefix = '{} {} 2 {} -1'.format(num_pairs,self.name, self.barcode)
        fragstr = prefix + fragstr
        return fragstr

# read in a HapCUT2 format fragment file into a list of fragment objects
def read_fragment_matrix(frag_matrix, vcf_file, chrom_filter=None):

    snp_ix = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            chrom = el[0]
            genomic_pos = int(el[1])-1

            if chrom_filter == None or chrom_filter == chrom:
                vcf_dict[snp_ix] = genomic_pos

            snp_ix += 1

    flist = []

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            name = el[1]

            dtype = int(el[2])

            if not dtype == 2:
                print("Input to LinkFragments should be unlinked 10X fragments (datatype 2), obtained by running extractHAIRS with --10X 1 option. Current datatype is {}.".format(dtype))
                exit(1)

            barcode = el[3]

            if barcode == 'NULL':
                barcode = None

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

            skip = False
            for a,b in call_list2:
                if a not in vcf_dict:
                    skip = True
            if skip:
                continue
            #qlist = [(ord(q) - 33) for q in qlist]

            alist= [(a,vcf_dict[a],b,c) for ((a,b),c) in zip(call_list2,qlist)]

            frag = fragment(alist,name,barcode,dtype=2)
            flist.append(frag)

    flist.sort(key=lambda x: x.seq[0][0])

    return flist

# print out a list of fragment objects to a HapCUT2 format fragment file
def write_fragment_matrix(flist,outfile, single_SNP_frags=False):
    lines = []

    for f in flist:
        if not single_SNP_frags and len(f.seq) < 2:
            continue

        if len(f.seq) == 0:
            continue

        firstpos = f.seq[0][0]
        lines.append((firstpos, str(f)))

    lines.sort()

    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)

def parse_bedfile(input_file):

    boundaries = []
    with open(input_file,'r') as inf:
        for line in inf:

            if len(line) < 3:
                continue

            el = line.strip().split('\t')

            chrom = el[0]
            start = int(el[1])
            stop  = int(el[2])
            barcode = el[3]
            yield (chrom, start, stop, barcode)

def link_fragments(hairs_file, vcf_file, bam_file, outfile, dist, single_SNP_frags):

    lines = []

    if not (os.path.isfile(bam_file+'.bai')):
        print("Bam file must be indexed.")
        exit(1)

    chroms = []
    chrom_set = set()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            chrom = el[0]
            if chroms == [] or chrom != chroms[-1]:
                assert(chrom not in chrom_set)
                chroms.append(chrom)
                chrom_set.add(chrom)


    for curr_chrom in chroms:

        print("Linking 10X fragments on chromosome: {}".format(curr_chrom))
        flist = read_fragment_matrix(hairs_file,vcf_file,chrom_filter=curr_chrom)

        barcode_to_flist = defaultdict(list)

        for f in flist:

            read_id = f.name if f.name[-3:] != '_MP' else f.name[:-3]

            if f.barcode != None:
                barcode_to_flist[f.barcode].append(f)

        print("  reading bedfile...")
        print("  generating new fragments for HAIRs in boundaries...")

        dup_snp_cover = 0

        # start is 0-indexed start point and stop is 1 past the last residue, 0-indexed
        for chrom,start,stop,barcode in get_molecules(bam_file, curr_chrom, dist=dist):
            #print("molecule: {} {} {} {}".format(chrom,start,stop,barcode))
            if chrom != curr_chrom:
                continue

            if barcode not in barcode_to_flist:
                continue

            barcode_flist = barcode_to_flist[barcode]

            seen_snps = defaultdict(int)
            bad_snps = set()
            new_fseq = []
            for f in barcode_flist:

                if f.used:
                    continue

                used = False
                for (snp_ix, genome_ix, allele_call, qual) in f.seq:

                    if snp_ix in bad_snps:
                        used = True
                        seen_snps[snp_ix] += 1
                        continue # this snp has a disagreement between fragments

                    if genome_ix >= start and genome_ix < stop:

                        used = True

                        if snp_ix in seen_snps:
                            seen_snps[snp_ix] += 1

                            # find location in the fragment we're building,
                            # that has this double-covered SNP
                            # to resolve the conflict
                            for i in range(len(new_fseq)):
                                if new_fseq[i][0] == snp_ix:

                                    if new_fseq[i][2] == allele_call:
                                        q1 = ord(qual) - 33
                                        q2 = ord(new_fseq[i][3]) - 33

                                        Q = q1 + q2 # combined quality score
                                        if Q > 93:
                                            Q = 93

                                        Q_char = chr(33 + Q)

                                        new_fseq[i] = (new_fseq[i][0], new_fseq[i][1], new_fseq[i][2], Q_char)

                                    else:

                                        del new_fseq[i]
                                        bad_snps.add(snp_ix)

                                    break

                        else:

                            new_fseq.append((snp_ix, genome_ix, allele_call, qual))
                            seen_snps[snp_ix] += 1

                f.used = used

            for k,v in seen_snps.items():
                if v > 1:
                    dup_snp_cover += v - 1

            new_fseq.sort(key=lambda x: x[0])
            prev_snp_ix = -1
            for (snp_ix, genome_ix, allele_call, qual) in new_fseq:
                if snp_ix <= prev_snp_ix:
                    #import pdb
                    #pdb.set_trace()
                    print("ERROR",file=sys.stderr)
                    exit(1)
                prev_snp_ix = snp_ix

            new_id = "{}:{}-{}:{}".format(chrom, start+1, stop+1, barcode)
            f = fragment(new_fseq,new_id,None,dtype=0)
            if (len(f.seq) == 0) or (not single_SNP_frags and len(f.seq) < 2):
                continue

            firstpos = f.seq[0][0]
            lines.append((firstpos, str(f)))

        linked_count = 0
        null_count = 0
        unlinked_count = 0
        for f in flist:

            if not f.used:
                if f.barcode == None:
                    null_count += 1
                else:
                    unlinked_count += 1
            else:
                if f.barcode == None:
                    print("linked fragment with no barcode")
                    exit(1)
                linked_count += 1

            if not f.used and (len(f.seq) >= 2 or (single_SNP_frags and len(f.seq) >= 1)):
                firstpos = f.seq[0][0]
                f.dtype = 0
                lines.append((firstpos, str(f)))
                f.used = True



        #print("  {} fragments linked to larger molecules".format(linked_count))
        #print("  {} unlinked fragments with barcodes".format(unlinked_count))
        #print("  {} unlinked fragments without barcodes".format(null_count))
        #print("  {} duplicate snp-cover".format(dup_snp_cover))

        del flist

    # write new fragments to file

    lines.sort()
    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)

def parseargs():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--fragments', nargs='?', type = str, help='file with unlinked hapcut2 fragments (generate using --10X 1 option in extractHAIRS)')
    parser.add_argument('-v', '--vcf_file', nargs='?', type = str, help='vcf file for phasing')
    parser.add_argument('-b', '--bam_file', nargs='?', type = str, help='bam file with barcoded reads')
    parser.add_argument('-o', '--outfile', nargs='?', type = str, help='output file with linked fragments')
    parser.add_argument('-d', '--distance', nargs='?', type = int, help='distance in base pairs that delineates separate 10X molecules',default=20000)

    parser.add_argument('-s', '--single_SNP_frags', action='store_true', help='whether to keep fragments overlapping only one SNP', default=False)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# main function
# parse input and run function to call alleles
if __name__ == '__main__':
    args = parseargs()

    link_fragments(args.fragments,args.vcf_file,args.bam_file, args.outfile, args.distance, args.single_SNP_frags)
