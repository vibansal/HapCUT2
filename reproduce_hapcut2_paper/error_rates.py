# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:07:17 2015

@author: peter
"""

# imports
from collections import defaultdict
import fileIO
import os
import random
#import itertools


def create_genomic_ix(hapblock_list, vcf_file):

    snp_ix = 0
    approx_len = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

            approx_len = genomic_pos

    new_hapblock_list = []

    for blk in hapblock_list:
        new_blk = [(snp_ix, vcf_dict[snp_ix], a1, a2) for (snp_ix, a1, a2) in blk]
        new_hapblock_list.append(new_blk)

    return new_hapblock_list, approx_len


# assume hapblock list only uses genomic index
# add on a SNP index to it for the computation of AN50
def create_SNP_ix(hapblock_list, vcf_file):

    snp_ix = 0
    approx_len = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[genomic_pos] = snp_ix
            snp_ix += 1
            approx_len = genomic_pos

    new_hapblock_list = []

    for blk in hapblock_list:
        new_blk = [(vcf_dict[genomic_pos], genomic_pos, a1, a2) for (genomic_pos, a1, a2) in blk]
        new_hapblock_list.append(new_blk)

    return new_hapblock_list, approx_len

# compute the MEC score of two haplotypes
# RECYCLES PARSING CODE UNNECESSARILY, (not worth updating)?
def compute_MEC(frag_file, hapblock_file, num_snps):

    # it's simplest just to reparse the hapblock file here as a list
    hap1 = ['-']*num_snps
    hap2 = ['-']*num_snps
    try:
        with open(hapblock_file, 'r') as hbf:

            for line in hbf:

                elements = line.strip().split('\t')
                if 'BLOCK' in line or len(elements) < 3:
                    continue

                pos = int(elements[0])-1
                hap1[pos] = elements[1]
                hap2[pos] = elements[2]
    except FileNotFoundError:
        # most of the time, this should mean that the program timed out and therefore didn't produce a phase.
        pass

    MEC = 0

    with open(frag_file, 'r') as ff:
        for line in ff:
            el = line.strip().split()        # line elements
            if len(el) == 0:
                continue # empty line, such as at end of file
            num_blks = int(el[0])
            assert len(el) >= 3 + 2*num_blks     # technically should have qual score, but for our purposes do not enforce

            for i in range(0, num_blks):
                offset   = int(el[2+2*i]) - 1 # convert to 0-index
                seq      = el[3+2*i]
                l        = len(seq)
                hap1_ref = hap1[offset:(offset+l)]   # relevant slice of hap1
                hap2_ref = hap2[offset:(offset+l)]   # relevant slice of hap1
                hap1_EC  = 0
                hap2_EC  = 0
                for x,y in zip(seq, hap1_ref):
                    if x == '1' and y == '0' or x == '0' and y == '1':
                        hap1_EC += 1
                for x,y in zip(seq, hap2_ref):
                    if x == '1' and y == '0' or x == '0' and y == '1':
                        hap2_EC += 1
                MEC += min(hap1_EC, hap2_EC)
    return MEC

# returns a binary list representing which positions are covered or not.
def find_covered_positions(frag_file, num_snps):

    covered = [0]*num_snps

    with open(frag_file, 'r') as infile:
        for line in infile:
            el = line.strip().split()
            num_blks_oldformat = int((len(el)-3)/2)
            num_blks_newformat = int((len(el)-6)/2)

            if num_blks_oldformat == int(el[0]):
                # old format
                for i in range(0,num_blks_oldformat):
                    pos = int(el[2*i+2])-1 # SNPs are 1-indexed
                    seq = el[2*i+3]
                    for j, seq_base in enumerate(seq):
                        snp_ix = pos + j
                        covered[snp_ix] = 1
                        # move along on qual string
            elif num_blks_newformat == int(el[0]):
                # new format

                for i in range(0,num_blks_newformat):
                    pos = int(el[2*i+5])-1 # SNPs are 1-indexed
                    seq = el[2*i+6]
                    for j, seq_base in enumerate(seq):
                        snp_ix = pos + j
                        covered[snp_ix] = 1
                        # move along on qual string

    return covered

# this function is needed for "counting ahead" at beginning of blocks.
# error_rate() needs to properly combine switch errors into mismatches
# such that switch errors are minimized (basically, if a block begins
# with an odd number of consecutive switch errors, it should assume
# that this is a sequence of all mismatches and not a "1-less" sequence of mismatches
# with a switch error at the end of it.
def count_consecutive_switches(t_dict, hap, allele):
    count = 0
    first_SNP = True
    switched = False

    for i, a1, a2 in hap:
        x = t_dict[i]                        # base in true haplotype
        y = a1 if allele == 0 else a2   # base in assembled haplotype
        if x == '-' or y == '-':
            if first_SNP:
                continue
            else:
                break
        elif first_SNP:
            switched = (t_dict[i] != y)
            first_SNP = False
        elif (x != y and not switched) or (x == y and switched):
            count += 1
            switched = not switched
        else:
            break
    return count

# compute a list of the distances between adjacent switch errors

def distance_between_switches(t_blocklist, a_blocklist):

    # iterate over SNPs in the true and assembled haplotypes in parallel
    # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.

    distance_list = []

    for t_block in t_blocklist:
        # convert t_block to a dict for convenience
        t_dict = defaultdict(lambda: '-')
        for i, a1, a2 in t_block:
            t_dict[i] = a1

        for a_block in a_blocklist:

            # first we see which complement of the haplotype is most correct
            a1count = 0 # count of how many correct SNPs
            for i, a1, a2 in a_block:
                y = a1
                x = t_dict[i]

                if x == '-' or y == '-':
                    continue

                if x == y:
                    a1count += 1

            a2count = 0 # count of how many correct SNPs
            for i, a1, a2 in a_block:
                y = a2
                x = t_dict[i]

                if x == '-' or y == '-':
                    continue

                if x == y:
                    a2count += 1


            using_a1 = True if a1count > a2count else False
            switched = False
            distance_since_switch = 0

            for i, a1, a2 in a_block:
                y = a1 if using_a1 else a2 # select our correct haplotype
                x = t_dict[i]

                if x == '-' or y == '-':
                    continue

                if x == y:
                    if switched: # we were switched, now we're flipping back
                        distance_list.append(distance_since_switch)
                        distance_since_switch = 0
                        switched = False
                else:
                    switched = True
                    distance_since_switch += 1

    return distance_list

# combine two dicts
def merge_dicts(d1,d2):
    d3 = d2.copy()
    for k, v in d1.items():
        assert(k not in d3)
        d3[k] = v
    return d3

# the "error_result" abstraction and its overloaded addition operator are handy
# for combining results for the same chromosome across blocks (when the "ground truth"
# is a set of blocks rather than trio), and combining results across different chromosomes
# into genome wide stats
class error_result():
    def __init__(self, tool_name=None,dataset_name=None,ref=None,switch_count=None,poss_sw=None,mismatch_count=None,poss_mm=None,flat_count=None,poss_flat=None,
                 phased_count=None,num_covered=None,num_snps=None,maxblk_snps=None,approx_len=None,runtime=None,
                 AN50_spanlst=None,N50_spanlst=None,switch_loc=None,mismatch_loc=None,missing_loc=None):

        def create_dict(val,d_type):
            new_dict = defaultdict(d_type)
            if ref != None and val != None:
                new_dict[ref] = val
            return new_dict

        self.ref          = set() # set of references in this result (e.g. all chromosomes)
        if ref != None:
            self.ref.add(ref)

        self.tool_name = tool_name
        self.dataset_name = dataset_name

        # these are things that can be summed for the same reference,
        # e.g. switch counts for separate blocks are additive
        self.switch_count   = create_dict(switch_count,   int)
        self.poss_sw        = create_dict(poss_sw,        int)
        self.mismatch_count = create_dict(mismatch_count, int)
        self.poss_mm        = create_dict(poss_mm,        int)
        self.flat_count     = create_dict(flat_count,     int)
        self.poss_flat      = create_dict(poss_flat,      int)
        self.phased_count   = create_dict(phased_count,   int)
        self.AN50_spanlst   = create_dict(AN50_spanlst,   list)
        self.N50_spanlst    = create_dict(N50_spanlst,    list)

        # these are things that are non-additive properties, because they
        # refer to the whole reference and would be double-counted
        # e.g. if we combine errors for two blocks, on same chromosome, we add their errors
        # but we can't just add "num_snps", their chromosomes' total snp counts
        # so we use dictionaries to make sure these properties aren't duplicated

        self.num_covered  = create_dict(num_covered, int)
        self.num_snps     = create_dict(num_snps,    int)
        self.maxblk_snps  = create_dict(maxblk_snps, int)
        self.approx_lens  = create_dict(approx_len,  int)
        self.runtime      = create_dict(runtime,     int)

        self.switch_loc   = create_dict(switch_loc,   list)
        self.mismatch_loc = create_dict(mismatch_loc, list)
        self.missing_loc  = create_dict(missing_loc,  list)

    # combine two error rate results
    def __add__(self,other):
        new_err                = error_result()

        if self.tool_name == None and other.tool_name == None:
            new_err.tool_name = None
        elif self.tool_name == None and other.tool_name != None:
            new_err.tool_name = other.tool_name
        elif self.tool_name != None and other.tool_name == None:
            new_err.tool_name = self.tool_name
        elif self.tool_name != None and other.tool_name != None:
            assert(self.tool_name == other.tool_name)
            new_err.tool_name = self.tool_name

        if self.dataset_name == None and other.dataset_name == None:
            new_err.dataset_name = None
        elif self.dataset_name == None and other.dataset_name != None:
            new_err.dataset_name = other.dataset_name
        elif self.dataset_name != None and other.dataset_name == None:
            new_err.dataset_name = self.dataset_name
        elif self.dataset_name != None and other.dataset_name != None:
            assert(self.dataset_name == other.dataset_name)
            new_err.dataset_name = self.dataset_name

        new_err.ref            = self.ref.union(other.ref)
        new_err.switch_count   = merge_dicts(self.switch_count,   other.switch_count)
        new_err.poss_sw        = merge_dicts(self.poss_sw,        other.poss_sw)
        new_err.mismatch_count = merge_dicts(self.mismatch_count, other.mismatch_count)
        new_err.poss_mm        = merge_dicts(self.poss_mm,        other.poss_mm)
        new_err.flat_count     = merge_dicts(self.flat_count,     other.flat_count)
        new_err.poss_flat      = merge_dicts(self.poss_flat,      other.poss_flat)
        new_err.phased_count   = merge_dicts(self.phased_count,   other.phased_count)
        new_err.AN50_spanlst   = merge_dicts(self.AN50_spanlst,   other.AN50_spanlst)
        new_err.N50_spanlst    = merge_dicts(self.N50_spanlst,    other.N50_spanlst)
        new_err.num_covered    = merge_dicts(self.num_covered,    other.num_covered)
        new_err.num_snps       = merge_dicts(self.num_snps,       other.num_snps)
        new_err.maxblk_snps    = merge_dicts(self.maxblk_snps,    other.maxblk_snps)
        new_err.approx_lens    = merge_dicts(self.approx_lens,    other.approx_lens)
        new_err.runtime        = merge_dicts(self.runtime,        other.runtime)
        new_err.switch_loc     = merge_dicts(self.switch_loc,     other.switch_loc)
        new_err.mismatch_loc   = merge_dicts(self.mismatch_loc,   other.mismatch_loc)
        new_err.missing_loc    = merge_dicts(self.missing_loc,    other.missing_loc)
        return new_err

    def update_runtime(self, ref, runtime_file):
        if os.path.isfile(runtime_file):
            self.runtime[ref] = fileIO.parse_runtime_file(runtime_file)
        else:
            self.runtime[ref] = 0

    def get_num_covered(self):
        return sum(self.num_covered.values())

    def get_switch_count(self):
        return sum(self.switch_count.values())

    def get_mismatch_count(self):
        return sum(self.mismatch_count.values())

    def get_flat_count(self):
        return sum(self.flat_count.values())

    def get_poss_sw(self):
        return sum(self.poss_sw.values())

    def get_poss_mm(self):
        return sum(self.poss_mm.values())

    def get_poss_flat(self):
        return sum(self.poss_flat.values())

    def get_num_snps(self):
        return sum(self.num_snps.values())

    def get_approx_len(self):
        return sum(self.approx_lens.values())

    def get_phased_count(self):
        return sum(self.phased_count.values())

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        poss_sw = self.get_poss_sw()

        if poss_sw > 0:
            return switch_count/poss_sw
        else:
            return 0

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return mismatch_count/poss_mm
        else:
            return 0


    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return (self.get_switch_count() + self.get_mismatch_count())/poss_mm
        else:
            return 0

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return flat_count/poss_flat
        else:
            return 0

    def get_missing_rate(self):
        num_cov = self.get_num_covered()
        if num_cov > 0:
            return 1-sum(self.phased_count.values())/num_cov
        else:
            return 0

    def get_AN50(self):
        AN50 = 0
        AN50_spanlst = sum(self.AN50_spanlst.values(),[])
        AN50_spanlst.sort(reverse=True)
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > self.get_num_snps()/2:
                AN50 = span
                break
        return AN50

    def get_N50(self):
        N50  = 0
        N50_spanlst = sum(self.N50_spanlst.values(),[])
        N50_spanlst.sort(reverse=True)

        total = 0
        for span in N50_spanlst:
            total += span
            if total > self.get_approx_len()/2:
                N50 = span
                break
        return N50

    def get_max_blk_snp_percent(self):
        snps_in_max_blks = sum(self.maxblk_snps.values())
        sum_all_snps     = self.get_num_snps()

        if sum_all_snps > 0:
            return snps_in_max_blks / sum_all_snps
        else:
            return 0

    def get_runtime(self):
        return sum(self.runtime.values())

    def __str__(self):
        s = ('''
tool:            {}
dataset:         {}
switch rate:     {}
mismatch rate:   {}
flat rate:       {}
missing rate:    {}
switch errors:   {}
poss. switch:    {}
mismatch errors: {}
poss. mismatch:  {}
flat errors:     {}
poss. flat:      {}
phased count:    {}
num covered:     {}
AN50:            {}
N50:             {}
max blk snp %:   {}
runtime:         {}
        '''.format(self.tool_name, self.dataset_name, self.get_switch_rate(), self.get_mismatch_rate(),
                   self.get_flat_error_rate(), self.get_missing_rate(),
                   self.get_switch_count(), self.get_poss_sw(),
                   self.get_mismatch_count(), self.get_poss_mm(),
                   self.get_flat_count(), self.get_poss_flat(),
                   self.get_phased_count(), self.get_num_covered(),
                   self.get_AN50(),self.get_N50(),self.get_max_blk_snp_percent(),self.get_runtime()))
        return s

# compute error rates by using another haplotype block file as ground truth
def hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file=None, use_SNP_index=True, tool_name=None, dataset_name=None):

    # parse and get stuff to compute error rates
    t_blocklist = fileIO.parse_hapblock_file(truth_file,use_SNP_index)
    a_blocklist = fileIO.parse_hapblock_file(assembly_file,use_SNP_index)
    # compute error result object, update the runtime and AN50 / completeness
    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, runtime_file, use_SNP_index, tool_name=tool_name, dataset_name=dataset_name)
    return(err)

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, runtime_file=None, use_SNP_index=True, tool_name=None, dataset_name=None,largest_blk_only=False):

    # parse and get stuff to compute error rates
    t_blocklist = fileIO.parse_vcf_phase(vcf_file,use_SNP_index)
    a_blocklist = fileIO.parse_hapblock_file(assembly_file,use_SNP_index)
    # compute error result object, update the runtime and AN50 / completeness

    if largest_blk_only:
        largest_blk = []
        for blk in a_blocklist:
            if len(blk) > len(largest_blk):
                largest_blk = blk

        a_blocklist = [largest_blk]

    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, runtime_file, use_SNP_index,tool_name=tool_name, dataset_name=dataset_name)
    return(err)

def hapblock_vcf_error_rate_COMMON_ALL_CHROM(tool_list,A_list,R_list,frag_files,vcf_files,dataset_name):
    # there are m chrom and n tools
    # A_list and R_list are m long, and each index has a list of files n long for each tool on that chrom
    # A is assemblies, R is runtimes
    # frag_files and vcf_files are m long (one per chrom)
    curr_results = dict()
    for tool in tool_list:
        curr_results[tool] = error_result()

    for assembly_list, runtime_files,frag_file,vcf_file in zip(A_list,R_list,frag_files,vcf_files):
        result = hapblock_vcf_error_rate_COMMON(assembly_list,runtime_files,frag_file,vcf_file,tool_names=tool_list, dataset_name=dataset_name)
        for res, tool in zip (result,tool_list):
            curr_results[tool] += res
    print("=====================================================")
    print("ERROR RESULTS ACROSS ALL CHROM, FOR COMMON PHASED SNPS")
    print("=====================================================")

    for tool in tool_list:
        print(curr_results[tool])

    #lst = [curr_results[tool] for tool in tool_list]
    #return lst

    return curr_results

def hapblock_vcf_error_rate_COMMON(assembly_list,runtime_files,frag_file,vcf_file,tool_names=None, dataset_name=None):
    phase_set = set() # set of positions phased
    # tool list is list of strings n long, and there are m chrom
    a_blocklist_list = [fileIO.parse_hapblock_file(a,use_SNP_index=True) for a in assembly_list]

    for i, a_blocklist in enumerate(a_blocklist_list):
        # conglomerate all the phased spots in first file
        if i == 0:
            for blk in a_blocklist:
                for pos, a1, a2 in blk:
                    if a1 in ['0','1'] and a1 != a2:
                        phase_set.add(pos)
        else:
            # take the intersection with all the remaining files
            curr_phased = set()
            for blk in a_blocklist:
                for pos, a1, a2 in blk:
                    if a1 in ['0','1'] and a1 != a2:
                        curr_phased.add(pos)

            phase_set = phase_set.intersection(curr_phased)

    res = []

    t_blocklist = fileIO.parse_vcf_phase(vcf_file,use_SNP_index=True)

    for a_blocklist, runtime_file,tool in zip(a_blocklist_list,runtime_files,tool_names):
        err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, runtime_file, use_SNP_index=True, phase_set=phase_set, tool_name=tool, dataset_name=dataset_name)
        res.append(err)

    return res

def read_fragment_matrix(frag_matrix):

    names = []
    flist = []

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            name = el[1]

            call_list  = el[2:(2+2*num_blks)]              # extract base call part of line
            call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            #qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            alist= [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)]

            flist.append(alist)
            names.append(name)

    zipped = zip(names,flist)
    sorted_zipped = sorted(zipped,key=lambda x: x[1][0][0])
    sorted_names, sorted_flist = [list(z) for z in zip(*sorted_zipped)]

    return sorted_flist, sorted_names


# num_covered should be the number of SNPs with coverage in the fragment matrix file
# debug returns extra lists with indexes of errors
def error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file=None, runtime_file=None, use_SNP_index=True, phase_set=None, tool_name=None, dataset_name=None):

    ref_name    = fileIO.get_ref_name(vcf_file)
    num_snps = fileIO.count_SNPs(vcf_file)
    num_covered = sum(find_covered_positions(frag_file, num_snps)) if frag_file != None else 0

    if use_SNP_index:
        a_blocklist_double_index, approx_len = create_genomic_ix(a_blocklist, vcf_file)
    else:
        a_blocklist_double_index, approx_len = create_SNP_ix(a_blocklist, vcf_file)

    switch_count   = 0
    mismatch_count = 0
    poss_sw        = 0 # count of possible positions for switch errors
    poss_mm        = 0 # count of possible positions for mismatches
    phased_count   = 0
    maxblk_snps    = 0
    switch_loc     = []
    mismatch_loc   = []
    missing_loc    = []
    AN50_spanlst   = []
    N50_spanlst    = []

    for blk in a_blocklist_double_index:

        first_pos  = -1
        last_pos   = -1
        first_SNP  = -1
        last_SNP   = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2 in blk:

            if a1 == '-' or (phase_set != None and snp_ix not in phase_set):
                missing_loc.append(pos)
            else:
                phased_count += 1

                blk_phased+=1
                if first_pos == -1:
                    first_pos = pos
                    first_SNP = snp_ix
                last_pos = pos
                last_SNP = snp_ix

        blk_total = last_SNP - first_SNP + 1

        AN50_spanlst.append(((last_pos-first_pos)*(blk_phased/blk_total), blk_phased))
        N50_spanlst.append((last_pos-first_pos))

        if blk_phased > maxblk_snps:
            maxblk_snps = blk_phased

    for t_block in t_blocklist:

        switched       = False
        last_base_was_switch = False

        # convert t_block to a dict for convenience
        t_dict = defaultdict(lambda: '-')
        for i, a1, a2 in t_block:
            t_dict[i] = a1

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

            blk_switches    = [0,0]
            blk_mismatches  = [0,0]
            blk_switchlist  = [[],[]]
            blk_mmlist      = [[],[]]
            for a in [0,1]: # choose which allele to score. this only makes a difference for minimizing switch errors vs mismatches in corner cases.

                first_SNP = True
                for blk_ix, (i, a1, a2) in enumerate(a_block):
                    y = a1 if a == 0 else a2
                    x = t_dict[i]

                    if x == '-' or y == '-' or (phase_set != None and i not in phase_set):
                        continue

                    if first_SNP:
                        switched = (x != y)
                        if count_consecutive_switches(t_dict, a_block[blk_ix:], a) % 2 == 1:
                            last_base_was_switch = True
                        else:
                            last_base_was_switch = False
                        first_SNP = False
                        continue

                    # if there is a mismatch against the true haplotype and we are in a normal state,
                    # or if there is a "match" that isn't a match because we are in a switched state,
                    # then we need to flip the state again and iterate the count
                    if (x != y and not switched) or (x == y and switched): # current base is mismatched, implying a switch
                        switched = not switched                  # flip the "switched" status

                        if last_base_was_switch:                 # if last base was a switch then this is actually a single-base mismatch
                            # count the 2 switches as a single-base mismatch instead
                            blk_mismatches[a] += 1
                            blk_mmlist[a].append(i)
                            blk_switches[a] -= 1      # undo count from last base switch
                            if len(blk_switchlist[a]) > 0:
                                blk_switchlist[a].pop()
                            if (blk_switches[a] < 0):
                                blk_switches[a] = 0
                            last_base_was_switch = False

                        else:

                            blk_switches[a] += 1
                            blk_switchlist[a].append(i)
                            last_base_was_switch = True

                    else: # current base is not mismatched
                        last_base_was_switch = False

                # special case for switch on last base of previous a_block; should count as a mismatch
                if last_base_was_switch:
                    # count the switch as a single-base mismatch instead
                    blk_mismatches[a] += 1
                    blk_mmlist[a].append(i)
                    blk_switches[a] -= 1
                    if len(blk_switchlist[a]) > 0:
                        blk_switchlist[a].pop()

                    if (blk_switches[a] < 0):
                        blk_switches[a] = 0

            if blk_switches[0] < blk_switches[1]:
                switch_count   += blk_switches[0]
                mismatch_count += blk_mismatches[0]
                switch_loc += blk_switchlist[0]
                mismatch_loc += blk_mmlist[0]

            else:
                switch_count   += blk_switches[1]
                mismatch_count += blk_mismatches[1]
                switch_loc     += blk_switchlist[1]
                mismatch_loc   += blk_mmlist[1]

        assert(len(switch_loc) == switch_count)
        assert(len(mismatch_loc) == mismatch_count)

        # tally up how many possible positions there are for switch errors and mismatches
        # count how many phased SNPs there are so we can calculate a rate of pruned SNPs



        for blk in a_blocklist:
            phased_known = 0
            for i, a1, a2 in blk:
                if t_dict[i] != '-' and a1 != '-' and (phase_set == None or i in phase_set):
                    phased_known += 1

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between SNPs counts, and -2 for the two ends.
            if phased_known >= 4:
                poss_sw += (phased_known - 3)
            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                poss_mm += phased_known

        poss_flat  = poss_mm
        flat_count = 0

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

            flat_count1 = 0
            flat_count2 = 0

            for i, a1, a2 in a_block:

                if a1 == '-' or a2 == '-' or t_dict[i] == '-' or (phase_set != None and i not in phase_set):
                    continue

                if (a1 != t_dict[i]):
                    flat_count1 += 1
                if (a2 != t_dict[i]):
                    flat_count2 += 1

            if flat_count1 < flat_count2:
                flat_count += flat_count1
            else:
                flat_count += flat_count2

    runtime = -1
    if runtime_file != None:
        runtime = fileIO.parse_runtime_file(runtime_file)

    total_error = error_result(ref=ref_name,tool_name=tool_name,dataset_name=dataset_name,
             switch_count=switch_count,poss_sw=poss_sw, mismatch_count=mismatch_count,
             poss_mm=poss_mm,flat_count=flat_count,poss_flat=poss_flat,
             phased_count=phased_count,num_covered=num_covered,num_snps=num_snps,
             maxblk_snps=maxblk_snps,approx_len=approx_len,runtime=runtime,
             AN50_spanlst=AN50_spanlst,N50_spanlst=N50_spanlst,switch_loc=switch_loc,
             mismatch_loc=mismatch_loc,missing_loc=missing_loc)

    return total_error



class PDER_result:

    def __init__(self,scores, counts, binsize, sample_frac, max_dist):
        self.scores = scores
        self.counts = counts
        self.binsize = binsize
        self.sample_frac = sample_frac
        self.max_dist = max_dist

    def get_error_stats(self):
        result = []
        for s, c in zip(self.scores, self.counts):
            if c > 0:
                result.append(s/c)
            else:
                result.append(0)

        return result

    def get_bin_starts(self):
        d = int(self.max_dist / self.binsize) # get the max possible distance bin
        return [i * self.binsize for i in range(0, d+1)]

        # combine two error rate results
    def __add__(self,other):
        assert(self.binsize == other.binsize)
        assert(self.sample_frac == other.sample_frac)
        max_dist = max(self.max_dist, other.max_dist)
        l = max(len(self.scores), len(other.scores))
        scores = [0.0]*l
        for i, s in enumerate(self.scores):
            scores[i] += s
        for i, s in enumerate(other.scores):
            scores[i] += s

        counts = [0.0]*l
        for i, c in enumerate(self.counts):
            counts[i] += c
        for i, c in enumerate(other.counts):
            counts[i] += c

        return PDER_result(scores,counts,self.binsize, self.sample_frac, max_dist)

# new error rate calculation 4/18/2016

def per_distance_error_rate(assembly_file, vcf_file, binsize=100000, sample_frac=0.01):

    t_blocklist = fileIO.parse_vcf_phase(vcf_file,use_SNP_index=False)
    a_blocklist = fileIO.parse_hapblock_file(assembly_file,use_SNP_index=False)

    return per_distance_error_rate_general(t_blocklist, a_blocklist, binsize, sample_frac)

# new error rate calculation 4/18/2016

def per_distance_error_rate_general(t_blocklist, a_blocklist, binsize=100000, sample_frac=0.05):

    # get maximum distance between any pair of SNPs (length of largest block)
    max_dist = 0

    for blk in a_blocklist:
        dist = blk[-1][0] - blk[0][0]
        if dist > max_dist:
            max_dist = dist


    d = int(max_dist / binsize) # get the max possible distance bin

    scores = [0]*(d+1)
    counts = [0]*(d+1)

    # for each block of the ground truth haplotype
    for t_block in t_blocklist:

        # convert truth block to a dict for convenience
        t_dict = defaultdict(lambda: '-')
        for i, a1, a2 in t_block:
            t_dict[i] = a1

        # for each block of the assembled haplotype
        for a_block in a_blocklist:
            l = len(a_block)

            # for each distance-between-SNPs that we can observe

            # for each pair of SNP positions in our sample
            samplesize = int(l**2 * sample_frac)

            for k in range(0,samplesize):
                i = random.randrange(l)
                j = random.randrange(l)
                if i > j:
                    temp = i
                    i = j
                    j = temp


                pos1, a1, foo1 = a_block[i]
                pos2, a2, foo2 = a_block[j]

                # ground truth phase is not known between these SNPs so skip
                if t_dict[pos1] == '-' or t_dict[pos2] == '-':
                    continue
                if a1 == '-' or a2 == '-':                         # assembled phase not called between these SNPs
                    continue

                dist = int((pos2 - pos1)/binsize)
                # count[] holds the possible max score


                if a1 == '-' or a2 == '-':                         # assembled phase not called between these SNPs
                    continue

                counts[dist] += 1


                if (t_dict[pos1] == t_dict[pos2]) == (a1 == a2): # phase matches ground truth
                    scores[dist] += 1

    return PDER_result(scores, counts, binsize, sample_frac, max_dist)

# assess completeness of a haplotype
# returns (AN50, % total SNPs in largest block)
def frac_SNPs_per_num_blks(hapblock_list, vcf_file, use_SNP_index=True):

    num_snps = fileIO.count_SNPs(vcf_file)

    if use_SNP_index:
        hapblock_list, approx_len = create_genomic_ix(hapblock_list, vcf_file)
    else:
        hapblock_list, approx_len = create_SNP_ix(hapblock_list, vcf_file)

    phased_count_list = []

    for BLK in hapblock_list:
        blk_phased = 0

        for snp_ix, pos, a1, a2 in BLK:

            if a1 != '-':
                blk_phased+=1

        phased_count_list.append(blk_phased)

    phased_count_list.sort(reverse=True)

    res_list = []
    num_blks = list(range(1,len(phased_count_list)+1))

    for x in num_blks:
        val = sum(phased_count_list[:x])/num_snps
        res_list.append(val)

    print("{} SNPs in largest block vs. {} SNPs in all others".format(res_list[0],res_list[-1]-res_list[0]))

    return num_blks, res_list
