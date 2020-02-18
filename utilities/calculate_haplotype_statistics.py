#!/usr/bin/env python3

from __future__ import print_function
# test
# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

## edited 02/10/2018 to print 'fraction of SNVs phased and to allow for chr7 to match with '7' in different VCF files

# imports
from collections import defaultdict
import argparse
import sys
import statistics

desc = '''
Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.
Error rates for an assembled haplotype (specified by -v1 and optionally -h1 arguments)
are computed with respect to a "reference" haplotype (specified by -v2 and optionally -h2 arguments).
All files must contain information for one chromosome only!
To compute aggregate statistics across multiple chromosomes, provide files for
each chromosome/contig as an ordered list, using the same chromosome order between flags.

Note: Triallelic variants are supported, but variants with more than 2 alternative alleles
are currently NOT supported. These variants are ignored. Also, variants where the ref and alt
alleles differ between the test haplotype and reference haplotype are skipped.
'''

def parse_args():

    parser = argparse.ArgumentParser(description=desc)
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-v1', '--vcf1', nargs='+', type = str, help='A phased, single sample VCF file to compute haplotype statistics on.')
    parser.add_argument('-v2', '--vcf2', nargs='+', type = str, help='A phased, single sample  VCF file to use as the "ground truth" haplotype.')
    parser.add_argument('-h1', '--haplotype_blocks1', nargs='+', type = str, help='Override the haplotype information in "-v1" with the information in this HapCUT2-format haplotype block file. If this option is used, then the VCF specified with -v1 MUST be the same VCF used with HapCUT2 (--vcf) to produce the haplotype block file!')
    parser.add_argument('-h2', '--haplotype_blocks2', nargs='+', type = str, help='Override the haplotype information in "-v2" with the information in this HapCUT2-format haplotype block file. If this option is used, then the VCF specified with -v2 MUST be the same VCF used with HapCUT2 (--vcf) to produce the haplotype block file!')
    parser.add_argument('-i', '--indels', action="store_true", help='Use this flag to consider indel variants. Default: Indels ignored.',default=False)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def parse_hapblock_file(hapblock_file,vcf_file,indels=False):

    snp_ix = 0
    vcf_dict = dict()
    CHROM = None
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            chrom = el[0]

            if CHROM == None:
                CHROM = chrom
            elif chrom != CHROM:
                print("ERROR: VCFs should contain one chromosome per file")
                print("VCF file:             " + vcf_file)
                exit(1)

            consider = True

            a0 = el[3]
            a1 = el[4]
            a2 = None

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    consider = False


            genotype = el[9].split(':')[0]

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if consider and (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                consider = False

            genomic_pos = int(el[1])-1
            if consider:
                vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

    blocklist = [] # data will be a list of blocks, each with a list tying SNP indexes to haplotypes

    with open(hapblock_file, 'r') as hbf:
        for line in hbf:
            if len(line) < 3: # empty line
                continue
            if 'BLOCK' in line:
                blocklist.append([])
                continue

            el = line.strip().split('\t')
            if len(el) < 3: # not enough elements to have a haplotype
                continue

            snp_ix = int(el[0])-1

            if snp_ix not in vcf_dict:
                continue

            pos = vcf_dict[snp_ix]

            if len(el) >= 5:

                chrom = el[3]
                if chrom != CHROM and chrom.lstrip('chr') != CHROM.lstrip('chr'):
                    print("ERROR: Chromosome in haplotype block file doesn't match chromosome in VCF")
                    print("Haplotype block file: " + hapblock_file)
                    print("VCF file:             " + vcf_file)
                    exit(1)

                pos2 = int(el[4])-1
                assert(pos == pos2) # if present, genomic position in haplotype block file should match that from VCF

            allele1 = el[1]
            allele2 = el[2]
            ref = el[5]
            alt1 = el[6]
            alt2 = None

            if ',' in alt1:
                alt_lst = alt1.split(',')
                if len(alt_lst) == 2:
                    alt1,alt2 = alt_lst
                else:
                    continue

            blocklist[-1].append((snp_ix, pos, allele1, allele2, ref, alt1, alt2))

    return [b for b in blocklist if len(b) > 1]

def parse_vcf_phase(vcf_file, CHROM, indels = False):

    #block = []
    PS_index = None
    blocks = defaultdict(list)

    with open(vcf_file, 'r') as vcf:

        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue
            if len(el) != 10:
                print("VCF file must be single-sample.")
                exit(1)

            # get the index where the PS information is
            for i,f in enumerate(el[8].split(':')):
                if i == 0:
                    assert(f == 'GT')
                if f == 'PS':
                    if PS_index == None:
                        PS_index = i
                    else:
                        assert(PS_index == i)
                    break

    if PS_index == None:
        print("WARNING: PS flag is missing from VCF. Assuming that all phased variants are in the same phase block.")

    with open(vcf_file, 'r') as vcf:

        snp_ix = 0

        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue
            if len(el) != 10:
                print("VCF file must be single-sample.")
                exit(1)

            consider = True

            phase_data = el[9]

            a0 = el[3]
            a1 = el[4]
            a2 = None

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    consider = False

            dat = el[9].split(':')
            genotype = dat[0]

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if consider and (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                consider = False

            ps = None
            if PS_index == None:
                ps = 1    # put everything in one block
            elif consider and len(dat) > PS_index:
                ps = dat[PS_index]
                if ps == '.':
                    consider = False

            chrom = el[0]

            if chrom != CHROM and chrom.lstrip('chr') != CHROM.lstrip('chr'):
                print("ERROR: Chromosome in reference haplotype VCF doesn't match chromosome in VCF used for phasing")
                print("reference haplotype VCF: " + vcf_file)
                print("{} != {}".format(CHROM, chrom))
                exit(1)

            pos = int(el[1])-1
            if ps != None and consider and phase_data[1] == '|':
                blocks[ps].append((snp_ix, pos, phase_data[0:1], phase_data[2:3], a0, a1, a2))

            snp_ix += 1

    return [v for k,v in sorted(list(blocks.items())) if len(v) > 1]

# given a VCF file, simply count the number of heterozygous SNPs present.
def count_SNPs(vcf_file,indels=False):

    count = 0

    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            a0 = el[3]
            a1 = el[4]
            a2 = None

            if ',' in a1:
                alt_lst = a1.split(',')
                if len(alt_lst) == 2:
                    a1,a2 = alt_lst
                else:
                    continue

            genotype = el[9][:3]
            #print ("pre: ",line.strip())

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                continue

            if genotype[0] == genotype[2]:
                continue

            if (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                continue

            count += 1

    return count

# given a VCF file from a single chromosome return the name of that chromosome
# will almost always be from a single chromosome but don't assume that
def get_ref_name(vcf_file):

    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[0] == '#':
                continue
            elif len(line.strip().split('\t')) < 5:
                continue
            return line.strip().split('\t')[0]
    # default
    print("ERROR")
    exit(1)

# this function is needed for "counting ahead" at beginning of blocks.
# error_rate() needs to properly combine switch errors into mismatches
# such that switch errors are minimized (basically, if a block begins
# with an odd number of consecutive switch errors, it should assume
# that this is a sequence of all mismatches and not a "1-less" sequence of mismatches
# with a switch error at the end of it.
def count_consecutive_switches(t1_dict, hap, allele):
    count = 0
    first_SNP = True
    switched = False

    for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in hap:
        x = t1_dict[pos]                        # base in true haplotype
        y = a1 if allele == 0 else a2   # base in assembled haplotype
        if x == '-' or y == '-':
            if first_SNP:
                continue
            else:
                break
        elif first_SNP:
            switched = (t1_dict[pos] != y)
            first_SNP = False
        elif (x != y and not switched) or (x == y and switched):
            count += 1
            switched = not switched
        else:
            break
    return count

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
    def __init__(self, ref=None,switch_count=None,poss_sw=None,mismatch_count=None,poss_mm=None,flat_count=None,poss_flat=None,
                 phased_count=None,num_snps=None,maxblk_snps=None,
                 AN50_spanlst=None,N50_spanlst=None,switch_loc=None,mismatch_loc=None):

        def create_dict(val,d_type):
            new_dict = defaultdict(d_type)
            if ref != None and val != None:
                new_dict[ref] = val
            return new_dict

        self.ref          = set() # set of references in this result (e.g. all chromosomes)
        if ref != None:
            self.ref.add(ref)

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

        self.num_snps     = create_dict(num_snps,    int)
        self.maxblk_snps  = create_dict(maxblk_snps, int)

        self.switch_loc   = create_dict(switch_loc,   list)
        self.mismatch_loc = create_dict(mismatch_loc, list)


    # combine two error rate results
    def __add__(self,other):
        new_err = error_result()

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
        new_err.num_snps       = merge_dicts(self.num_snps,       other.num_snps)
        new_err.maxblk_snps    = merge_dicts(self.maxblk_snps,    other.maxblk_snps)
        new_err.switch_loc     = merge_dicts(self.switch_loc,     other.switch_loc)
        new_err.mismatch_loc   = merge_dicts(self.mismatch_loc,   other.mismatch_loc)

        return new_err

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

    def get_phased_count(self):
        return sum(self.phased_count.values())

    # error rate accessor functions
    def get_switch_rate(self):
        switch_count = self.get_switch_count()
        poss_sw = self.get_poss_sw()
        if poss_sw > 0:
            return float(switch_count)/poss_sw
        else:
            return 0

    def get_mismatch_rate(self):
        mismatch_count = self.get_mismatch_count()
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(mismatch_count)/poss_mm
        else:
            return 0


    def get_switch_mismatch_rate(self):
        poss_mm = self.get_poss_mm()

        if poss_mm > 0:
            return float(self.get_switch_count() + self.get_mismatch_count())/poss_mm
        else:
            return 0

    def get_flat_error_rate(self):
        flat_count = self.get_flat_count()
        poss_flat = self.get_poss_flat()
        if poss_flat > 0:
            return float(flat_count)/poss_flat
        else:
            return 0

    def get_AN50(self):
        AN50 = 0
        AN50_spanlst = sum(self.AN50_spanlst.values(),[])
        AN50_spanlst.sort(reverse=True)
        phased_sum = 0
        for span, phased in AN50_spanlst:
            phased_sum += phased
            if phased_sum > self.get_num_snps()/2.0:
                AN50 = span
                break
        return AN50

    def get_N50_phased_portion(self):
        N50  = 0
        N50_spanlst = sum(self.N50_spanlst.values(),[])
        N50_spanlst.sort(reverse=True)

        L = sum(N50_spanlst)

        total = 0
        for span in N50_spanlst:
            total += span
            if total > L/2.0:
                N50 = span
                break
        return N50

    def get_median_block_length(self):
        spanlst = sum(self.N50_spanlst.values(),[])
        return statistics.median(spanlst)

#    def get_max_blk_snp_percent(self):
#        snps_in_max_blks = sum(self.maxblk_snps.values())
#        sum_all_snps     = self.get_num_snps()
#
#        if sum_all_snps > 0:
#            return float(snps_in_max_blks) / sum_all_snps
#        else:
#            return 0

    def __str__(self):

        s = ('''
switch rate:        {}
mismatch rate:      {}
flat rate:          {}
phased count:       {}
AN50:               {}
N50:                {}
num snps max blk:   {}
            '''.format(self.get_switch_rate(), self.get_mismatch_rate(),
                   self.get_flat_error_rate(), self.get_phased_count(),
                   self.get_AN50(), self.get_N50_phased_portion(), sum(self.maxblk_snps.values())))

        return s

# compute error rates by using another haplotype block file as ground truth
def hapblock_hapblock_error_rate_multiple(assembly_files, assembly_vcf_files, truth_files, truth_vcf_files, indels):

    err = error_result()
    for assembly_file, assembly_vcf_file, truth_file, truth_vcf_file in zip(assembly_files, assembly_vcf_files, truth_files, truth_vcf_files):
        err += hapblock_hapblock_error_rate(assembly_file, assembly_vcf_file, truth_file, truth_vcf_file, indels)

    return err

# compute error rates by using another haplotype block file as ground truth
def hapblock_hapblock_error_rate(assembly_file, assembly_vcf_file, truth_file, truth_vcf_file, indels):

    # parse and get stuff to compute error rates
    t_blocklist = parse_hapblock_file(truth_file,truth_vcf_file,indels)
    a_blocklist = parse_hapblock_file(assembly_file,assembly_vcf_file,indels)
    # compute error result object
    err = error_rate_calc(t_blocklist, a_blocklist, assembly_vcf_file, indels)
    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def hapblock_vcf_error_rate_multiple(assembly_files, vcf_files, phased_vcf_files, indels, largest_blk_only=False):

    err = error_result()
    for assembly_file, vcf_file, phased_vcf_file in zip(assembly_files, vcf_files, phased_vcf_files):
        err += hapblock_vcf_error_rate(assembly_file, vcf_file, phased_vcf_file, indels)

    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have phase information
def hapblock_vcf_error_rate(assembly_file, vcf_file, phased_vcf_file, indels, largest_blk_only=False):

    # parse and get stuff to compute error rates
    CHROM = get_ref_name(vcf_file)
    t_blocklist = parse_vcf_phase(phased_vcf_file, CHROM, indels)
    a_blocklist = parse_hapblock_file(assembly_file,vcf_file,indels)
    # compute error result object
    if largest_blk_only:
        largest_blk = []
        for blk in a_blocklist:
            if len(blk) > len(largest_blk):
                largest_blk = blk

        a_blocklist = [largest_blk]

    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, indels)
    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def vcf_hapblock_error_rate_multiple(vcf_files, truth_files, truth_vcf_files, indels, largest_blk_only=False):

    err = error_result()
    for vcf_file, truth_file, truth_vcf_file in zip(vcf_files, truth_files, truth_vcf_files):
        err += vcf_hapblock_error_rate(vcf_file, truth_file, truth_vcf_file, indels)

    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have phase information
def vcf_hapblock_error_rate(vcf_file, truth_file, truth_vcf_file, indels, largest_blk_only=False):

    # parse and get stuff to compute error rates
    CHROM = get_ref_name(vcf_file)
    a_blocklist = parse_vcf_phase(vcf_file, CHROM, indels)
    t_blocklist = parse_hapblock_file(truth_file,truth_vcf_file,indels)
    # compute error result object
    if largest_blk_only:
        largest_blk = []
        for blk in a_blocklist:
            if len(blk) > len(largest_blk):
                largest_blk = blk

        a_blocklist = [largest_blk]

    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, indels)
    return err

# compute haplotype error rates between 2 VCF files
def vcf_vcf_error_rate_multiple(assembled_vcf_files, reference_vcf_files, indels, largest_blk_only=False):

    err = error_result()
    for assembled_vcf_file, reference_vcf_file in zip(assembled_vcf_files, reference_vcf_files):
        err += vcf_vcf_error_rate(assembled_vcf_file, reference_vcf_file, indels)

    return err

# compute haplotype error rates between 2 VCF files
def vcf_vcf_error_rate(assembled_vcf_file, reference_vcf_file, indels, largest_blk_only=False):

    # parse and get stuff to compute error rates
    CHROM = get_ref_name(assembled_vcf_file)
    t_blocklist = parse_vcf_phase(reference_vcf_file, CHROM, indels)
    a_blocklist = parse_vcf_phase(assembled_vcf_file, CHROM, indels)
    # compute error result object
    if largest_blk_only:
        largest_blk = []
        for blk in a_blocklist:
            if len(blk) > len(largest_blk):
                largest_blk = blk

        a_blocklist = [largest_blk]

    err = error_rate_calc(t_blocklist, a_blocklist, assembled_vcf_file, indels)
    return err

def error_rate_calc(t_blocklist, a_blocklist, vcf_file, indels=False, phase_set=None):

    ref_name    = get_ref_name(vcf_file)
    num_snps = count_SNPs(vcf_file,indels)

    switch_count   = 0
    mismatch_count = 0
    poss_sw        = 0 # count of possible positions for switch errors
    poss_mm        = 0 # count of possible positions for mismatches
    poss_flat      = 0
    flat_count     = 0
    phased_count   = 0
    maxblk_snps    = 0
    different_alleles = 0
    switch_loc     = []
    mismatch_loc   = []
    AN50_spanlst   = []
    N50_spanlst    = []

    for blk in a_blocklist:

        first_pos  = -1
        last_pos   = -1
        first_SNP  = -1
        last_SNP   = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

            #print('{}\t{}\t{}\t{}'.format(snp_ix, pos, a1, a2))
            if not (a1 == '-' or (phase_set != None and snp_ix not in phase_set)):

                phased_count += 1

                blk_phased+=1
                if first_pos == -1:
                    first_pos = pos
                    first_SNP = snp_ix
                last_pos = pos
                last_SNP = snp_ix

        blk_total = last_SNP - first_SNP + 1

        AN50_spanlst.append(((last_pos-first_pos)*(float(blk_phased)/blk_total), blk_phased))
        N50_spanlst.append((last_pos-first_pos))

        if blk_phased > maxblk_snps:
            maxblk_snps = blk_phased
    for t_block in t_blocklist:

        switched       = False
        last_base_was_switch = False

        # convert t_block to a dict for convenience
        t1_dict = defaultdict(lambda: '-')
        t2_dict = defaultdict(lambda: '-')
        a_dict = defaultdict(lambda: ('-','-','-'))
        for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in t_block:
            t1_dict[pos] = a1
            t2_dict[pos] = a2
            a_dict[pos] = (ref_str,alt1_str,alt2_str)

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.
        for a_block in a_blocklist:

            blk_switches    = [0,0]
            blk_mismatches  = [0,0]
            blk_switchlist  = [[],[]]
            blk_mmlist      = [[],[]]
            for a in [0,1]: # choose which allele to score. this only makes a difference for minimizing switch errors vs mismatches in corner cases.

                switched       = False
                last_base_was_switch = False
                first_SNP = True
                for blk_ix, (snp_ix, pos, a1, a2,  ref_str, alt1_str, alt2_str) in enumerate(a_block):
                    y = a1 if a == 0 else a2
                    x = t1_dict[pos]

                    if x == '-' or y == '-' or (phase_set != None and pos not in phase_set):
                        continue

                    #print("({},{}) == ({},{})".format(ref_str,alt_str,*a_dict[pos]))
                    if {t1_dict[pos],t2_dict[pos]} != {a1,a2} or (ref_str,alt1_str,alt2_str) != a_dict[pos]:
                        if a == 0:
                            different_alleles += 1
                        continue

                    if first_SNP:
                        switched = (x != y)
                        if count_consecutive_switches(t1_dict, a_block[blk_ix:], a) % 2 == 1:
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
                            blk_mmlist[a].append(pos)
                            blk_switches[a] -= 1      # undo count from last base switch
                            if len(blk_switchlist[a]) > 0:
                                blk_switchlist[a].pop()
                            if (blk_switches[a] < 0):
                                blk_switches[a] = 0
                            last_base_was_switch = False

                        else:

                            blk_switches[a] += 1
                            blk_switchlist[a].append(pos)
                            last_base_was_switch = True

                    else: # current base is not mismatched
                        last_base_was_switch = False

                # special case for switch on last base of previous a_block; should count as a mismatch
                if last_base_was_switch:
                    # count the switch as a single-base mismatch instead
                    blk_mismatches[a] += 1
                    blk_mmlist[a].append(pos)
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
            for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in blk:

                if {t1_dict[pos],t2_dict[pos]} != {a1,a2} or (ref_str,alt1_str,alt2_str) != a_dict[pos]:
                    continue

                if t1_dict[pos] != '-' and a1 != '-' and (phase_set == None or pos in phase_set):
                    phased_known += 1

            # a switch error is only possible in blocks len 4 or greater
            # this is because switches on the ends are counted as mismatches.
            # the -3 term: -1 because only between SNPs counts, and -2 for the two ends.
            if phased_known >= 4:
                poss_sw += (phased_known - 3)
            # a mismatch can happen in any block length 2 or more, in any position.
            if phased_known >= 2:
                poss_mm += phased_known

        # iterate over SNPs in the true and assembled haplotypes in parallel
        # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.

        for a_block in a_blocklist:

            flat_count1 = 0
            flat_count2 = 0
            #print("*******************")
            for snp_ix, pos, a1, a2, ref_str, alt1_str, alt2_str in a_block:


                if {t1_dict[pos],t2_dict[pos]} != {a1,a2} or (ref_str,alt1_str,alt2_str) != a_dict[pos]:
                    continue

                if a1 == '-' or a2 == '-' or t1_dict[pos] == '-' or (phase_set != None and pos not in phase_set):
                    continue

                if (a1 != t1_dict[pos]):
                    flat_count1 += 1
                if (a2 != t1_dict[pos]):
                    flat_count2 += 1

            if flat_count1 < flat_count2:
                flat_count += flat_count1
            else:
                flat_count += flat_count2

    if different_alleles > 0:
        print("WARNING: {} positions had different ref,alt pairs and were skipped.".format(different_alleles))

    poss_flat  = poss_mm

    if poss_sw == 0 and poss_mm == 0:
        print('WARNING: Possible switch positions and possible mismatch positions are both 0, it is likely that something is very wrong.',file=sys.stderr)

    total_error = error_result(ref=ref_name,
             switch_count=switch_count,poss_sw=poss_sw, mismatch_count=mismatch_count,
             poss_mm=poss_mm,flat_count=flat_count,poss_flat=poss_flat,
             phased_count=phased_count,num_snps=num_snps,
             maxblk_snps=maxblk_snps,
             AN50_spanlst=AN50_spanlst,N50_spanlst=N50_spanlst,switch_loc=switch_loc,
             mismatch_loc=mismatch_loc)

    return total_error

if __name__ == '__main__':

    args = parse_args()

    if (args.vcf1 == None or args.vcf2 == None):
        print("ERROR: Missing required arguments.\n--vcf1 and --vcf2 options are required", file=sys.stderr)
        sys.exit(1)

    if args.haplotype_blocks1 == None and args.haplotype_blocks2 == None:
        print(vcf_vcf_error_rate_multiple(args.vcf1, args.vcf2, args.indels))
    elif args.haplotype_blocks1 != None and args.haplotype_blocks2 == None:
        print(hapblock_vcf_error_rate_multiple(args.haplotype_blocks1, args.vcf1, args.vcf2, args.indels))
    elif args.haplotype_blocks1 != None and args.haplotype_blocks2 != None:
        print(hapblock_hapblock_error_rate_multiple(args.haplotype_blocks1, args.vcf1, args.haplotype_blocks2, args.vcf2, args.indels))
    elif args.haplotype_blocks1 == None and args.haplotype_blocks2 != None:
        print(vcf_hapblock_error_rate_multiple(args.vcf1, args.haplotype_blocks2, args.vcf2, args.indels))
