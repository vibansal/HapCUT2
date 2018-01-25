from __future__ import print_function
# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

# imports
from collections import defaultdict
import argparse
import sys

desc = '''
Calculate statistics on haplotypes assembled using HapCUT2 or similar tools.
Error rates for an assembled haplotype (specified by -h1,-v1,-f1 arguments)
are computed with respect to a "reference" haplotype (specified by -h2, -v2 arguments or -pv argument).
All files must contain information for one chromosome only (except --contig_size_file)!
To compute aggregate statistics across multiple chromosomes, provide files for
each chromosome/contig as an ordered list, using the same chromosome order between flags.
'''

def parse_args():

    parser = argparse.ArgumentParser(description=desc)
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-h1', '--haplotype_blocks', nargs='+', type = str, help='haplotype block file(s) to compute statistics on')
    parser.add_argument('-v1', '--vcf', nargs='+', type = str, help='VCF file(s) that was used to generate h1 haplotype fragments and phase h1 haplotype (--vcf in extractHAIRS and HapCUT2)')
    parser.add_argument('-f1', '--fragments', nargs='+', type = str, help='HapCUT2 format fragment file(s) used to generate input haplotype block file (-h1)')
    parser.add_argument('-pv', '--phased_vcf', nargs='*', type = str, help='compute errors with respect to this phased single-individual VCF file(s). NOTE: Files must be separated by contig/chromosome! (Use with no arguments to use same VCF(s) from --vcf.)')
    parser.add_argument('-h2', '--reference_haplotype_blocks', nargs='+', type = str, help='compute errors with respect to this haplotype block file(s)')
    parser.add_argument('-v2', '--reference_vcf', nargs='*', type = str, help='VCF file(s) that was used to generate h2 haplotype fragments and phase h2 haplotype (--vcf in extractHAIRS and HapCUT2). Use with no arguments to use same VCF(s) from --vcf.')
    parser.add_argument('-i', '--indels', action="store_true", help='Use this flag to consider indel variants. Default: Indels ignored.',default=False)
    parser.add_argument('-c', '--contig_size_file', nargs='?', type = str, help='Tab-delimited file with size of contigs (<contig>\\t<size>). If not provided, N50 will not be calculated.')

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

            a0 = el[3]
            a1 = el[4]
            a2 = None
            if ',' in a1:
                a1,a2 = a1.split(',')

            genotype = el[9].split(':')[0]
            consider = True
            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if (not indels) and (('0' in genotype and len(a0) != 1) or
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
                if chrom != CHROM:
                    print("ERROR: Chromosome in haplotype block file doesn't match chromosome in VCF")
                    print("Haplotype block file: " + hapblock_file)
                    print("VCF file:             " + vcf_file)
                    exit(1)

                pos2 = int(el[4])-1
                assert(pos == pos2) # if present, genomic position in haplotype block file should match that from VCF

            allele1 = el[1]
            allele2 = el[2]

            blocklist[-1].append((snp_ix, pos, allele1, allele2))

    return blocklist

def parse_vcf_phase(vcf_file, CHROM, indels = False):

    block = []

    with open(vcf_file, 'r') as vcf:

        snp_ix = 0

        for line in vcf:
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            if len(el) < 10:
                continue

            phase_data = el[9]

            a0 = el[3]
            a1 = el[4]
            a2 = None
            if ',' in a1:
                a1,a2 = a1.split(',')

            genotype = el[9].split(':')[0]
            consider = True

            if not (len(genotype) == 3 and genotype[0] in ['0','1','2'] and
                    genotype[1] in ['/','|'] and genotype[2] in ['0','1','2']):
                consider = False

            if genotype[0] == genotype[2]:
                consider = False

            if (not indels) and (('0' in genotype and len(a0) != 1) or
                ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                consider = False

            chrom = el[0]

            if chrom != CHROM:
                print("ERROR: Chromosome in reference haplotype VCF doesn't match chromosome in VCF used for phasing")
                print("reference haplotype VCF: " + vcf_file)

                exit(1)
            pos = int(el[1])-1
            if consider and phase_data[0:3] == '1|0' or phase_data[0:3] == '0|1':
                block.append((snp_ix, pos, phase_data[0:1], phase_data[2:3]))

            snp_ix += 1

    return [block] # we return a list containing the single block so format consistent with hapblock file format

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
                a1,a2 = a1.split(',')

            genotype = el[9][:3]

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

def count_covered_positions(frag_file, vcf_file, indels):

    covered = set()

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
                        covered.add(snp_ix)
                        # move along on qual string
            elif num_blks_newformat == int(el[0]):
                # new format

                for i in range(0,num_blks_newformat):
                    pos = int(el[2*i+5])-1 # SNPs are 1-indexed
                    seq = el[2*i+6]
                    for j, seq_base in enumerate(seq):
                        snp_ix = pos + j
                        covered.add(snp_ix)
                        # move along on qual string

    snv_ix = 0
    if not indels:
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
                    a1,a2 = a1.split(',')

                genotype = el[9][:3]

                if (('0' in genotype and len(a0) != 1) or
                    ('1' in genotype and len(a1) != 1) or ('2' in genotype and len(a2) != 1)):
                    if snv_ix in covered:
                        covered.remove(snv_ix)

                snv_ix += 1 
                        
    return len(covered)

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

    for snp_ix, pos, a1, a2 in hap:
        x = t_dict[pos]                        # base in true haplotype
        y = a1 if allele == 0 else a2   # base in assembled haplotype
        if x == '-' or y == '-':
            if first_SNP:
                continue
            else:
                break
        elif first_SNP:
            switched = (t_dict[pos] != y)
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
                 phased_count=None,num_covered=None,num_snps=None,maxblk_snps=None,
                 AN50_spanlst=None,N50_spanlst=None,switch_loc=None,mismatch_loc=None,missing_loc=None,contig_size_file=None):

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

        self.num_covered  = create_dict(num_covered, int)
        self.num_snps     = create_dict(num_snps,    int)
        self.maxblk_snps  = create_dict(maxblk_snps, int)

        self.switch_loc   = create_dict(switch_loc,   list)
        self.mismatch_loc = create_dict(mismatch_loc, list)
        self.missing_loc  = create_dict(missing_loc,  list)

        if contig_size_file == None:
            self.contig_sizes = None
        else:
            self.contig_sizes = dict()
            with open(contig_size_file,'r') as inf:
                for line in inf:
                    if len(line) < 2:
                        continue
                    el = line.strip().split('\t')
                    self.contig_sizes[el[0]] = int(el[1])

    # combine two error rate results
    def __add__(self,other):
        new_err = error_result()

        if self.contig_sizes == None and other.contig_sizes == None:
            new_err.contig_sizes = None
        elif self.contig_sizes == None and other.contig_sizes != None:
            new_err.contig_sizes = other.contig_sizes
        elif self.contig_sizes != None and other.contig_sizes == None:
            new_err.contig_sizes = self.contig_sizes
        elif self.contig_sizes != None and other.contig_sizes != None:
            assert(self.contig_sizes == other.contig_sizes)
            new_err.contig_sizes = self.contig_sizes

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
        new_err.switch_loc     = merge_dicts(self.switch_loc,     other.switch_loc)
        new_err.mismatch_loc   = merge_dicts(self.mismatch_loc,   other.mismatch_loc)
        new_err.missing_loc    = merge_dicts(self.missing_loc,    other.missing_loc)

        return new_err

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

    def get_len(self):
        return sum([self.contig_sizes[contig] for contig in self.N50_spanlst.keys()])

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

    def get_missing_rate(self):
        num_cov = self.get_num_covered()
        if num_cov > 0:
            return 1.0-sum(self.phased_count.values())/float(num_cov)
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

    def get_N50(self):
        N50  = 0
        N50_spanlst = sum(self.N50_spanlst.values(),[])
        N50_spanlst.sort(reverse=True)

        if self.contig_sizes == None:
            return 'Not calculated'

        L = self.get_len()

        total = 0
        for span in N50_spanlst:
            total += span
            if total > L/2.0:
                N50 = span
                break
        return N50

    def get_max_blk_snp_percent(self):
        snps_in_max_blks = sum(self.maxblk_snps.values())
        sum_all_snps     = self.get_num_snps()

        if sum_all_snps > 0:
            return float(snps_in_max_blks) / sum_all_snps
        else:
            return 0

    def __str__(self):

        s = ('''
switch rate:          {}
mismatch rate:        {}
flat rate:            {}
missing rate:         {}
phased count:         {}
AN50:                 {}
N50:                  {}
max block snp frac:   {}
            '''.format(self.get_switch_rate(), self.get_mismatch_rate(),
                   self.get_flat_error_rate(), self.get_missing_rate(),
                   self.get_phased_count(),
                   self.get_AN50(),self.get_N50(),self.get_max_blk_snp_percent()))

        return s

# compute error rates by using another haplotype block file as ground truth
def hapblock_hapblock_error_rate_multiple(truth_files, truth_vcf_files, assembly_files, frag_files, vcf_files, contig_size_file,indels):

    err = error_result()
    for truth_file, truth_vcf_file, assembly_file, frag_file, vcf_file in zip(truth_files, truth_vcf_files, assembly_files, frag_files, vcf_files):
        err += hapblock_hapblock_error_rate(truth_file, truth_vcf_file, assembly_file, frag_file, vcf_file, contig_size_file, indels)

    return err

# compute error rates by using another haplotype block file as ground truth
def hapblock_hapblock_error_rate(truth_file, truth_vcf_file, assembly_file, frag_file, vcf_file, contig_size_file,indels):

    # parse and get stuff to compute error rates
    t_blocklist = parse_hapblock_file(truth_file,truth_vcf_file,indels)
    a_blocklist = parse_hapblock_file(assembly_file,vcf_file,indels)
    # compute error result object
    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, contig_size_file,indels)
    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have trio phase information
def hapblock_vcf_error_rate_multiple(assembly_files, frag_files, vcf_files, phased_vcf_files, contig_size_file, indels, largest_blk_only=False):

    err = error_result()
    for assembly_file, frag_file, vcf_file, phased_vcf_file in zip(assembly_files, frag_files, vcf_files, phased_vcf_files):
        err += hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, phased_vcf_file, contig_size_file, indels)

    return err

# compute error rates by using phase data in a VCF as ground truth
# requires VCF to have phase information
def hapblock_vcf_error_rate(assembly_file, frag_file, vcf_file, phased_vcf_file, contig_size_file, indels, largest_blk_only=False):

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

    err = error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, contig_size_file, indels)
    return err

# num_covered should be the number of SNPs with coverage in the fragment matrix file
def error_rate_calc(t_blocklist, a_blocklist, vcf_file, frag_file, contig_size_file, indels=False, phase_set=None):

    ref_name    = get_ref_name(vcf_file)
    num_snps = count_SNPs(vcf_file,indels)
    num_covered = count_covered_positions(frag_file, vcf_file, indels)

    switch_count   = 0
    mismatch_count = 0
    poss_sw        = 0 # count of possible positions for switch errors
    poss_mm        = 0 # count of possible positions for mismatches
    poss_flat      = 0
    flat_count     = 0
    phased_count   = 0
    maxblk_snps    = 0
    switch_loc     = []
    mismatch_loc   = []
    missing_loc    = []
    AN50_spanlst   = []
    N50_spanlst    = []

    for blk in a_blocklist:

        first_pos  = -1
        last_pos   = -1
        first_SNP  = -1
        last_SNP   = -1
        blk_phased = 0

        for snp_ix, pos, a1, a2 in blk:

            #print('{}\t{}\t{}\t{}'.format(snp_ix, pos, a1, a2))
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

        AN50_spanlst.append(((last_pos-first_pos)*(float(blk_phased)/blk_total), blk_phased))
        N50_spanlst.append((last_pos-first_pos))

        if blk_phased > maxblk_snps:
            maxblk_snps = blk_phased
    for t_block in t_blocklist:

        switched       = False
        last_base_was_switch = False

        # convert t_block to a dict for convenience
        t_dict = defaultdict(lambda: '-')
        for snp_ix, pos, a1, a2 in t_block:
            t_dict[pos] = a1

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
                for blk_ix, (snp_ix, pos, a1, a2) in enumerate(a_block):
                    y = a1 if a == 0 else a2
                    x = t_dict[pos]

                    if x == '-' or y == '-' or (phase_set != None and pos not in phase_set):
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
            for snp_ix, pos, a1, a2 in blk:
                if t_dict[pos] != '-' and a1 != '-' and (phase_set == None or pos in phase_set):
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
            for snp_ix, pos, a1, a2 in a_block:


                if a1 == '-' or a2 == '-' or t_dict[pos] == '-' or (phase_set != None and pos not in phase_set):
                    continue

                if (a1 != t_dict[pos]):
                    flat_count1 += 1
                if (a2 != t_dict[pos]):
                    flat_count2 += 1

            if flat_count1 < flat_count2:
                flat_count += flat_count1
            else:
                flat_count += flat_count2

    poss_flat  = poss_mm

    total_error = error_result(ref=ref_name,
             switch_count=switch_count,poss_sw=poss_sw, mismatch_count=mismatch_count,
             poss_mm=poss_mm,flat_count=flat_count,poss_flat=poss_flat,
             phased_count=phased_count,num_covered=num_covered,num_snps=num_snps,
             maxblk_snps=maxblk_snps,
             AN50_spanlst=AN50_spanlst,N50_spanlst=N50_spanlst,switch_loc=switch_loc,
             mismatch_loc=mismatch_loc,missing_loc=missing_loc,contig_size_file=contig_size_file)

    return total_error

if __name__ == '__main__':

    args = parse_args()

    if (args.haplotype_blocks == None or args.vcf == None or args.fragments == None):
        print("ERROR: Missing required arguments.\n--haplotype_blocks, --vcf, and --fragments options are required", file=sys.stderr)
        sys.exit(1)

    if (args.phased_vcf == None) and (args.reference_haplotype_blocks == None or args.reference_vcf == None):
        print("ERROR: Missing reference haplotype to compute error against.\nProvide either --phased_vcf, or both --reference_haplotype_blocks and --reference_vcf", file=sys.stderr)
        sys.exit(1)

    if (args.phased_vcf != None) and (args.reference_haplotype_blocks != None or args.reference_vcf != None):
        print("ERROR: Incompatible reference haplotype arguments.\nProvide either --phased_vcf, or both --reference_haplotype_blocks and --reference_vcf", file=sys.stderr)
        sys.exit(1)

    if (args.reference_haplotype_blocks != None and args.reference_vcf == None) or (args.reference_haplotype_blocks == None and args.reference_vcf != None):
        print("ERROR: Please provide both --reference_haplotype_blocks and --reference_vcf if reference is in haplotype block format.", file=sys.stderr)
        sys.exit(1)

    if args.contig_size_file == None:
        print("WARNING: Contig size file (-c) not provided. N50 will not be calculated.", file=sys.stderr)

    reference_vcf = args.vcf if args.reference_vcf == [] else args.reference_vcf
    phased_vcf = args.vcf if args.phased_vcf == [] else args.phased_vcf

    if args.phased_vcf != None:
        print(hapblock_vcf_error_rate_multiple(args.haplotype_blocks, args.fragments, args.vcf, phased_vcf, args.contig_size_file, args.indels))
    elif args.reference_haplotype_blocks != None:
        print(hapblock_hapblock_error_rate_multiple(args.reference_haplotype_blocks, reference_vcf, args.haplotype_blocks, args.fragments, args.vcf, args.contig_size_file, args.indels))
