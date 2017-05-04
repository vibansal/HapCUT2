#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

import sys
from math import log10

def parse_hapblock_file(hapblock_file,use_SNP_index=True):

    blocklist = [] # data will be a list of blocks, each with a list tying SNP indexes to haplotypes

    try:
        with open(hapblock_file, 'r') as hbf:

            for line in hbf:
                if len(line) < 3: # empty line
                    continue
                if 'BLOCK' in line:
                    blocklist.append([])
                    continue

                elements = line.strip().split('\t')
                if len(elements) < 3: # not enough elements to have a haplotype
                    continue

                pos = int(elements[4])-1 if not use_SNP_index else int(elements[0])-1
                allele1 = elements[1]
                allele2 = elements[2]

                blocklist[-1].append((pos, allele1, allele2))

    except FileNotFoundError:
        # most of the time, this should mean that the program timed out and therefore didn't produce a phase.
        print("File {} was missing. Error calculation will continue without it.".format(hapblock_file), file=sys.stderr)
        pass

    return blocklist

def parse_vcf_phase(vcf_file,use_SNP_index=False):

    block = []

    with open(vcf_file, 'r') as vcf:
        i = 0
        for line in vcf:
            if line[0] == '#':
                continue

            elements = line.strip().split('\t')
            if len(elements) < 10:
                continue

            phase_data = elements[9]
            pos = int(elements[1])-1 if not use_SNP_index else i
            if phase_data[0:3] == '1|0' or phase_data[0:3] == '0|1':
                block.append((pos, phase_data[0:1], phase_data[2:3]))

            i += 1

    return [block] # we return a list containing the single block so format consistent with hapblock file format

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

def count_frags(frag_file):
    count = 0
    with open(frag_file, 'r') as infile:
        for line in infile:
            if len(line.split()) >= 5:
                count += 1
    return count


def parse_runtime_file(runtime_file):
    with open(runtime_file, 'r') as rf:
        return float(rf.readline().strip())

def prune_hapblock_file(hapblock_file, output_file, snp_conf_cutoff, split_conf_cutoff, use_refhap_heuristic):

    snp_conf_cutoff = -10*log10(1-snp_conf_cutoff)
    split_conf_cutoff = -10*log10(1-split_conf_cutoff)

    with open(hapblock_file,'r') as inf, open(output_file,'w') as of:
        blk_count = 0
        for line in inf:
            if 'BLOCK' in line:
                blk_count = 0
            if len(line) < 3 or 'BLOCK' in line or '****' in line:
                print(line,file=of,end='')
                continue

            el = line.strip().split()
            pruned_refhap_heuristic = int(el[8])
            split_conf = float(el[9]) if el[9] != '.' else 100
            snp_conf   = float(el[10]) if el[10] != '.' else 100

            if split_conf < split_conf_cutoff and blk_count >= 2:
                print('******** ',file=of)
                print('BLOCK: (from split)',file=of)

            if (use_refhap_heuristic and pruned_refhap_heuristic) or (snp_conf < snp_conf_cutoff):
                el[1] = '-'
                el[2] = '-'
                print('\t'.join(el),file=of)
            else:
                print(line,file=of,end='')

            blk_count += 1

# convert the new fragment file format to the old format
# assumes normal fragments (no special modeling)

def prune_probhap_file(hapblock_file, output_file, emission_cutoff, split_cutoff):

    with open(hapblock_file,'r') as inf, open(output_file,'w') as of:
        blk_count = 0
        had_posterior_split = False
        for line in inf:
            if 'BLOCK' in line:
                blk_count = 0
                had_posterior_split = False
            if len(line) < 3 or 'BLOCK' in line or '****' in line:
                print(line,file=of,end='')
                continue

            el = line.strip().split()
            transition_prob = float(el[3])
            posterior_prob  = float(el[4])
            emission_prob   = float(el[5])

            if (((transition_prob < split_cutoff) or
                (posterior_prob < split_cutoff and not had_posterior_split))
                and blk_count >= 2):
                had_posterior_split = True
                print('******** ',file=of)
                print('BLOCK: (from split)',file=of)

            if emission_prob < emission_cutoff:
                el[1] = '-'
                el[2] = '-'
                print('\t'.join(el),file=of)
            else:
                print(line,file=of,end='')

            blk_count += 1

def new_to_old_format(inputfile,outputfile):

    with open(inputfile,'r') as i, open(outputfile,'w') as o:
        for line in i:
            if len(line) < 3:
                continue

            el = line.strip().split()
            edited_el = el[0:2] + el[5:]
            new_line = ' '.join(edited_el)
            print(new_line,file=o)

# convert the old fragment file format to the new format
# assumes normal fragments (no special modeling)
def old_to_new_format(inputfile,outputfile):

    with open(inputfile,'r') as i, open(outputfile,'w') as o:
        for line in i:
            if len(line) < 3:
                continue
            el = line.strip().split()
            new_fields = ['0','-1','-1']
            edited_el = el[0:2] + new_fields + el[2:]
            new_line = ' '.join(edited_el)
            print(new_line,file=o)
