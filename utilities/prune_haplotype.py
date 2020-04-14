#!/usr/bin/env python3

import argparse
import sys
import os

KEEP_PRUNED=False

## fix the 'BLOCK' line of each haplotype block in the pruned haplotype file, also removes pruned variants
def fix_block_header(hapblock_file, output_file):

    with open(hapblock_file,'r') as inf, open(output_file,'w') as of:
        varlines = []
        for line in inf:
            if 'BLOCK' in line: 
                var = line.split()
                offset=-1; first = -1; last =0; phased=0; firstp=-1; lastp=0; fragments = 0;
            elif '****' in line:  ## end of block 
                if len(varlines) >= 2: 
                    print('BLOCK: offset:',offset,'len:',last-first+1,'phased:',phased,'SPAN:',lastp-firstp+1,'fragments',fragments,file=of)
                    for var in varlines: print(var,file=of,end='')
                    print(line,file=of,end='')
                varlines= []
            else: 
                var = line.split()
                if offset < 0: offset = int(var[0]);
                if first < 0: first = int(var[0]);
                if firstp < 0: firstp = int(var[4]);
                lastp  = int(var[4]);
                last = int(var[0]); 
                if var[1] != '-' and var[2] != '-':   
                    phased +=1; 
                    varlines.append(line)
                #elif KEEP_PRUNED: varlines.append(line)

        if len(varlines) >= 2: ## last block left
            print('BLOCK: offset:',offset,'len:',last-first+1,'phased:',phased,'SPAN:',lastp-firstp,'fragments',fragments,file=of)
            for var in varlines: print(var,file=of,end='')
  

def prune_hapblock_file(hapblock_file, output_file, snp_conf_cutoff, split_conf_cutoff, use_refhap_heuristic):

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


def parseargs():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_haplotype_blocks', nargs='?', type = str, help='input haplotype blocks to prune')
    parser.add_argument('-o', '--output_haplotype_blocks', nargs='?', type = str, help='output (pruned) haplotype blocks')
    parser.add_argument('-mq', '--min_mismatch_qual', nargs='?', type = float, help='minimum mismatch quality to output. default: 30', default = 30)
    parser.add_argument('-sq', '--min_switch_qual', nargs='?', type = float, help='split block at positions with switch quality below this value. default: 30', default = 30)
    parser.add_argument('-d', '--discrete_pruning', action='store_true', help='also prune positions that fail by a discrete pruning method. default: False', default = False)

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

    prune_hapblock_file(args.input_haplotype_blocks, args.output_haplotype_blocks+'.temp', args.min_mismatch_qual, args.min_switch_qual, args.discrete_pruning)
    fix_block_header(args.output_haplotype_blocks+'.temp',args.output_haplotype_blocks);
    os.remove(args.output_haplotype_blocks+'.temp');
