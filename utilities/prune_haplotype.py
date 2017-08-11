import argparse
import sys

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

    prune_hapblock_file(args.input_haplotype_blocks, args.output_haplotype_blocks, args.min_mismatch_qual, args.min_switch_qual, args.discrete_pruning)
