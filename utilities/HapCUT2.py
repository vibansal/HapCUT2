import multiprocessing
import os
import time
import argparse
import sys
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

SPACER = "---------------------------------------------------------------------------------"

# parse HapCUT2 program input.
def parseargs():

    desc = '''
    HapCUT2 is a maximum-likelihood-based tool for assembling haplotypes from DNA sequence reads, designed to "just work" with excellent speed and accuracy.
    '''

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-v', '--vcf', nargs='?', type = str, help='VCF file with variants to phase',required=True)
    parser.add_argument('-b', '--bams', nargs='+', type = str, help='generic bamfile(s)', default=[])
    parser.add_argument('-o', '--out', nargs='?', type = str, help='output VCF file',required=True)
    parser.add_argument('-p', '--processes', nargs='?', type = str, help='number of processes to use', default=4)
    parser.add_argument('-il', '--illumina_bams', nargs='+', type = str, help='illumina bamfile(s)', default=[])
    parser.add_argument('-hic', '--hic_bams', nargs='+', type = str, help='Hi-C bamfile(s)', default=[])
    parser.add_argument('-tenX', '--tenX_bams', nargs='+', type = str, help='10x genomics bamfile(s)', default=[])
    parser.add_argument('-pb', '--pacbio_bams', nargs='+', type = str, help='PacBio bamfile(s)', default=[])
    parser.add_argument('-ont', '--ont_bams', nargs='+', type = str, help='Oxford Nanopore bamfile(s)', default=[])
    parser.add_argument('--ref', nargs='?', type = str, help='reference sequence file (in fasta format, gzipped is okay), optional but required for indels, should be indexed',default=0)
    parser.add_argument('--indels', action='store_true', help='phase indel variants')
    parser.add_argument('--converge', nargs='?', type = int, help='cut off iterations (global or maxcut) after this many iterations with no improvement.',default=5)
    parser.add_argument('--verbose', action='store_true', help='verbose mode: print extra information to stdout and stderr.')

    parser.add_argument('--qvoffset', nargs='?', type = int, help='quality value offset, 33/64 depending on how quality values were encoded', default=33)
    parser.add_argument('--mbq', nargs='?', type = int, help='minimum base quality to consider a base for haplotype fragment',default=13)
    parser.add_argument('--maxIS', nargs='?', type = int, help='maximum insert size for a paired-end read to be considered as a single fragment for phasing',default=1000)
    parser.add_argument('--minIS', nargs='?', type = int, help='minimum insert size for a paired-end read to be considered as a single fragment for phasing',default=0)
    parser.add_argument('--tenX_distance', nargs='?', type = int, help='distance in base pairs that delineates separate 10X molecules',default=20000)
    parser.add_argument('--PEonly', action='store_true', help='do not use single end reads')
    parser.add_argument('--noquality', nargs='?', type = int, help='if the bam file does not have quality string, this value will be used as the uniform quality value',default=0)
    parser.add_argument('--ep', action='store_true', help='set this flag to estimate HMM parameters from aligned reads (only with long reads)')


    # optional nicknames for the genomes that the samfiles map to (recommended)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# run a shell command
# prints the PID and the command before running
# pid: the PID of the current process
# cmd: the command to run
def run_command(cmd):
    pid = os.getpid()

    prog = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #except subprocess.CalledProcessError as e:
    #    print("ERROR: PID {} returned non-zero return code {}".format(pid, e.returncode), file=sys.stderr)
    #    print(e.output,file=sys.stderr)

    out, err = prog.communicate()
    o = out.decode("ISO-8859-1")
    e = err.decode("ISO-8859-1")

    print(SPACER)
    print("[PID {}] {}".format(pid, cmd))
    for line in o.split('\n'):
        print("[PID {}] {}".format(pid, line),file=sys.stdout)
    for line in e.split('\n'):
        print("[PID {}] {}".format(pid, line),file=sys.stderr)

    if prog.returncode != 0:

        #print(SPACER)
        print("ERROR: Job on PID {} ended with non-zero return code {}".format(pid,prog.returncode))
        exit(1)

    return o+'\n'+e

# build a string representing an extractHAIRS command
# args: the argparse object holding the command line options
# chrom: the chromosome to extract fragments for
# bam: the bamfile to extract fragments for
# out: the fragment file to write to
# mode_str: a string with sequencing-technology-specific flags
def build_extractHAIRS_command(args, chrom, bam, out, mode_str):
    cmd = ("./build/extractHAIRS {} --region {} --bam {} --vcf {} --out {} --ref {} --indels {} --qvoffset {} --mbq {} --maxIS {} --minIS {} --PEonly {} --noquality {} --ep {}"
    .format(mode_str, chrom, bam, args.vcf, out, args.ref, int(args.indels), args.qvoffset, args.mbq, args.maxIS,
    args.minIS, int(args.PEonly), args.noquality, int(args.ep)))
    return cmd

# main function
# parse command line arguments,
# starts up a multiprocessing queue,
# passes off extractHAIRS jobs to the queue,
# then merges the extractHAIRS fragment files per chromosome
# and phases them with HAPCUT2
def main():

    # parse command line arguments
    args = parseargs()

    # check that at least one BAM file was provided
    bams_len = (len(args.bams) + len(args.illumina_bams) + len(args.hic_bams)
               + len(args.tenX_bams) + len(args.pacbio_bams) + len(args.ont_bams))
    if bams_len == 0:
        print("Error: at least one BAM file must be specified",file=sys.stderr)
        sys.exit(1)

    # build the ordered list of chromosomes in the source VCF
    chromset = set()
    chroms = []
    with open(args.vcf, 'r') as infile:
        for line in infile:
            if line[0] == '#':
                continue
            chrom = line.strip().split('\t')[0]
            if chrom not in chromset:
                chromset.add(chrom)
                chroms.append(chrom)

    task_list = []
    # frag_files[chrom] contains the list of fragment files for that chomosome
    frag_files = defaultdict(list)
    # phased_vcf_files[chrom] contains the HapCUT2 phased VCF file for that chromosome
    phased_vcf_files = dict()
    combined_fragment_files = dict()

    files_to_remove = []

    def process_tasks(cmd_lst):
        futures = []
        with ProcessPoolExecutor(max_workers=args.processes) as executor:
            for cmd in cmd_lst:
                futures.append(executor.submit(run_command, cmd))

        job_failed = False
        for f in futures:
            if f.exception() != None:
                job_failed = True

        if job_failed:
            print("Exception encountered in subprocess. Exiting...", file=sys.stderr)
            # clean up files
            for file in files_to_remove:
                os.remove(file)
            for file in phased_vcf_files.values():
                os.remove(file)
            for file in combined_fragment_files.values():
                os.remove(file)
            for file_lst in frag_files.values():
                for file in file_lst:
                    os.remove(file)
            exit(1)


    # for each chromosome and each bam,
    # use extractHAIRS to extract haplotype fragments into temp files
    for chrom in chroms:

        # generic BAMs
        for bam in args.bams:
            (handle, frag_filename) = tempfile.mkstemp()
            os.close(handle)
            frag_files[chrom].append(frag_filename)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename, "--nf 1");
            task_list.append(cmd)

        # illumina BAMs
        for bam in args.illumina_bams:
            (handle, frag_filename) = tempfile.mkstemp()
            os.close(handle)
            frag_files[chrom].append(frag_filename)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename, "--nf 1");
            task_list.append(cmd)

        # 10X genomics BAMs
        for bam in args.tenX_bams:
            (handle, frag_filename_unlinked) = tempfile.mkstemp()
            os.close(handle)
            (handle, frag_filename_linked) = tempfile.mkstemp()
            os.close(handle)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename_unlinked, "--10X 1");

            cmd += "; python3 LinkFragments.py -f {} -v {} -b {} -d {} -o {}".format(frag_filename_unlinked, args.vcf, bam, args.tenX_distance, frag_filename_linked)
            files_to_remove.append(frag_filename_unlinked)

            frag_files[chrom].append(frag_filename_linked)

            task_list.append(cmd)

        # HiC BAMs
        for bam in args.hic_bams:
            (handle, frag_filename) = tempfile.mkstemp()
            os.close(handle)
            frag_files[chrom].append(frag_filename)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename, "--hic 1");
            task_list.append(cmd)

        # PacBio BAMs
        for bam in args.pacbio_bams:
            (handle, frag_filename) = tempfile.mkstemp()
            os.close(handle)
            frag_files[chrom].append(frag_filename)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename, "--pacbio 1 --nf 1");
            task_list.append(cmd)

        # Oxford Nanopore BAMs
        for bam in args.ont_bams:
            (handle, frag_filename) = tempfile.mkstemp()
            os.close(handle)
            frag_files[chrom].append(frag_filename)
            cmd = build_extractHAIRS_command(args, chrom, bam, frag_filename, "--ont 1 --nf 1");
            task_list.append(cmd)

    # process extractHAIRS jobs
    process_tasks(task_list)
    task_list.clear()

    # for each chromosome,
    # merge the fragment files for that chromosome into one big file
    for chrom in chroms:

        # combine the fragment files into one big fragment file
        (handle, combined_frag_filename) = tempfile.mkstemp()
        os.close(handle)
        cmd = "cat " + " ".join(frag_files[chrom]) + " > " + combined_frag_filename
        combined_fragment_files[chrom] = combined_frag_filename
        task_list.append(cmd)

    # process cat jobs
    process_tasks(task_list)
    task_list.clear()

    # for each chromosome,
    # phase the fragment files with HapCUT2. Write the output to the final output VCF.
    for chrom in chroms:

        # temp file to write HapCUT2 output to
        (handle, out_vcf_filename) = tempfile.mkstemp()
        os.close(handle)
        phased_vcf_files[chrom] = out_vcf_filename

        # parameters
        long_reads = int(len(args.pacbio_bams)+len(args.ont_bams)+len(args.tenX_bams) > 1)
        hic = int(len(args.hic_bams) > 1)

        # run HAPCUT2
        cmd = ("./build/HAPCUT2 --f {} --vcf {} --out {} --only_chrom {} --helper_script_mode 1 --converge {} --verbose {} --qo {} --long_reads {} --hic {} --nf 1"
               .format(combined_fragment_files[chrom], args.vcf, out_vcf_filename, chrom, args.converge, int(args.verbose), args.qvoffset, long_reads, hic))

        task_list.append(cmd)

    # process HapCUT2 jobs
    process_tasks(task_list)
    task_list.clear()

    print(SPACER)
    print("Combining VCFs from each chromosome into a single VCF...")
    # combine the HapCUT2 phased VCFs into one big VCF
    with open(args.out, 'w') as outfile:
        for chrom in chroms:
            with open(phased_vcf_files[chrom], 'r') as inf:
                for line in inf:
                    print(line.strip(),file=outfile)

    # clean up files
    for file in files_to_remove:
        os.remove(file)
    for file in phased_vcf_files.values():
        os.remove(file)
    for file in combined_fragment_files.values():
        os.remove(file)
    for file_lst in frag_files.values():
        for file in file_lst:
            os.remove(file)

    print("Done! Phased VCF written to: {}".format(args.out))


if __name__ == "__main__":
    main()
