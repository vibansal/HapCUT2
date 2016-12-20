# This file contains wrapper functions to run various haplotype assembly tools

import time
import sys
import os
import fileIO
import subprocess
import shutil

# Runs (cmd) as a process, kills it after (timeout) seconds
# kills process if memory hits 8 GB
def run_process(cmd, timeout=None, memlimit=7700000):
    cmd += ' 2>&1'
    print(cmd)
    # credit to roland smith for method of retrieving stdout and stderr: http://stackoverflow.com/questions/14059558/why-is-python-no-longer-waiting-for-os-system-to-finish
    t1 = time.time()
    if timeout == None:
        cmd = "ulimit -v {}; {}".format(memlimit, cmd)
        prog = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    else:
        timeout_cmd = 'ulimit -v {}; timeout {} {}'.format(memlimit, timeout, cmd)
        prog = subprocess.Popen(timeout_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = prog.communicate()
    t2 = time.time()
    runtime = t2 -t1
    o = out.decode("ISO-8859-1")
    print(o)
    return runtime

# a wrapper function for calling the HapCUT2 program
def run_hapcut2(path_to_bin, frag_file, vcf_file, output_file, converge, threshold, hapcut_xtra_args,timeout=None):
    cmd = "{} --fragments {} --vcf {} --output {} --converge {} --t {} {}".format(path_to_bin, frag_file, vcf_file, output_file, converge, threshold, hapcut_xtra_args)
    return run_process(cmd,timeout)
'''
def run_hapcut2(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads, MEC, converge, threshold, hapcut_xtra_args,timeout=None):
    cmd = "{} --fragments {} --vcf {} --output {} --maxiter {} --maxcutiter {} --longreads {} --maxmem 32000 --MEC {} --converge {} --t {} {}".format(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads, MEC, converge, threshold, hapcut_xtra_args)
    return run_process(cmd,timeout)
'''
# a wrapper function for calling the HapCUT program
def run_hapcut(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads=0, timeout=None):
    cmd = "{} --fragments {} --VCF {} --output {} --maxiter {} --maxcutiter {} --longreads {} --maxmem 32000".format(path_to_bin, frag_file, vcf_file, output_file, maxiter, maxcutiter, longreads)
    return run_process(cmd,timeout)

def run_refhap(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

def run_dgs(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH -a DGS {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

def run_fasthare(path_to_jar, frag_file, output_file,timeout=None):
    cmd = "java -cp {} mpg.molgen.sih.main.SIH -a FastHare {} {}".format(path_to_jar, frag_file, output_file)
    return run_process(cmd,timeout)

# unfortunately probhap is implemented in python2, so it has to be run in the same
# ugly fashion as HapCUT and RefHAP with a system call
def run_probhap(path_to_probhap, path_to_python2, frag_file, vcf_file, output_file, probhap_frag_file,timeout=None):

    num_reads = fileIO.count_frags(frag_file)
    num_snps  = fileIO.count_SNPs(vcf_file)

    # probhap requires top line in fragment file to specify matrix dimensions, so we need to edit the file...

    with open (probhap_frag_file, 'w') as pff:
        print('{} {}'.format(num_reads, num_snps), file=pff)
        with open (frag_file, 'r') as ff:
            for line in ff:
                print(line.strip(), file=pff)

    parsed_reads = "{}.parsed_fragments".format(output_file)
    phase = "{}.uncorrected".format(output_file)
    assignments = "{}.assignments".format(output_file)
    cmd = "{0} {1} --reads {2} --parsed-reads {3} --phase {4} --assignments {5}; python2 {6}-postprocess.py --filtered-reads {3} --assignments {5} --blocks {4} --corrected-blocks {7}".format(path_to_python2, path_to_probhap, probhap_frag_file, parsed_reads, phase, assignments, path_to_probhap[:-3], output_file)

    return run_process(cmd,timeout)

# a wrapper function for calling the MixSIH program
def run_mixsih(path_to_bin, frag_file, vcf_file, output_file, mixsih_frag_file,timeout=None):

    num_snps  = fileIO.count_SNPs(vcf_file)

    # mixsih fragment file has no qual scores in the last column
    # it also has the num SNPs at top
    # we will simultaneously estimate average miscall rate and remove the qual scores
    total_qual = 0
    num_qual   = 0
    mixsih_frag_file_unsorted = mixsih_frag_file + ".unsorted"
    with open (mixsih_frag_file_unsorted, 'w') as mff:
        print(str(num_snps), file=mff)
        with open(frag_file, 'r') as ff:
            for line in ff:
                el = line.strip().split()
                if len(el) > 2:
                    # tally up the present qual scores
                    for q in el[-1]:
                        num_qual += 1
                        total_qual += 10**((ord(q)-33)/-10)

                    # remove qual string
                    el = el[:-1]
                    print(' '.join(el), file=mff)
                else:
                    print(line.strip(), file=mff)

    # mixsih requires we input a miscall rate, so we take the average
    miscall_rate = total_qual/num_qual
    # need to be sorted by value of third column...
    os.system("sort -n -k 3 {} > {}".format(mixsih_frag_file_unsorted, mixsih_frag_file))

    cmd = "{0} -a {1} {2} {3}.profile {3}".format(path_to_bin, miscall_rate, mixsih_frag_file, output_file)
    return run_process(cmd,timeout)

# a wrapper function for calling the Haptree program
def run_haptree(path_to_bin, frag_file, vcf_file, output_file, haptree_vcf_file,timeout=None):

    # haptree vcf requires no header
    with open (haptree_vcf_file, 'w') as ht_vcf:
        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line[:1] != '#':
                    print(line.strip(), file=ht_vcf)

    haptree_output_dir = output_file + ".dir"

    cmd = "{} {} {} {}".format(path_to_bin, frag_file, haptree_vcf_file, haptree_output_dir)
    haptree_solution = os.path.join(haptree_output_dir,"HapTreeSolution")
    if os.path.isfile(haptree_solution):
        shutil.move(haptree_solution, output_file)
        os.rmdir(haptree_output_dir)

    return run_process(cmd,timeout)
