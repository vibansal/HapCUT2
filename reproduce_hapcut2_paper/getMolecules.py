# https://github.com/RCollins13/10XWGS/blob/master/getMolecules.py

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Tweaks by Peter Edge <pedge@eng.ucsd.edu> on 10/12/2016
# - removed a print statement that didn't do anything for python3 compatibility
# - function added to be called from python modules rather than shell

"""
Estimate original molecule sizes and coordinates from 10X linked-read WGS barcodes
"""
import argparse
from collections import defaultdict, Counter, namedtuple
from subprocess import call
import pysam

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
        #Save 10X barcode as gem
        gem = read.get_tag('RX')
	#print read.reference_name,read.reference_start,read.reference_end

        #If the read is from a new contig/chromosome, write out all molecules held in memory
        #then add read to emptied dictionary
        if read.reference_name != current_chr and current_chr is not None:
            for barcode in gemcodes:
                yield molecule(gemcodes[barcode][0].chr, min([pos for chr,
                               pos, end in gemcodes[barcode]]),
                               max([pos for chr, pos, end in gemcodes[barcode]]),
                               barcode, len(gemcodes[barcode]))

            gemcodes = defaultdict(list)

            gemcodes[gem].append(coords(read.reference_name, read.reference_start, read.reference_end))  # added read.reference_end -- pedge

        #If barcode has been seen previously and new read is colinear but beyond
        #dist, yield old barcode as interval before adding new read to list
        elif gem in gemcodes and read.reference_start - gemcodes[gem][-1].pos > dist:
            yield molecule(gemcodes[gem][0].chr,
                           min([pos for chr, pos, end in gemcodes[gem]]),
                           max([pos for chr, pos, end in gemcodes[gem]]),
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
                       max([pos for chr, pos, end in gemcodes[barcode]]),
                       barcode, len(gemcodes[barcode]))

# Peter Edge
# this is a slight f the main function meant to be called from another Python3 module
# it requires 'ref', a chromosome to output molecules for

def getMolecules(bam,outfile,ref=None,dist=50000):

    #Open outfile
    bamf = pysam.AlignmentFile(bam,"rb");
    if ref == None:
        ibam = bamf
    else:
        ibam = bamf.fetch(reference=ref)

    fout = open(outfile, 'w')

    #Get gemcode regions
    for bed in get_gemcode_regions(ibam, dist):
        #Turn molecule object into string
        bed_str = '{0}\t{1}\t{2}\t{3}\t{4}'.format(bed.chr, bed.start, bed.end,
                                                   bed.barcode, bed.readcount)

        #Write to file
        fout.write(bed_str + '\n')

    #Close outfile
    fout.close()

    #Sort outfile by chr then start (system call)
    call('sort -Vk1,1 -k2,2n -o ' + outfile + ' ' + outfile, shell=True)

#Main function
def main():
    #Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ibam', type=pysam.AlignmentFile,
                        help='Input bam')
    parser.add_argument('outfile', help='Output bed file')
    parser.add_argument('-d', '--dist', type=int, default=50000,
                        help='Molecule partitioning distance in bp (default: 50000)')
    args = parser.parse_args()

    #Open outfile
    fout = open(args.outfile, 'w')

    #Get gemcode regions
    for bed in get_gemcode_regions(args.ibam, args.dist):
        #Turn molecule object into string
        bed_str = '{0}\t{1}\t{2}\t{3}\t{4}'.format(bed.chr, bed.start, bed.end,
                                                   bed.barcode, bed.readcount)

        #Write to file
        fout.write(bed_str + '\n')

    #Close outfile
    fout.close()

    #Sort outfile by chr then start (system call)
    call('sort -Vk1,1 -k2,2n -o ' + args.outfile + ' ' + args.outfile, shell=True)


#Main block
if __name__ == '__main__':
    main()
