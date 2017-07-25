
import pysam
from LinkFragments import fragment, link_fragments
import filecmp
import os
import unittest
import shutil

# variant_locations should be a list of (chr, pos, het) tuples
# "het" is true if heteroygous else false

read_len = 100
# for simplicity, at generated heteroygous positions we just assume the ref is example_a1 and variant is example_a2
example_a1 = 'A'
example_a2 = 'T'

chroms = ['chr1','chr2','chr3']
chrom_to_id = {'chr1':0,'chr2':1,'chr3':2}
EXTRACTHAIRS = '../build/extractHAIRS'

def link_fragments_test_helper(output_dir, variant_locations, read_infos, single_reads=True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # generate VCF file
    vcf_filename = os.path.join(output_dir,'test.vcf')
    with open(vcf_filename,'w') as vcf_out:
        for chrom,pos,het in variant_locations:
            genotype = '0/1' if het else '0/0'
            if het:
                a1,a2 = example_a1,example_a2
            else:
                a1,a2 = 'N','N'

            line = '{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}'.format(chrom,pos,a1,a2,genotype)
            print(line,file=vcf_out)

    frag_file = os.path.join(output_dir,'test.vcf')

    # generate bam file
    header = dict()
    header['SQ'] = []
    for chrom in chroms:
        header['SQ'].append({'SN':chrom, 'LN':1000000})

    bam_filename = os.path.join(output_dir,'test.bam')
    bam_file = pysam.AlignmentFile(bam_filename,'wb',header=header)

    # variants is simply a list of indices into variant_locations, if it has variant at that loc
    for read_id, read_chrom, read_start1, barcode, flag, variants, mate_info in read_infos:

        read_start = read_start1 - 1 # 0 index
        read_seqlst = ['N']*read_len

        for var_chrom, var_pos1, var_het in variant_locations:
            var_pos = var_pos1 - 1 # 0 index
            if var_chrom == read_chrom and var_pos >= read_start and var_pos < read_start + read_len:
                read_seqlst[var_pos - read_start] = example_a1

        for var_ix1 in variants:
            var_ix = var_ix1 - 1
            var_chrom, var_pos1, var_het = variant_locations[var_ix]
            var_pos = var_pos1 - 1 # 0 index
            assert(var_chrom == read_chrom and var_pos >= read_start and var_pos < read_start + read_len)
            read_seqlst[var_pos - read_start] = example_a2
                # variant is in this read

        ar = pysam.AlignedSegment()
        ar.query_sequence = ''.join(read_seqlst)
        ar.reference_id = chrom_to_id[read_chrom]
        ar.qname = read_id
        ar.flag = flag # need to do paired-end later...
        ar.reference_start = read_start
        ar.mapping_quality = 60
        ar.cigartuples = [(0,read_len)]
        ar.qual = ''.join(['5']*read_len)
        ar.set_tag('BX',barcode)

        if mate_info != None:

            mate_ix, tlen = mate_info
            ar.next_reference_start = read_infos[mate_ix][2] - 1
            ar.next_reference_id    = chrom_to_id[read_infos[mate_ix][1]]
            ar.template_length = tlen

        bam_file.write(ar)

    bam_file.close()
    pysam.index(bam_filename)

    unlinked_hairs_filename = os.path.join(output_dir,'unlinked.frag')
    cmd = '{} --vcf {} --bam {} --out {} --10x 1'.format(EXTRACTHAIRS, vcf_filename, bam_filename, unlinked_hairs_filename)
    print(cmd)
    os.system(cmd)

    linked_hairs_filename = os.path.join(output_dir,'linked.frag')
    link_fragments(unlinked_hairs_filename, vcf_filename, bam_filename, linked_hairs_filename, 20000, single_reads)

    lines = []
    with open(linked_hairs_filename,'r') as inf:

        for line in inf:
            line = line.strip()
            lines.append(line)

    return lines

class test_basic_functionality(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_basic_functionality/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr1',153000,1),
                       ('chr1',157000,1),
                       ('chr1',160000,1),
                       ('chr1',165000,1),
                       ('chr1',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr1',152950,'AAA',0,(),None],
                    ['r4','chr1',156950,'AAA',0,(),None],
                    ['r5','chr1',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr1',165000,'AAA',0,(8,),None], # variant on the first position of read
                    ['r7','chr1',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 chr1:152950-170100:AAA 0 -1 -1 5 00111 55555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_basic_functionality/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr2',153000,1),
                       ('chr2',157000,1),
                       ('chr3',160000,1),
                       ('chr3',165000,1),
                       ('chr3',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr2',152950,'AAA',0,(),None],
                    ['r4','chr2',156950,'AAA',0,(),None],
                    ['r5','chr3',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr3',165000,'AAA',0,(),None], # variant on the first position of read
                    ['r7','chr3',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 chr2:152950-157050:AAA 0 -1 -1 5 00 55',
                          '1 chr3:159901-170100:AAA 0 -1 -1 7 101 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

class test_singleton_fragments(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_singleton_fragments/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',53000,1),
                       ('chr1',107000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',52950,'GGG',0,(2,),None],
                    ['r3','chr1',106950,'AAA',0,(3,),None]]

        expected_lines = ['1 chr1:1950-2050:GGG 0 -1 -1 1 1 5',
                          '1 chr1:52950-53050:GGG 0 -1 -1 2 1 5',
                          '1 chr1:106950-107050:AAA 0 -1 -1 3 1 5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_singleton_fragments/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr2',53000,1),
                       ('chr2',74000,1),
                       ('chr3',1000,1),
                       ('chr3',11000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(),None],
                    ['r2','chr2',53000,'GGG',0,(2,),None],
                    ['r3','chr2',74000,'GGG',0,(),None],
                    ['r4','chr3',950,'AAA',0,(4,),None],
                    ['r5','chr3',2000,'AAT',0,(),None],
                    ['r6','chr3',11000,'AAA',0,(),None]]

        expected_lines = ['1 chr2:53000-53100:GGG 0 -1 -1 2 1 5',
                          '1 chr2:74000-74100:GGG 0 -1 -1 3 0 5',
                          '1 chr3:950-11100:AAA 0 -1 -1 4 10 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

class test_overlapping_fragments(unittest.TestCase):

    def setUp(self):
        pass

    # overlapping reads agree on variant
    def test1(self):
        out_dir = 'tests/test_overlapping_fragments/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 111 5I5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads agree on reference
    def test2(self):
        out_dir = 'tests/test_overlapping_fragments/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 101 5I5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant (call ignored)
    def test3(self):
        out_dir = 'tests/test_overlapping_fragments/test3'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant, 3 reads (call ignored)
    def test4(self):
        out_dir = 'tests/test_overlapping_fragments/test4'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',4950,'GGG',0,(2,),None],
                    ['r5','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant, 3 reads different read loc (call ignored)
    def test5(self):
        out_dir = 'tests/test_overlapping_fragments/test5'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4960,'GGG',0,(),None],
                    ['r4','chr1',4970,'GGG',0,(2,),None],
                    ['r5','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # test that the quality caps out at a valid character
    def test6(self):
        out_dir = 'tests/test_overlapping_fragments/test6'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4950,'GGG',0,(2,),None],
                    ['r5','chr1',4950,'GGG',0,(2,),None],
                    ['r6','chr1',4950,'GGG',0,(2,),None],
                    ['r7','chr1',4950,'GGG',0,(2,),None],
                    ['r8','chr1',4950,'GGG',0,(2,),None],
                    ['r9','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 111 5~5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)


    # test multiple overlaps in the same reads
    def test7(self):
        out_dir = 'tests/test_overlapping_fragments/test7'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',5010,1),
                       ('chr1',5020,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4960,'AAA',0,(),None],
                    ['r5','chr1',9950,'GGG',0,(5,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 11001 5III5',
                          '1 chr1:4960-5060:AAA 0 -1 -1 2 000 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # test multiple overlaps (and disagreements) in the same reads
    def test8(self):
        out_dir = 'tests/test_overlapping_fragments/test8'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',5010,1),
                       ('chr1',5020,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,3),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4960,'AAA',0,(),None],
                    ['r5','chr1',9950,'GGG',0,(5,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 11 4 01 5II5',
                          '1 chr1:4960-5060:AAA 0 -1 -1 2 000 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # test multiple overlaps (and disagreements) in the same reads, multiple chromosomes present
    def test8(self):
        out_dir = 'tests/test_overlapping_fragments/test8'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr2',5000,1),
                       ('chr2',5010,1),
                       ('chr2',5020,1),
                       ('chr2',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr2',4950,'GGG',0,(2,3),None],
                    ['r3','chr2',4950,'GGG',0,(2,),None],
                    ['r4','chr2',4960,'AAA',0,(),None],
                    ['r5','chr2',5010,'AAA',0,(3,),None],
                    ['r6','chr3',9950,'GGG',0,(),None]]

        expected_lines = ['1 chr1:1950-2050:GGG 0 -1 -1 1 1 5',
                          '2 chr2:4950-5050:GGG 0 -1 -1 2 1 4 0 II',
                          '2 chr2:4960-5110:AAA 0 -1 -1 2 0 4 0 5I']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

class test_basic_functionality_no_single_reads(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_basic_functionality_no_single_reads/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr1',153000,1),
                       ('chr1',157000,1),
                       ('chr1',160000,1),
                       ('chr1',165000,1),
                       ('chr1',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr1',152950,'AAA',0,(),None],
                    ['r4','chr1',156950,'AAA',0,(),None],
                    ['r5','chr1',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr1',165000,'AAA',0,(8,),None], # variant on the first position of read
                    ['r7','chr1',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 chr1:152950-170100:AAA 0 -1 -1 5 00111 55555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_basic_functionality_no_single_reads/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr2',153000,1),
                       ('chr2',157000,1),
                       ('chr3',160000,1),
                       ('chr3',165000,1),
                       ('chr3',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr2',152950,'AAA',0,(),None],
                    ['r4','chr2',156950,'AAA',0,(),None],
                    ['r5','chr3',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr3',165000,'AAA',0,(),None], # variant on the first position of read
                    ['r7','chr3',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 chr2:152950-157050:AAA 0 -1 -1 5 00 55',
                          '1 chr3:159901-170100:AAA 0 -1 -1 7 101 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

class test_singleton_fragments_no_single_reads(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_singleton_fragments_no_single_reads/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',53000,1),
                       ('chr1',107000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',52950,'GGG',0,(2,),None],
                    ['r3','chr1',106950,'AAA',0,(3,),None]]

        expected_lines = []

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_singleton_fragments_no_single_reads/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr2',53000,1),
                       ('chr2',74000,1),
                       ('chr3',1000,1),
                       ('chr3',11000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(),None],
                    ['r2','chr2',53000,'GGG',0,(2,),None],
                    ['r3','chr2',74000,'GGG',0,(),None],
                    ['r4','chr3',950,'AAA',0,(4,),None],
                    ['r5','chr3',2000,'AAT',0,(),None],
                    ['r6','chr3',11000,'AAA',0,(),None]]

        expected_lines = ['1 chr3:950-11100:AAA 0 -1 -1 4 10 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

class test_overlapping_fragments_no_single_reads(unittest.TestCase):

    def setUp(self):
        pass

    # overlapping reads agree on variant
    def test1(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 111 5I5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads agree on reference
    def test2(self):
        out_dir = 'tests/test_overlapping_fragments/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 101 5I5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant (call ignored)
    def test3(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test3'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant, 3 reads (call ignored)
    def test4(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test4'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(),None],
                    ['r4','chr1',4950,'GGG',0,(2,),None],
                    ['r5','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # overlapping reads disagree ref/variant, 3 reads different read loc (call ignored)
    def test5(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test5'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4960,'GGG',0,(),None],
                    ['r4','chr1',4970,'GGG',0,(2,),None],
                    ['r5','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 1 3 1 55']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # test that the quality caps out at a valid character
    def test6(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test6'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4950,'GGG',0,(2,),None],
                    ['r5','chr1',4950,'GGG',0,(2,),None],
                    ['r6','chr1',4950,'GGG',0,(2,),None],
                    ['r7','chr1',4950,'GGG',0,(2,),None],
                    ['r8','chr1',4950,'GGG',0,(2,),None],
                    ['r9','chr1',9950,'GGG',0,(3,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 111 5~5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)


    # test multiple overlaps in the same reads
    def test7(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test7'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',5010,1),
                       ('chr1',5020,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4960,'AAA',0,(),None],
                    ['r5','chr1',9950,'GGG',0,(5,),None]]

        expected_lines = ['1 chr1:1950-10050:GGG 0 -1 -1 1 11001 5III5',
                          '1 chr1:4960-5060:AAA 0 -1 -1 2 000 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # test multiple overlaps (and disagreements) in the same reads
    def test8(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test8'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',5000,1),
                       ('chr1',5010,1),
                       ('chr1',5020,1),
                       ('chr1',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr1',4950,'GGG',0,(2,3),None],
                    ['r3','chr1',4950,'GGG',0,(2,),None],
                    ['r4','chr1',4960,'AAA',0,(),None],
                    ['r5','chr1',9950,'GGG',0,(5,),None]]

        expected_lines = ['2 chr1:1950-10050:GGG 0 -1 -1 1 11 4 01 5II5',
                          '1 chr1:4960-5060:AAA 0 -1 -1 2 000 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)

    # test multiple overlaps (and disagreements) in the same reads, multiple chromosomes present
    def test8(self):
        out_dir = 'tests/test_overlapping_fragments_no_single_reads/test8'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr2',5000,1),
                       ('chr2',5010,1),
                       ('chr2',5020,1),
                       ('chr2',10000,1)]

        read_info = [['r1','chr1',1950,'GGG',0,(1,),None],
                    ['r2','chr2',4950,'GGG',0,(2,3),None],
                    ['r3','chr2',4950,'GGG',0,(2,),None],
                    ['r4','chr2',4960,'AAA',0,(),None],
                    ['r5','chr2',5010,'AAA',0,(3,),None],
                    ['r6','chr3',9950,'GGG',0,(),None]]

        expected_lines = ['2 chr2:4950-5050:GGG 0 -1 -1 2 1 4 0 II',
                          '2 chr2:4960-5110:AAA 0 -1 -1 2 0 4 0 5I']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info,single_reads=False)
        self.assertEqual(output_lines,expected_lines)


class test_no_barcode(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_no_barcode/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr1',153000,1),
                       ('chr1',157000,1),
                       ('chr1',160000,1),
                       ('chr1',165000,1),
                       ('chr1',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr1',152950,None,0,(),None],
                    ['r4','chr1',156950,None,0,(),None],
                    ['r5','chr1',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr1',165000,'AAA',0,(8,),None], # variant on the first position of read
                    ['r7','chr1',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 r3 0 -1 -1 5 0 5',
                          '1 r4 0 -1 -1 6 0 5',
                          '1 chr1:159901-170100:AAA 0 -1 -1 7 111 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_no_barcode/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,0),
                       ('chr1',3000,1),
                       ('chr1',7000,1),
                       ('chr1',150000,0),
                       ('chr1',153000,1),
                       ('chr1',153010,1),
                       ('chr1',160000,1),
                       ('chr1',165000,1),
                       ('chr1',170000,1)]

        read_info = [['r1','chr1',2950,'GGG',0,(2,),None],
                    ['r2','chr1',6950,'GGG',0,(3,),None],
                    ['r3','chr1',152950,None,0,(),None],
                    ['r4','chr1',156950,None,0,(),None],
                    ['r5','chr1',159901,'AAA',0,(7,),None], # variant on the last position of read
                    ['r6','chr1',165000,'AAA',0,(8,),None], # variant on the first position of read
                    ['r7','chr1',170000,'AAA',0,(9,),None]]

        expected_lines = ['1 chr1:2950-7050:GGG 0 -1 -1 2 11 55',
                          '1 r3 0 -1 -1 5 00 55',
                          '1 chr1:159901-170100:AAA 0 -1 -1 7 111 555']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

class test_paired_end(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        out_dir = 'tests/test_paired_end/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',2300,1),
                       ('chr1',53000,1),
                       ('chr1',53500,1),
                       ('chr1',107000,1)]

        read_info = [['r1','chr1',1950,'GGG',67,(1,),(1,445)],
                    ['r1','chr1',2295,'GGG',131,(),(0,-445)],
                    ['r2','chr1',52950,'GGG',67,(3,),(3,650)],
                    ['r2','chr1',53500,'GGG',131,(4,),(2,-650)],
                    ['r3','chr1',106950,'AAA',67,(5,),(5,100)],
                    ['r3','chr1',106950,'AAA',131,(5,),(4,-100)]]

        expected_lines = ['1 chr1:1950-2395:GGG 0 -1 -1 1 10 55',
                          '1 chr1:52950-53600:GGG 0 -1 -1 3 11 55',
                          '1 chr1:106950-107050:AAA 0 -1 -1 5 1 5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    def test2(self):
        out_dir = 'tests/test_paired_end/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr1',2000,1),
                       ('chr1',2300,1),
                       ('chr1',13000,1),
                       ('chr1',13500,1),
                       ('chr1',27000,1)]

        read_info = [['r1','chr1',1950,'GGG',67,(1,),(1,445)],
                    ['r1','chr1',2295,'GGG',131,(),(0,-445)],
                    ['r2','chr1',12950,'GGG',67,(3,),(3,650)],
                    ['r2','chr1',13500,'GGG',131,(4,),(2,-650)],
                    ['r3','chr1',26950,'AAA',67,(5,),(5,100)],
                    ['r3','chr1',26950,'AAA',131,(5,),(4,-100)]]

        expected_lines = ['1 chr1:1950-13600:GGG 0 -1 -1 1 1011 5555',
                          '1 chr1:26950-27050:AAA 0 -1 -1 5 1 5']

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)


class test_bad_reads(unittest.TestCase):

    # not primary alignment
    def test1(self):
        out_dir = 'tests/test_bad_reads/test1'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr3',160000,1),
                       ('chr3',170000,1),
                       ('chr3',189000,1)]

        read_info = [['r1','chr3',160000,'AAA',0,(1,),None],
                    ['r2','chr3',170000,'AAA',256,(2,),None],
                    ['r3','chr3',189000,'AAA',0,(3,),None]]

        expected_lines = ['1 chr3:160000-160100:AAA 0 -1 -1 1 1 5',
                          '1 chr3:189000-189100:AAA 0 -1 -1 3 1 5',]

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # read is QC fail
    def test2(self):
        out_dir = 'tests/test_bad_reads/test2'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr3',160000,1),
                       ('chr3',170000,1),
                       ('chr3',189000,1)]

        read_info = [['r1','chr3',160000,'AAA',0,(1,),None],
                    ['r2','chr3',170000,'AAA',512,(2,),None],
                    ['r3','chr3',189000,'AAA',0,(3,),None]]

        expected_lines = ['1 chr3:160000-160100:AAA 0 -1 -1 1 1 5',
                          '1 chr3:189000-189100:AAA 0 -1 -1 3 1 5',]

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

    # read is PCR or optical duplicate
    def test3(self):
        out_dir = 'tests/test_bad_reads/test3'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

        variant_loc = [('chr3',160000,1),
                       ('chr3',170000,1),
                       ('chr3',189000,1)]

        read_info = [['r1','chr3',160000,'AAA',0,(1,),None],
                    ['r2','chr3',170000,'AAA',1024,(2,),None],
                    ['r3','chr3',189000,'AAA',0,(3,),None]]

        expected_lines = ['1 chr3:160000-160100:AAA 0 -1 -1 1 1 5',
                          '1 chr3:189000-189100:AAA 0 -1 -1 3 1 5',]

        output_lines = link_fragments_test_helper(out_dir, variant_loc, read_info)
        self.assertEqual(output_lines,expected_lines)

if __name__ == '__main__':
    if not os.path.exists(EXTRACTHAIRS):
        print("ERROR: Please build extractHAIRS executable and place it in {} or specify correct executable at top of this script.".format(EXTRACTHAIRS))
        exit(1)
    unittest.main()
