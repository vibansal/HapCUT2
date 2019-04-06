
chroms = ['chr{}'.format(x) for x in range(1,23)] + ['chrX']

for c in chroms:
    duitama_chrom_pos = []
    with open("data/temp/duitama_haplotypes/{}.real_refhap.phase".format(c),'r') as duitama_file:
        for line in duitama_file:
            el = line.strip().split()
            pos = int(el[0])
            duitama_chrom_pos.append((c, pos))

    NA12878_chrom_pos = []
    with open("data/NA12878_hg18_VCFs/{}.vcf".format(c),'r') as NA12878_vcf:
        for line in NA12878_vcf:
            # header lines
            if line[0] == '#':
                continue

            el = line.strip().split('\t')
            chrom = el[0] # add chr label
            pos = int(el[1])

            NA12878_chrom_pos.append((chrom,pos))

    print("checking that indices in Duitama haplotype with {} elements matches indices in NA12878 haplotype with {} elements...".format(len(duitama_chrom_pos), len(NA12878_chrom_pos)))
    assert(duitama_chrom_pos == NA12878_chrom_pos)
    print("...PASSED")
