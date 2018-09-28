
header = ''
prev_chrom = None
chroms = ['chr{}'.format(x) for x in range(1,23)] + ['chrX']
output_filenames = ['data/NA12878_hg18_VCFs/{}.vcf'.format(c) for c in chroms]
output_files = {chrom : open(f,'w') for (chrom,f) in zip(chroms,output_filenames)}
with open("data/temp/NA12878_trio_genotypes_original.vcf",'r') as infile:
    for line in infile:
        # header lines
        if line[0] == '#':
            if line[:6] == '#CHROM': # need to remove parent sample labels
                el = line.strip().split('\t')
                line = '\t'.join(el[:9] + [el[11]])
            header += line # add line to header so we can print it to separate chrom files
            continue

        el = line.strip().split('\t')
        assert(len(el) == 12) # VCF line with 3 individuals has 12 elements
        el = el[:9] + [el[11]] # remove parents

        el[0] = 'chr' + el[0]
        chrom = el[0] # add chr label
        if chrom != prev_chrom: # we're on to a new chromosome
            print(header, file=output_files[chrom]) # print header

        new_line = '\t'.join(el)
        # heterozygous variants for NA12878 only
        if el[9][:3] in ['0/1','1/0','0|1','1|0']:
            print(new_line, file=output_files[chrom])

        prev_chrom = chrom

# close new VCFs
for f in output_files.values():
    f.close()
