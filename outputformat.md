
#### Description of phased haplotype file (block format). 

Haplotype blocks are printed to the output file (haplotype_output_file). For a given block, column 2 represents
the allele on one chromosome copy (0 for reference, 1 for variant), while column 3 represents
the allele on the other copy.

Each block starts with a block header with the following format:

BLOCK: offset: \<SNV offset\> len: \<SNV span of block\> phased: \<\# SNVs phased\> SPAN: \<base pair span of block\> fragments \<\# of fragments in block\>

Following the header, there is one line per SNV with the following tab-delimited fields:

1. VCF file index (1-based index of the line in the input VCF describing variant)
2. allele on haploid chromosome copy A (0 means reference allele, 1 means variant allele, - means an unphased variant)
3. allele on haploid chromosome copy B (0 means reference allele, 1 means variant allele, - means an unphased variant)
4. chromosome
5. position
6. reference allele (allele corresponding to 0 in column 2 or 3)
7. variant allele (allele corresponding to 1 in column 2 or 3)
8. VCF genotype field (unedited, directly from original VCF)
9. discrete pruning status (1 means pruned, 0 means phased)
10. switch quality: phred-scaled estimated probability that there is a switch error starting at this SNV (0 means switch error is likely, 100 means switch is unlikely)
11. mismatch quality: phred-scaled estimated probability that there is a mismatch [single SNV] error at this SNV (0 means SNV is low quality, 100 means SNV is high quality)

Each block ends with a line with '********'. 

Some variants can be left unphased by HapCUT2 or unphased after post-processing. Such variants have the two alleles as '-'  in columns 2 and 3. 

Field 9 describes the status of the SNV under the discrete SNV pruning method introduced by RefHap (an SNV is pruned if there are equal reads consistent and inconsistent with the phase), with the slight difference that reads are assigned to haplotypes using likelihoods in our implementation. Use the option "--discrete_pruning 1" to automatically prune SNPs ('- -' phasing) based on the value of this field.

The values in field 10 and field 11 are quality scores that range from 0 to 100 (less confident to more confident). Field 10 is useful for controlling switch errors, especially on data types such as error-prone SMRT reads. It is empty by default (".") unless switch error scores are computed using "--error_analysis_mode 1" (compute switch error confidence but leave blocks intact and all SNVs phased for manual pruning later).

Important note: flag "--split_blocks 1" (compute switch error confidence and automatically split blocks using the value of --split_threshold) is currently broken, for the time being use "--error_analysis_mode 1" and manually split using field 10.

Field 11 is useful for controlling mismatch (single SNV) haplotype errors, similarly to field 9. The default behavior of HapCUT2 is to prune individual SNVs for which this confidence is less than 6.98 (probability of error 0.2), as these are highly likely to be errors.

#### sample of output format 

BLOCK: offset: 46 len: 215 phased: 167 SPAN: 81447 fragments 240
46      1       0       chr20   152792  A       C       1|0:1071:157,192:257,285:495:0/1:.:PATMAT       0       .       100.00  3
48      1       0       chr20   153259  C       G       1|0:1171:198,192:316,293:495:0/1:.:PATMAT       0       .       100.00  25
49      0       1       chr20   153376  A       G       0|1:1085:171,190:282,289:495:0/1:.:PATMAT       0       .       100.00  23
52      1       0       chr20   153759  G       A       1|0:879:176,201:140,164:448:0/1:.:PATMAT        0       .       100.00  22
56      1       0       chr20   153898  A       T       1|0:1125:205,194:283,249:643:0/1:.:PATMAT       0       .       100.00  40
60      1       0       chr20   156080  T       C       1|0:1024:190,181:263,241:495:0/1:.:PATMAT       0       .       100.00  3
62      1       0       chr20   156473  A       G       1|0:1143:190,207:287,297:495:0/1:.:PATMAT       0       .       100.00  4
63      0       1       chr20   156486  A       G       0|1:1124:193,174:293,276:527:0/1:.:PATMAT       0       .       100.00  4
********
