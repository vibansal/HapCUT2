
## 10X Genomics Linked-Reads

10X Genomics Linked Reads require an extra step to link short reads together into barcoded molecules.

NOTE: It is required that the BAM reads have the BX (corrected barcode) tag.

(1) use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information. This is a necessary precursor step to running HapCUT2.
```
./build/extractHAIRS --10X 1 --bam reads.sorted.bam --VCF variants.VCF --out unlinked_fragment_file
```
(2) use LinkFragments to link fragments into barcoded molecules:
```
python3 utilities/LinkFragments.py --bam reads.sorted.bam --VCF variants.VCF --fragments unlinked_fragment_file --out linked_fragment_file
```
(3) use HAPCUT2 to assemble fragment file into haplotype blocks.
```
./build/HAPCUT2 --nf 1 --fragments linked_fragment_file --VCF variantcalls.vcf --output haplotype_output_file
```

It is assumed that reads with the same barcode occurring within 20 kb of another belong to the same molecule, and reads separated by more than this distance are assigned to separate molecules. This distance can be controlled using the ```-d``` parameter in LinkFragments.

