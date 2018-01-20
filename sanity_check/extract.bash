./hapcut/extractHAIRS \
    --VCF file.vcf --bam file.bam \
    > hairs 2> log

./HapCUT2/build/extractHAIRS \
    --VCF file.vcf --bam file.bam \
    --out hairs2 \
    > log2 2> err2
