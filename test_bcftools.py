import bcftools


def test_1():
    bcftools.run("data/GCF_000005845.2_ASM584v2_genomic.fna", "data/simulated.bam")
