import analysis


def test_1():
    analysis.run("data/GCF_000005845.2_ASM584v2_genomic.fna", "data/simulated.bam", "mRNA")

test_1()
