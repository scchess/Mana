import bcftools
import subprocess


def test_1():
    x = bcftools.run("plasmid", "data/2022_17_483289_F8.fasta", "data/mrna17_plasmid_allpassedreads_sorted.bam")
    assert(x["consensus_path"] == ".mana/mrna17_plasmid_allpassedreads_sorted.bam_consensus.fa")
    assert(x["variants_path"] == ".mana/mrna17_plasmid_allpassedreads_sorted.bam_variants.vcf.gz")
    x = subprocess.run(["md5sum", ".mana/mrna17_plasmid_allpassedreads_sorted.bam_consensus.fa"], capture_output=True, text=True).stdout
    assert "358d717cfc38ea954e028aee466916b3" in x
