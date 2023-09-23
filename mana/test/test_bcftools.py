#!/usr/bin/env python
import mana.bcftools as bcftools
import os
import subprocess


'''
1. Download the consensus from https://www.dropbox.com/personal/files%20for%20ted%20and%20mana/files%20and%20results%20plasmid%20testing%20mana%20vs%20manual/manual%20script%20files%20and%20results%20plasmid/vac%2017%20manual%20results%20and%20files

2. md5sum mrna17_plasmid_allpassedreads_sorted_BCF_consensus.fa
> 358d717cfc38ea954e028aee466916b3

3. Download the mana consensus from https://www.dropbox.com/personal/files%20for%20ted%20and%20mana/files%20and%20results%20plasmid%20testing%20mana%20vs%20manual/mana%20files%20and%20results%20plasmid/vac17%20plasmid%20mana

4. md5sum mrna17_plasmid_allpassedreads_sorted.bam_consensus.fa
> a5fc24521c988eae4dec54d7c2838659

5. md5sum data/2022_17_483289_F8.fasta
> 73a79de441ee7bffdfdc7c151af0a90c  data/2022_17_483289_F8.fasta

md5sum data/mrna17_plasmid_allpassedreads_sorted.bam
> bbe7d8c191140d195b44dda43100fb3b  data/mrna17_plasmid_allpassedreads_sorted.bam

# Run variant calling...
bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f data/2022_17_483289_F8.fasta \
	data/mrna17_plasmid_allpassedreads_sorted.bam | bcftools call -c -M --ploidy 1 -Oz -o A.vcf.gz

md5sum A.vcf.gz
> b7eed5709df1cabc0cd9524edc9826e0  A.vcf.gz

bcftools index A.vcf.gz
bcftools consensus -a - -f data/2022_17_483289_F8.fasta A.vcf.gz > A.fa

> md5sum A.fa
358d717cfc38ea954e028aee466916b3  A.fa

> md5sum mrna17_plasmid_allpassedreads_sorted_BCF_consensus.fa

python3 mana.py --plasmid -b data/mrna17_plasmid_allpassedreads_sorted.bam -f data/2022_17_483289_F8.fasta -o data/manatestvac17_plasmid -o out
'''


dataDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "data")

def test_1():
    x = bcftools.run("plasmid", os.path.join(dataDir, "2022_17_483289_F8.fasta"), os.path.join(dataDir, "mrna17_plasmid_allpassedreads_sorted.bam"))
    assert(x["consensus_path"] == ".mana/mrna17_plasmid_allpassedreads_sorted.bam_consensus.fa")
    assert(x["variants_path"] == ".mana/mrna17_plasmid_allpassedreads_sorted.bam_variants.vcf.gz")
    x = subprocess.run(["md5sum", ".mana/mrna17_plasmid_allpassedreads_sorted.bam_consensus.fa"], capture_output=True, text=True).stdout
    assert "358d717cfc38ea954e028aee466916b3" in x


def test_2():
    x = bcftools.run("mRNA", os.path.join(dataDir, "2022_17_483289_F8.fasta"), os.path.join(dataDir, "2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam"), cached=False, targets="2022-14:1290-1908")
    assert(x["consensus_path"] == ".mana/2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam_consensus.fa")
    assert(x["variants_path"] == ".mana/2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam_variants.vcf.gz")
    x = subprocess.run(["md5sum", ".mana/2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam_consensus.fa"], capture_output=True, text=True).stdout
    assert "dddfe987e9f91cbbda458e8e40fc4ff4" in x

test_1()
test_2()
