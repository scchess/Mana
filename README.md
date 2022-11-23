## Summary

Command-line tool that analyses ONT alignments (.BAM) to report quality control statistics.

## Installation

This program requires the following dependencies:

* [SAMTool 1.5.1](http://www.htslib.org)
* [BCFtools 1.16](http://www.htslib.org)
* [Pysam 0.20.0](https://pysam.readthedocs.io/en/latest/api.html)
* [bedtools 2.30.0](https://bedtools.readthedocs.io/en/latest)

The dependencies can be installed manually (please check their listed websites), however, an
alternative is `docker`. Building this program on docker is easy and straightforward. Docker
will automatically work out the dependencies.

    git clone https://github.com/scchess/Mana.git
    cd Mana
    docker build -t mana .

The installation time will take about 30 minutes.

## Quick Start

Use `docker` to show command-line usage:

    docker run mana python3 mana.py

In order to run the program with your files, it is important to make sure the files are mounted correctly to docker.

Let's take an example, assume the following file paths:

* `/data/my_alignments/alignnment.bam`
* `/home/my_files/my_references/reference.fasta`
* `/my_mana`

To run a plasmid analysis with the file paths, one would do:

    wget https://www.dropbox.com/s/umeo8q5j4epl1ip/2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1
    wget https://www.dropbox.com/s/qdvr1zg2t3j4t1p/2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1
    wget https://www.dropbox.com/s/2wh2dokdtc6pras/cDNA_Mod_37C_NEBT7_BaseGfpmRNA_1strun_passed_pychop_oldcdn_sorted.bam?dl=1
    wget https://www.dropbox.com/s/3zb49m7erdng40f/mrna17_plasmid_allpassedreads_sorted.bam?dl=1
    wget https://www.dropbox.com/s/qdvr1zg2t3j4t1p/2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1
    wget https://www.dropbox.com/s/umeo8q5j4epl1ip/2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1
    wget https://www.dropbox.com/s/o23nlmyo57s60lg/2022_14_482162_F4.fasta?dl=1

    wget https://www.dropbox.com/s/0y7dxqxkn019q5g/plasmid_gfp_ref.fasta?dl=1
    wget https://www.dropbox.com/s/o23nlmyo57s60lg/2022_14_482162_F4.fasta?dl=1
    wget https://www.dropbox.com/s/r9fcm76fsn1gi8l/2022_17_483289_F8.fasta?dl=1
    mv 2022_14_482162_F4.fasta?dl=1 2022_14_482162_F4.fasta

    mv 2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1 2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam
    mv 2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1 2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam
    mv cDNA_Mod_37C_NEBT7_BaseGfpmRNA_1strun_passed_pychop_oldcdn_sorted.bam?dl=1 cDNA_Mod_37C_NEBT7_BaseGfpmRNA_1strun_passed_pychop_oldcdn_sorted.bam
    mv mrna17_plasmid_allpassedreads_sorted.bam?dl=1 mrna17_plasmid_allpassedreads_sorted.bam
    mv 2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1 2022_17_483289_F8_covid_vac_base_passed_pychop_oldcdn_sorted.bam
    mv 2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam?dl=1 2022_14_482162_F4_covid_vac_base_passed_pychop_oldcdn_sorted.bam

    mv plasmid_gfp_ref.fasta?dl=1 plasmid_gfp_ref.fasta
    mv 2022_14_482162_F4.fasta?dl=1 2022_14_482162_F4.fasta
    mv 2022_17_483289_F8.fasta?dl=1 2022_17_483289_F8.fasta

    docker run -v ${PWD}:/src -i -t mana \
               python3 mana.py --mrna \
               -b /src/mrna17_plasmid_allpassedreads_sorted.bam \
               -f /src/2022_17_483289_F8.fasta \
               -o /src/manatestvac17_plasmid

    Note how the directories are mounted to docker with the "-v" option.

The analysis time will be about 10 minutes.
