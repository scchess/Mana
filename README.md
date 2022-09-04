## Summary

Command-line tool that analyses ONT alignments (.BAM) to report quality control statistics.

## Installation

This program requires the following dependencies:

* [SAMTool](https://samtool.org/)
* [Pysam](https://pysam.readthedocs.io/en/latest/api.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)

The dependencies can be installed manually, however, an alternative is `docker`. Building this program on
docker is easy and straightforward. Docker will automatically work out the dependencies.

    git clone https://github.com/scchess/Mana.git
    cd Mana
    docker build -t mana .

## Quick Start

Let's download a sample alignment file for plasmid analysis. This sample file came from a recent experiment from the Mercer lab.

    wget https://www.dropbox.com/s/25bvchax1wgvf5m/cDNA_UnMod_37C_NEBT7_BaseGfpmRNA_1strun_allpassedreads_sorted.bam?dl=1
    mv cDNA_UnMod_37C_NEBT7_BaseGfpmRNA_1strun_allpassedreads_sorted.bam?dl=1 plasmid.bam

Use `docker` to show command-line usage:

    docker run mana python3 mana.py

Run an plasmid analysis:

    docker run -v ${PWD}:/src -i -t mana python3 mana.py --plasmid -b plasmid.bam
    cat outputs/mrna_results.txt
