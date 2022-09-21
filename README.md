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

Use `docker` to show command-line usage:

    docker run mana python3 mana.py

Run a plasmid analysis:

    docker run -v ${PWD}:/src -i -t mana python3 mana.py --plasmid -b plasmid.bam -f plasmid.fasta
    cat outputs/report.txt
