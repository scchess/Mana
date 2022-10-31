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

In order to run the program with your files, it is important to make sure the files are mounted correctly to docker.

Let's take an example, assume the following file paths:

* /data/my_alignments/alignnment.bam
* /home/my_files/my_references/reference.fasta
* /my_mana

To run a plasmid analysis with the file paths, one would do:

    docker run -v /data/my_alignments:/data1
               -v /home/my_files/my_references:/data2
               -v /my_mana_outputs:/data3 ${PWD}:/src
               -i -t mana python3 mana.py
               --plasmid
               -b /data1/alignnment.bam
               -f /data2/reference.fasta
               -o /data3/mana_ouputs
    
Note how the directories are mounted to docker with the "-v" option. For example, the `/data/my_alignments/alignnment.bam` alignment file is broken into:

* `-v /data/my_alignments:/data1`
* `-b /data1/alignnment.bam`
