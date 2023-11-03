## Summary

Command-line tool that analyses ONT alignments (.BAM) to report quality control statistics.

## License

Academic Free License v3.0

## Installation

This program requires the following dependencies:

* [SAMTool 1.5.1](http://www.htslib.org)
* [BCFtools 1.16](http://www.htslib.org)
* [Pysam 0.20.0](https://pysam.readthedocs.io/en/latest/api.html)
* [bedtools 2.30.0](https://bedtools.readthedocs.io/en/latest)

The dependencies can be installed manually (please check their listed websites), however, an
alternative is `docker`. Building this program on docker is easy and straightforward. Docker
will automatically work out the dependencies. A pre-built docker image is available:

    git clone https://github.com/scchess/Mana.git
    cd Mana
    docker pull scchess/mana && docker tag scchess/mana mana

If you would like to build it yourself (>= 20 minutes and not recommended, a pre-built docker image is available), please run the following:  

    git clone https://github.com/scchess/Mana.git
    cd Mana
    docker build -t mana .

## Quick Start

Use `docker` to show command-line usage:

    docker run mana python3 mana.py

In order to run the program with your files, it is important to make sure the files are mounted correctly to docker.

To run an analysis with the file paths, one would do:

    docker run -v ${PWD}:/src -i -t mana \
               python3 mana.py --mrna \
               -b /src/alignment.bam \
               -f /src/file.fasta \
               -o /src/out

Note how the directories are mounted to docker with the "-v" option.

The analysis time will be about 10 minutes.
