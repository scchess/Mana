## Summary

Command-line tool that analyses ONT alignments (.bam) to report the following quality control statistics:

## Quickstart

This program requests a few dependencies. The easiest option is run it with `docker`. 

* Download docker from https://docs.docker.com/get-docker
* Open a terminal and type `docker`. Make sure `docker` is running.

Build Mana with docker:

    sudo docker build -t mana .

Download some sample files from Dropbox:

    wget https://www.dropbox.com/s/2jp6ej8i2z3g2o9/cDNA_Mod_37C_Pe7_BaseGfpmRNA_allpassedreads_sorted.bam?dl=1
    mv cDNA_Mod_37C_PefT7_BaseGfpmRNA_allpassedreads_sorted.bam?dl=1 cDNA_Mod_37C_PefT7_BaseGfpmRNA_allpassedreads_sorted.bam

Show usage:

    sudo docker run mana python3 mana.py


## Development

   sudo docker build -t mana .
   sudo docker run -v ${PWD}:/data -p 3100:3100 -i -t mana jupyter notebook --allow-root --ip=0.0.0.0 --port=3100
