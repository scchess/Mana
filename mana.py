#!/usr/bin/env python
import os
import mana.tools as tools
import mana.settings as settings
import mana.analysis as analysis
import argparse


def create_bed(fasta, out, p1, p2):
    assert(os.path.exists(fasta))
    p1 = int(p1)
    p2 = int(p2)
    if p1 >= p2:
        raise Exception("p2 must be greater than p1")
    with open(fasta) as r:
        for line in r:
            chrom = line.strip().split(" ")[0].replace(">", "")
            break
    tmp = out + os.sep + "mrna_target.bed"
    with open(tmp, "w") as w:
        w.write(chrom + "\t" + str(p1) + "\t" + str(p2))    
    return tmp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="mana - analysis of mRNA manufacture quality using Oxford Nanopore Sequencing")
    parser.add_argument("--plasmid", help="Plasmid analysis", action="store_true")
    parser.add_argument("--mrna", help="mRNA analsysis", action="store_true")
    parser.add_argument("-o", default="outputs")
    parser.add_argument("-b", help="Input BAM file containing ONT reads aligned to plasmid sequence", required=True)
    parser.add_argument("-f", help="Input reference FASTA file of the plasmid sequence used to generate mRNA", required=True)
    parser.add_argument("-ecoli", help="Input FASTA file of the bacterium used for plasmid propagation. Default will be E.coli K12 ASM584v2.", required=False)
    parser.add_argument("-p1", help="Start coordinate of mRNA", required=False)
    parser.add_argument("-p2", help="Last coordinate of mRNA before PolyA tract", required=False)
    args = parser.parse_args()

    if not args.mrna and not args.plasmid:
        raise Exception("Either --plasmid and --mrna must be provided.")
    elif args.mrna and args.plasmid:
        raise Exception("Options --plasmid and --mrna have both been provided. Only one of those can be provided.")
    elif not os.path.exists(args.b):
        raise Exception("File not existed: " + args.b)
    elif args.p1 is not None and args.p2 is None:
        raise Exception("Both -p1 and -p2 options must be provided.")
    elif args.p2 is not None and args.p1 is None:
        raise Exception("Both -p1 and -p2 options must be provided.")

    settings._OUT_PATH = args.o
    mode = "plasmid" if args.plasmid else "mRNA"

    if args.ecoli is None:
        args.ecoli = "data/GCF_000005845.2_ASM584v2_genomic.fna"

    if not os.path.exists(args.ecoli):
        raise Exception(args.ecoli + " not found")
    if not os.path.exists(args.f):
        raise Exception(args.f + " not found")

    if args.p1 is not None:
        p1 = args.p1
        p2 = args.p2
        if not p1.isdigit():
            raise Exception(p1 + " is not a number")
        elif not p2.isdigit():
            raise Exception(p2 + " is not a number")
        tmp = create_bed(args.f, settings.TMP_PATH(), p1, p2)
        assert(os.path.exists(tmp))
        settings._BED_PATH = tmp

    os.system("mkdir -p " + args.o)
    tools.info(mode)
    analysis.run(args.f, args.ecoli, args.b, mode)
