import os
import analysis
from argparse import ArgumentParser


def usage():
    with open("mana.txt") as r:
        return r.read() + "\n "


if __name__ == '__main__':
    parser = ArgumentParser(add_help=False, usage=usage())
    parser.add_argument("-plasmid", help="Plasmid analysis", required=True)
    parser.add_argument("-mrna", help="mRNA analsysis", required=True)
    parser.add_argument("-b", help="Input BAM file containing ONT reads aligned to plasmid sequence.", required=True)
    parser.add_argument("-p", help="Input FASTA file of plasmid template sequence (used for in vitro transcription reaction).", required=True)
    parser.add_argument("-m1", help="Start coordinate of mRNA (ie. first base)", required=True)
    parser.add_argument("-m2", help="Last coordinate of mRNA (ie. last base)", required=True)
    parser.add_argument("-p1", help="Start coordinate of polyA tract (ie. first base)", required=False)
    parser.add_argument("-p2", help="Last coordinate of polyA tract (ie. last base)", required=False)
    args = parser.parse_args()

    params = {"plasmid": args.plasmid,
              "mrna": args.mrna,
              "bam": args.b,
              "fasta": args.p,
              "start": args.m1,
              "end": args.m2}

    if not os.path.exists(args.b):
        raise Exception("File not existed: " + args.b)

    logger = open("/tmp/log.txt", "w")
    logger.write("Parameters: " + str(params) + "\n")
    analysis.run(params)
    logger.close()
