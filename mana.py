import os
import tools
import analysis
import argparse


def usage():
    with open("mana.txt") as r:
        return r.read() + "\n\n"


class MyArgumentParser(argparse.ArgumentParser):
    def format_usage(self):
        return usage()

    def format_help(self):
        return usage()


if __name__ == '__main__':
    parser = MyArgumentParser()
    parser.add_argument("--plasmid", help="Plasmid analysis", action="store_true")
    parser.add_argument("--mrna", help="mRNA analsysis", action="store_true")
    parser.add_argument("-o", default="outputs")
    parser.add_argument("-b", help="Input BAM file containing ONT reads aligned to plasmid sequence.", required=True)
    parser.add_argument("-f", help="Input FASTA file of plasmid template sequence (used for in vitro transcription reaction).", required=True)
    parser.add_argument("-m1", help="Start coordinate of mRNA (ie. first base)", required=True)
    parser.add_argument("-m2", help="Last coordinate of mRNA (ie. last base)", required=True)
    parser.add_argument("-p1", help="Start coordinate of polyA tract (ie. first base)", required=False)
    parser.add_argument("-p2", help="Last coordinate of polyA tract (ie. last base)", required=False)
    args = parser.parse_args()

    if not args.mrna and not args.plasmid:
        raise Exception("Either -plasmid and -mrna must be provided.")
    elif args.mrna and args.plasmid:
        raise Exception("Options -plasmid and -mrna have both been provided. Only one of those can be provided.")
    elif not os.path.exists(args.b):
        raise Exception("File not existed: " + args.b)

    tmp_path = ".mana"
    out_path = args.o  # Where output files are saved

    mode = "plasmid" if args.plasmid else "mRNA"
    log_path = out_path + os.sep + "mrna_log.txt"
    report_path = out_path + os.sep + ("mrna_results.txt" if mode == "mRNA" else "plasmid_results.txt")

    os.system("mkdir -p " + args.o)
    with open(log_path, "w") as w1:
        params = {"mode": mode,
                  "bam": args.b,
                  "fasta": args.f,
                  "logger": w1,
                  "bed": "mrna_target_bed.bed",
                  "path": tmp_path,
                  "log_path": log_path,
                  "report_path": report_path,
                  "m1": args.m1,
                  "m2": args.m2}
        tools.info(str(params), w1)
        txt = analysis.run(params)
        with open(report_path, "w") as w2:
            w2.write(txt)
        tools.info("Generated: " + report_path, w1)
