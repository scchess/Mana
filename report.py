import os
import sys
import settings
import pandas as pd
from tools import formatDP
from tabulate import tabulate
tabulate.PRESERVE_WHITESPACE = True


def run(samtools, pysamstats, bedtools, bcftools, mode):
    assert(mode == "mRNA" or mode == "plasmid")
    with open("report.txt") as r:
        txt = r.read()

    os.system("cp " + bcftools["consensus_path"] + " " + settings.OUT_PATH())

    pysam_stats = pysamstats["stats"]
    samtools_stats = samtools["stats"]
    samtools_flag = samtools["flag_stat"]
    realigned_flag_stat = samtools["realigned_flag_stat"]
    start_end_bam_flagstat = bedtools["start_end_bam_flagstat"]["total"]
    end_not_start_bam_flagstat = bedtools["end_not_start_bam_flagstat"]["total"]
    start_not_end_bam_flagstat = bedtools["start_not_end_bam_flagstat"]["total"]
    no_start_no_end_bam_flagstat = bedtools["no_start_no_end_bam_flagstat"]["total"]

    txt = txt.replace("@@Mode@@", mode)
    txt = txt.replace("@@Command@@", " ".join(sys.argv[:]))
    txt = txt.replace("@@Alignments@@", samtools["file"])
    txt = txt.replace("@@Reference@@", pysamstats["fasta"])
    txt = txt.replace("@@LogPath@@", settings.LOG_FILE())
    txt = txt.replace("@@ReadLengthHistogram@@", samtools["read_length_svg"])
    txt = txt.replace("@@ReportPath@@", settings.REPORT_FILE())
    txt = txt.replace("@@Consensus@@", settings.OUT_PATH() + os.sep + os.path.basename(bcftools["consensus_path"]))

    txt = txt.replace("@@TotalReads@@", str(samtools_flag["primary"]))
    txt = txt.replace("@@MappedReads@@", str(samtools_flag["primary_mapped"]))
    txt = txt.replace("@@UnmappedReads@@", str(samtools_flag["unmapped"]))
    txt = txt.replace("@@ReadLengthAverage@@", str(samtools_stats["read_average"]))
    txt = txt.replace("@@AverageCoverage@@", str(formatDP(pysam_stats["mean_coverage"])))
    txt = txt.replace("@@MinCoverage@@", str(pysam_stats["min_coverage"]))
    txt = txt.replace("@@MaxCoverage@@", str(pysam_stats["max_coverage"]))

    txt = txt.replace("@@E_TotalReads@@", str(realigned_flag_stat["total"]))
    txt = txt.replace("@@E_MappedReads@@", str(realigned_flag_stat["mapped"]))
    txt = txt.replace("@@E_UnmappedReads@@", str(realigned_flag_stat["unmapped"]))

    txt = txt.replace("@@MatchA@@", str(formatDP(pysam_stats["match_mean"])))
    txt = txt.replace("@@MatchMin@@", str(pysam_stats["match_min"]))
    txt = txt.replace("@@MatchMax@@", str(pysam_stats["match_max"]))
    txt = txt.replace("@@MismatchA@@", str(formatDP(pysam_stats["mismatches_mean"])))
    txt = txt.replace("@@MismatchMin@@", str(pysam_stats["mismatches_min"]))
    txt = txt.replace("@@MismatchMax@@", str(pysam_stats["mismatches_max"]))
    txt = txt.replace("@@DeleteA@@", str(formatDP(pysam_stats["deletions_mean"])))
    txt = txt.replace("@@DeleteMin@@", str(pysam_stats["deletions_min"]))
    txt = txt.replace("@@DeleteMax@@", str(pysam_stats["deletions_max"]))
    txt = txt.replace("@@InsertA@@", str(formatDP(pysam_stats["insertions_mean"])))
    txt = txt.replace("@@InsertMin@@", str(pysam_stats["insertions_min"]))
    txt = txt.replace("@@InsertMax@@", str(pysam_stats["insertions_max"]))
    txt = txt.replace("@@TotalA@@", str(formatDP(pysam_stats["total_mean"])))
    txt = txt.replace("@@TotalMin@@", str(pysam_stats["total_min"]))
    txt = txt.replace("@@TotalMax@@", str(pysam_stats["total_max"]))

    margin = ""
    x1 = ["Frequency %", "Match", "Insertion", "Mismatch", "Deletion", "Total Error"]
    x2 = [margin, margin, margin, margin, margin, margin]
    x3 = ["avg", str(formatDP(pysam_stats["match_mean"])), str(formatDP(pysam_stats["insertions_mean"])), str(formatDP(pysam_stats["mismatches_mean"])), str(formatDP(pysam_stats["deletions_mean"])), str(formatDP(pysam_stats["total_mean"]))]
    x4 = [margin, margin, margin, margin, margin, margin]
    x5 = ["max", str(formatDP(pysam_stats["match_max"])), str(formatDP(pysam_stats["insertions_max"])), str(formatDP(pysam_stats["mismatches_max"])), str(formatDP(pysam_stats["deletions_max"])), str(formatDP(pysam_stats["total_max"]))]
    x6 = [margin, margin, margin, margin, margin, margin]
    x7 = ["min", str(formatDP(pysam_stats["match_min"])), str(formatDP(pysam_stats["insertions_min"])), str(formatDP(pysam_stats["mismatches_min"])), str(formatDP(pysam_stats["deletions_min"])), str(formatDP(pysam_stats["total_min"]))]

    df = pd.DataFrame([x1, x2, x3, x4, x5, x6, x7])
    df = df.transpose()
    t1 = tabulate(df, showindex=False, tablefmt="plain")
    t1 = t1.replace("Frequency ", "Frequency              ")
    t1 = t1.replace("Match ", "Match              ")
    t1 = t1.replace("Insertion ", "Insertion              ")
    t1 = t1.replace("Mismatch ", "Mismatch              ")
    t1 = t1.replace("Deletion ", "Deletion              ")
    t1 = t1.replace("Total Error ", "Total Error              ")
    txt = txt.replace("@@Table1@@", t1)

    x1 = ["Full Length Reads", "Degraded from 5-prime Reads", "Degraded from 3-prime Reads", "Off-target Reads"]
    x3 = [str(start_end_bam_flagstat), str(end_not_start_bam_flagstat), str(start_not_end_bam_flagstat), str(no_start_no_end_bam_flagstat)]
    df = pd.DataFrame([x1, x3])
    df = df.transpose()
    t2 = tabulate(df, showindex=False, tablefmt="plain")
    t2 = t2.replace("Total ", "Total               ")
    t2 = t2.replace("Degraded from 5-prime Reads", "Degraded from 5-prime Reads    ")
    t2 = t2.replace("Degraded from 3-prime Reads", "Degraded from 3-prime Reads  ")
    t2 = t2.replace("Off-target Reads", "Off-target Reads      ")

    if mode == "mRNA":
        txt = txt.replace("@@Table2@@", t2)
    else:
        txt = txt.replace("@@Table2@@", "")

    return txt
