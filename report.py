import settings
import pandas as pd
from tools import formatDP
from tabulate import tabulate
tabulate.PRESERVE_WHITESPACE = True


def run(samtools, pysamstats, bedtools, mode):
    assert(mode == "mRNA" or mode == "plasmid")
    with open("report.txt") as r:
        txt = r.read()

    pysam_stats = pysamstats["stats"]
    samtools_stats = samtools["stats"]
    samtools_flag = samtools["flag_stat"]

    txt = txt.replace("@@Alignments@@", samtools["file"])
    txt = txt.replace("@@Reference@@", pysamstats["fasta"])
    txt = txt.replace("@@LogPath@@", settings.LOG_FILE)
    txt = txt.replace("@@ReportPath@@", settings.REPORT_FILE)

    txt = txt.replace("@@TotalReads@@", str(samtools_flag["total"]))
    txt = txt.replace("@@MappedReads@@", str(samtools_flag["mapped"]))
    txt = txt.replace("@@UnmappedReads@@", str(samtools_flag["unmapped"]))
    txt = txt.replace("@@ReadLengthAverage@@", str(samtools_stats["read_average"]))
    txt = txt.replace("@@AverageCoverage@@", str(formatDP(pysam_stats["mean_coverage"])))
    txt = txt.replace("@@MinCoverage@@", str(pysam_stats["min_coverage"]))
    txt = txt.replace("@@MaxCoverage@@", str(pysam_stats["max_coverage"]))

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
    x1 = ["Frequency %", "Match", "Mismatch", "Deletion", "Insertion", "Total Error"]
    x2 = [margin, margin, margin, margin, margin, margin]
    x3 = ["avg", str(formatDP(pysam_stats["match_mean"])), str(formatDP(pysam_stats["mismatches_mean"])), str(formatDP(pysam_stats["deletions_mean"])), str(formatDP(pysam_stats["insertions_mean"])), str(formatDP(pysam_stats["total_mean"]))]
    x4 = [margin, margin, margin, margin, margin, margin]
    x5 = ["max", str(formatDP(pysam_stats["match_max"])), str(formatDP(pysam_stats["mismatches_max"])), str(formatDP(pysam_stats["deletions_max"])), str(formatDP(pysam_stats["insertions_max"])), str(formatDP(pysam_stats["total_max"]))]
    x6 = [margin, margin, margin, margin, margin, margin]
    x7 = ["min", str(formatDP(pysam_stats["match_min"])), str(formatDP(pysam_stats["mismatches_min"])), str(formatDP(pysam_stats["deletions_min"])), str(formatDP(pysam_stats["insertions_min"])), str(formatDP(pysam_stats["total_min"]))]

    df = pd.DataFrame([x1, x2, x3, x4, x5, x6, x7])
    df = df.transpose()
    t1 = tabulate(df, showindex=False, tablefmt="plain")
    t1 = t1.replace("Frequency ", "Frequency              ")
    t1 = t1.replace("Match ", "Match              ")
    t1 = t1.replace("Mismatch ", "Mismatch              ")
    t1 = t1.replace("Deletion ", "Deletion              ")
    t1 = t1.replace("Insertion ", "Insertion              ")
    t1 = t1.replace("Total Error ", "Total Error              ")
    txt = txt.replace("@@Table1@@", t1)

    on_target_c = 1 # flag2["total"]
    of_target_c = 1 # flag3["total"]
    on_target_p = formatDP(100 * (on_target_c / (on_target_c + of_target_c)))
    of_target_p = formatDP(100 * (of_target_c / (on_target_c + of_target_c)))
    on_target_p = str(on_target_p) + "%"
    of_target_p = str(of_target_p) + "%"
    total_reads = on_target_c + of_target_c

    margin = " "
    x1 = ["Reads", "Total", "On-target", "Off-target"]
    x2 = [margin, margin, margin, margin, margin, margin]
    x3 = ["%", "100%", str(on_target_p), str(of_target_p)]
    x4 = [margin, margin, margin, margin, margin, margin]
    x5 = ["count", str(total_reads), str(on_target_c), str(of_target_c)]
    df = pd.DataFrame([x1, x2, x3, x4, x5])
    df = df.transpose()
    t2 = tabulate(df, showindex=False, tablefmt="plain")
    t2 = t2.replace("Reads ", "Reads               ")
    t2 = t2.replace("Total ", "Total               ")
    t2 = t2.replace("On-target ", "On-target               ")
    t2 = t2.replace("Off-target ", "Off-target               ")
    txt = txt.replace("@@Table2@@", t2)

    return txt
