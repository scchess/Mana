import pandas as pd
from tools import formatDP
from tabulate import tabulate
tabulate.PRESERVE_WHITESPACE = True


def template_plasmid():
    with open("report_plasmid.txt") as r:
        return r.read()


def template_mRNA():
    with open("report_mrna.txt") as r:
        return r.read()


def generate_plasmid(bam_path, fasta_path, report_path, log_path, pysam1, flag1, flag2, flag3, stats1):
    pass


def generate_mRNA(bam_path, fasta_path, report_path, log_path, pysam1, flag1, flag2, flag3, stats1):
    txt = template_mRNA()
    txt = txt.replace("@@MatchA@@", str(formatDP(pysam1["match_mean"])))
    txt = txt.replace("@@MatchMin@@", str(pysam1["match_min"]))
    txt = txt.replace("@@MatchMax@@", str(pysam1["match_max"]))
    txt = txt.replace("@@MismatchA@@", str(formatDP(pysam1["mismatches_mean"])))
    txt = txt.replace("@@MismatchMin@@", str(pysam1["mismatches_min"]))
    txt = txt.replace("@@MismatchMax@@", str(pysam1["mismatches_max"]))
    txt = txt.replace("@@DeleteA@@", str(formatDP(pysam1["deletions_mean"])))
    txt = txt.replace("@@DeleteMin@@", str(pysam1["deletions_min"]))
    txt = txt.replace("@@DeleteMax@@", str(pysam1["deletions_max"]))
    txt = txt.replace("@@InsertA@@", str(formatDP(pysam1["insertions_mean"])))
    txt = txt.replace("@@InsertMin@@", str(pysam1["insertions_min"]))
    txt = txt.replace("@@InsertMax@@", str(pysam1["insertions_max"]))
    txt = txt.replace("@@TotalA@@", str(formatDP(pysam1["total_mean"])))
    txt = txt.replace("@@TotalMin@@", str(pysam1["total_min"]))
    txt = txt.replace("@@TotalMax@@", str(pysam1["total_max"]))

    margin = ""
    x1 = ["Frequency %", "Match", "Mismatch", "Deletion", "Insertion", "Total Error"]
    x2 = [margin, margin, margin, margin, margin, margin]
    x3 = ["avg", str(formatDP(pysam1["match_mean"])), str(formatDP(pysam1["mismatches_mean"])), str(formatDP(pysam1["deletions_mean"])), str(formatDP(pysam1["insertions_mean"])), str(formatDP(pysam1["total_mean"]))]
    x4 = [margin, margin, margin, margin, margin, margin]
    x5 = ["max", str(formatDP(pysam1["match_max"])), str(formatDP(pysam1["mismatches_max"])), str(formatDP(pysam1["deletions_max"])), str(formatDP(pysam1["insertions_max"])), str(formatDP(pysam1["total_max"]))]
    x6 = [margin, margin, margin, margin, margin, margin]
    x7 = ["min", str(formatDP(pysam1["match_min"])), str(formatDP(pysam1["mismatches_min"])), str(formatDP(pysam1["deletions_min"])), str(formatDP(pysam1["insertions_min"])), str(formatDP(pysam1["total_min"]))]

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

    on_target_c = flag2["total"]
    of_target_c = flag3["total"]
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

    txt = txt.replace("@@Alignments@@", bam_path)
    txt = txt.replace("@@Reference@@", fasta_path)
    txt = txt.replace("@@LogPath@@", log_path)
    txt = txt.replace("@@ReportPath@@", report_path)

    txt = txt.replace("@@TotalReads@@", str(flag1["total"]))
    txt = txt.replace("@@MappedReads@@", str(flag1["mapped"]))
    txt = txt.replace("@@UnmappedReads@@", str(flag1["unmapped"]))
    txt = txt.replace("@@ReadLengthAverage@@", str(stats1["read_average"]))
    txt = txt.replace("@@AverageCoverage@@", str(formatDP(pysam1["mean_coverage"])))
    txt = txt.replace("@@MinCoverage@@", str(pysam1["min_coverage"]))
    txt = txt.replace("@@MaxCoverage@@", str(pysam1["max_coverage"]))

    return txt
