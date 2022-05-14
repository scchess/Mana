from tools import formatDP


def template_plasmid():
    with open("report_plasmid.txt") as r:
        return r.read()


def generate_plasmid(bam_path, fasta_path, report_path, log_path, pysam1, flag1, flag2, flag3, stats1):
    template_plasmid()
    pass


def template_mRNA():
    with open("report_mrna.txt") as r:
        return r.read()


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

    on_target_c = flag2["total"]
    of_target_c = flag3["total"]
    on_target_p = formatDP(100 * (on_target_c / (on_target_c + of_target_c)))
    of_target_p = formatDP(100 * (of_target_c / (on_target_c + of_target_c)))
    on_target_p = str(on_target_p) + "%"
    of_target_p = str(of_target_p) + "%"
    total_reads = on_target_c + of_target_c

    txt = txt.replace("@@TotalP@@", "100%")
    txt = txt.replace("@@TotalC@@", str(total_reads))
    txt = txt.replace("@@OnTargetP@@", str(on_target_p))
    txt = txt.replace("@@OnTargetC@@", str(on_target_c))
    txt = txt.replace("@@OffTargetP@@", str(of_target_p))
    txt = txt.replace("@@OffTargetC@@", str(of_target_c))

    txt = txt.replace("@@Alignments@@", bam_path)
    txt = txt.replace("@@Reference@@", fasta_path)
    txt = txt.replace("@@AnalysisLog@@", log_path)
    txt = txt.replace("@@Results@@", report_path)

    txt = txt.replace("@@TotalReads@@", str(flag1["total"]))
    txt = txt.replace("@@MappedReads@@", str(flag1["mapped"]))
    txt = txt.replace("@@UnmappedReads@@", str(flag1["secondary"]))
    txt = txt.replace("@@ReadLengthAverage@@", str(stats1["read_average"]))
    txt = txt.replace("@@AverageCoverage@@", str(pysam1["mean_coverage"]))
    txt = txt.replace("@@MinCoverage@@", str(pysam1["min_coverage"]))
    txt = txt.replace("@@MaxCoverage@@", str(pysam1["max_coverage"]))

    return txt
