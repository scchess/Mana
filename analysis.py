import os
# import tools
import report
import settings
import flagstat
import samtools
import bcftools
import integrity
import pysamstats


def run(ref, bam, mode, cached=False):
    assert(os.path.exists(ref))
    assert(os.path.exists(bam))
    assert(mode == "mRNA" or mode == "plasmid")

    x1 = samtools.run(bam, cached)
    x2 = pysamstats.run(ref, bam, cached)
    x3 = bcftools.run(ref, bam, cached)

    depth_path = x1["depth_path"]
    stats_path = x1["stats_path"]
    pysam_path = x2["pysam_path"]
    coverage_path = x1["coverage_path"]
    flagstat_path = x1["flagstat_path"]
    consensus_path = x3["consensus_path"]

    assert(os.path.exists(pysam_path))
    assert(os.path.exists(depth_path))
    assert(os.path.exists(stats_path))
    assert(os.path.exists(coverage_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(consensus_path))

    if mode == "mRNA":
        x4 = integrity.run(file)


    return {"pysam1": pysam1, "flag1": flag1, "flag2": flag2, "flag3": flag3, "stats1": stats1}


def run__(x):
    assert(x["mode"] == "mRNA" or x["mode"] == "plasmid")
    assert("bam" in x)
    assert("path" in x)
    assert("fasta" in x)
    assert("log_path" in x)
    assert("report_path" in x)

    if x["mode"] == "mRNA":
        data = mRNA(x)
        pysam1 = data["pysam1"]
        flag1 = data["flag1"]
        flag2 = data["flag2"]
        flag3 = data["flag3"]
        stats1 = data["stats1"]
        txt = report.generate_mRNA(x["bam"], x["fasta"], x["report_path"], x["log_path"], pysam1, flag1, flag2, flag3, stats1)
    else:
        pass
    return txt
