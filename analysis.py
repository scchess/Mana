import os
import report
import flagstat
import samtools
import pysamstats
from tools import system


def mRNA(x):
    assert("bam") in x
    assert("fasta") in x

    bam = x["bam"]
    bed = x["bed"]
    path = x["path"]
    fasta = x["fasta"]
    logger = x.get("logger")

    assert(os.path.exists(bed))
    assert(os.path.exists(bam))
    assert(os.path.exists(fasta))

    coverage_path = path + os.sep + os.path.basename(bam) + "_coverage.txt"
    flagstat_path = path + os.sep + os.path.basename(bam) + "_flagstat.txt"
    depth_path = path + os.sep + os.path.basename(bam) + "_depth.txt"
    pysam_path = path + os.sep + os.path.basename(bam) + "_pysam.txt"
    stats_path = path + os.sep + os.path.basename(bam) + "_stats.txt"
    on_target_path = path + os.sep + os.path.basename(bam) + "_on_target.bam"
    off_target_path = path + os.sep + os.path.basename(bam) + "_off_target.bam"
    on_target_flagstat_path = path + os.sep + os.path.basename(bam) + "_on_target_flagstat.txt"
    off_target_flagstat_path = path + os.sep + os.path.basename(bam) + "_off_target_flagstat.txt"

    system("mkdir -p " + path, logger)
    system("samtools index " + bam, logger)
    system("samtools stats " + bam + " > " + stats_path, logger)
    system("samtools coverage " + bam + " > " + coverage_path, logger)
    system("samtools flagstat " + bam + " > " + flagstat_path, logger)
    system("samtools depth " + bam + " > " + depth_path, logger)
    system("pysamstats --max-depth=3000000 --fasta " + fasta + " --type variation " + bam + " > " + pysam_path, logger)
    system("bedtools intersect -a " + bam + " -b " + bed + " > " + on_target_path, logger)
    system("bedtools intersect -a " + bam + " -v -b " + bed + " > " + off_target_path, logger)
    system("samtools index " + on_target_path, logger)
    system("samtools index " + off_target_path, logger)
    system("samtools flagstat " + on_target_path + " > " + on_target_flagstat_path, logger)
    system("samtools flagstat " + off_target_path + " > " + off_target_flagstat_path, logger)

    assert(os.path.exists(depth_path))
    assert(os.path.exists(pysam_path))
    assert(os.path.exists(stats_path))
    assert(os.path.exists(coverage_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(on_target_flagstat_path))
    assert(os.path.exists(off_target_flagstat_path))

    pysam1 = pysamstats.parse(pysam_path)
    flag1 = flagstat.parse(flagstat_path)
    flag2 = flagstat.parse(on_target_flagstat_path)
    flag3 = flagstat.parse(off_target_flagstat_path)
    stats1 = samtools.parse(stats_path)

    return {"pysam1": pysam1, "flag1": flag1, "flag2": flag2, "flag3": flag3, "stats1": stats1}


def run(x):
    assert(x["analysis"] == "mRNA" or x["analysis"] == "plasma")
    if x["analysis"] == "mRNA":
        data = mRNA(x)
        pysam1 = data["pysam1"]
        flag1 = data["flag1"]
        flag2 = data["flag2"]
        flag3 = data["flag3"]
        stats1 = data["stats1"]
        txt = report.generate_mRNA(x["bam"], x["fasta"], x["report"], x["log"], pysam1, flag1, flag2, flag3, stats1)
        with open(x["report"], "w") as w:
            w.write(txt)
    return txt
