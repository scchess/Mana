import os
import report
import flagstat
import samtools
import pysamstats
from tools import system


def mRNA(x, cached=False):
    assert("bam") in x
    assert("fasta") in x

    m1 = x["m1"]
    m2 = x["m2"]
    bed = x["bed"]
    bam = x["bam"]
    path = x["path"]
    fasta = x["fasta"]
    logger = x.get("logger")

    assert(os.path.exists(bam))
    assert(os.path.exists(fasta))

    bed_ = path + os.sep + "intersect.bed"
    sampled = path + os.sep + "sampled.bam"

    with open(fasta) as r:
        for line in r:
            ref = line.split(" ")[0].replace(">", "")
            with open(bed_, "w") as w:
                w.write(ref + "\t" + str(m1) + "\t" + str(m2))
                break

    system("bedtools intersect -a " + bam + " -b " + bed_ + " | samtools view -bS - > " + sampled, logger)
    bam = sampled

    coverage_path = path + os.sep + os.path.basename(bam) + "_coverage.txt"
    flagstat_path = path + os.sep + os.path.basename(bam) + "_flagstat.txt"
    depth_path = path + os.sep + os.path.basename(bam) + "_depth.txt"
    pysam_path = path + os.sep + os.path.basename(bam) + "_pysam.txt"
    stats_path = path + os.sep + os.path.basename(bam) + "_stats.txt"
    on_target_path = path + os.sep + os.path.basename(bam) + "_on_target.bam"
    off_target_path = path + os.sep + os.path.basename(bam) + "_off_target.bam"
    on_target_flagstat_path = path + os.sep + os.path.basename(bam) + "_on_target_flagstat.txt"
    off_target_flagstat_path = path + os.sep + os.path.basename(bam) + "_off_target_flagstat.txt"

    if not cached:
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
    assert(x["mode"] == "mRNA" or x["mode"] == "plasmid")
    assert("bam" in x)
    assert("bed" in x)
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
