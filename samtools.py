import os
import tools
import settings
import pandas as pd


DEPTH_PATH = settings.TMP_PATH() + os.sep + "{}_depth.txt"
STATS_PATH = settings.TMP_PATH() + os.sep + "{}_stats.txt"
COVERAGE_PATH = settings.TMP_PATH() + os.sep + "{}_coverage.txt"
FLAGSTAT_PATH = settings.TMP_PATH() + os.sep + "{}_flagstat.txt"
READ_LENGTH_PATH = settings.TMP_PATH() + os.sep + "{}_length.txt"
READ_LENGTH_R = settings.TMP_PATH() + os.sep + "{}_length.R"

UNMAPPED_PATH = settings.TMP_PATH() + os.sep + "{}_unmapped.fq"
REALIGNED_PATH = settings.TMP_PATH() + os.sep + "{}_religned.bam"
REALIGNED_DEPTH_PATH = settings.TMP_PATH() + os.sep + "{}_religned.depth.txt"
REALIGNED_STATS_PATH = settings.TMP_PATH() + os.sep + "{}_religned.stats.txt"
REALIGNED_COVERAGE_PATH = settings.TMP_PATH() + os.sep + "{}_religned.coverage.txt"
REALIGNED_FLAGSTAT_PATH = settings.TMP_PATH() + os.sep + "{}_religned.flagstat.txt"


def parse_stats(file):
    assert(os.path.exists(file))
    read_average = None
    with open(file) as r:
        for line in r:
            toks = line.strip().split("\t")
            if "average length" in line:
                read_average = int(toks[2])
    assert(read_average is not None)
    return {"read_average": read_average}


def parse_flagstat(file):
    assert(os.path.exists(file))
    reads = []
    with open(file) as r:
        for line in r:
            toks = line.strip().split(" ")
            reads.append(int(toks[0]))
    assert(len(reads) > 0)

    total = reads[0]
    primary = reads[1]
    secondary = reads[2]
    supplementary = reads[3]
    duplicates = reads[4]
    primary_duplicates = reads[5]
    mapped = reads[6]
    primary_mapped = reads[7]
    unmapped = total - mapped

    return {"file": file,
            "total": total,
            "primary": primary,
            "unmapped": unmapped,
            "secondary": secondary,
            "supplementary": supplementary,
            "duplicates": duplicates,
            "primary_duplicates": primary_duplicates,
            "mapped": mapped,
            "primary_mapped": primary_mapped}


def run(file, cached=False):
    READ_LENGTH_SVG = settings.OUT_PATH() + os.sep + "{}_length.svg"

    depth_path = DEPTH_PATH.format(os.path.basename(file))
    stats_path = STATS_PATH.format(os.path.basename(file))
    coverage_path = COVERAGE_PATH.format(os.path.basename(file))
    flagstat_path = FLAGSTAT_PATH.format(os.path.basename(file))
    read_length_path = READ_LENGTH_PATH.format(os.path.basename(file))
    read_length_R = READ_LENGTH_R.format(os.path.basename(file))
    read_length_svg = READ_LENGTH_SVG.format(os.path.basename(file))

    unmapped_path = UNMAPPED_PATH.format(os.path.basename(file))
    realigned_path = REALIGNED_PATH.format(os.path.basename(file))
    realigned_depth_path = REALIGNED_DEPTH_PATH.format(os.path.basename(file))
    realigned_stats_path = REALIGNED_STATS_PATH.format(os.path.basename(file))
    realigned_coverage_path = REALIGNED_COVERAGE_PATH.format(os.path.basename(file))
    realigned_flagstat_path = REALIGNED_FLAGSTAT_PATH.format(os.path.basename(file))

    ecoil = "data/GCF_000005845.2_ASM584v2_genomic.fna"
    assert(os.path.exists(ecoil))

    read_count = "read_length.R"
    assert(os.path.exists(read_count))

    with open(read_count) as r:
        R = r.read().replace("@@File@@", read_length_path).replace("@@Output@@", read_length_svg)

    with open(read_length_R, "w") as w:
        w.write(R)

    if not cached:
        tools.run("samtools index " + file)
        tools.run("samtools depth " + file + " > " + depth_path)
        tools.run("samtools stats " + file + " > " + stats_path)
        tools.run("samtools coverage " + file + " > " + coverage_path)
        tools.run("samtools flagstat " + file + " > " + flagstat_path)
        tools.run("samtools view -F 2048 " + file + " | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > " + read_length_path)
        tools.run("Rscript " + read_length_R)
        tools.run("samtools view -b -f 4 " + file + " | samtools fastq > " + unmapped_path)
        tools.run("minimap2 -ax map-ont " + ecoil + " " + unmapped_path + " | samtools sort > " + realigned_path)
        tools.run("samtools index " + realigned_path)
        tools.run("samtools depth " + realigned_path + " > " + realigned_depth_path)
        tools.run("samtools stats " + realigned_path + " > " + realigned_stats_path)
        tools.run("samtools coverage " + realigned_path + " > " + realigned_coverage_path)
        tools.run("samtools flagstat " + realigned_path + " > " + realigned_flagstat_path)

    assert(os.path.exists(depth_path))
    assert(os.path.exists(stats_path))
    assert(os.path.exists(coverage_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(unmapped_path))
    assert(os.path.exists(read_length_svg))
    assert(os.path.exists(realigned_depth_path))
    assert(os.path.exists(realigned_stats_path))
    assert(os.path.exists(realigned_coverage_path))
    assert(os.path.exists(realigned_flagstat_path))

    stats = parse_stats(stats_path)
    pd_depth = pd.read_csv(depth_path, sep="\t")
    flag_stat = parse_flagstat(flagstat_path)
    pd_coverage = pd.read_csv(coverage_path, sep="\t")
    realigned_stats = parse_stats(realigned_stats_path)
    realigned_pd_depth = pd.read_csv(realigned_depth_path, sep="\t")
    realigned_flag_stat = parse_flagstat(realigned_flagstat_path)
    realigned_pd_coverage = pd.read_csv(realigned_coverage_path, sep="\t")

    return {"file": file,
            "stats": stats,
            "pd_depth": pd_depth,
            "flag_stat": flag_stat,
            "depth_path": depth_path,
            "stats_path": stats_path,
            "pd_coverage": pd_coverage,
            "reliagned_stat": realigned_stats,
            "realigned_pd_depth": realigned_pd_depth,
            "realigned_flag_stat": realigned_flag_stat,
            "realigned_pd_coverage": realigned_pd_coverage,
            "read_length_svg": read_length_svg,
            "coverage_path": coverage_path,
            "flagstat_path": flagstat_path,
            "unmapped_path": unmapped_path}
