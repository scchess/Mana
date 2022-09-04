import os
import tools
import settings
import pandas as pd


DEPTH_PATH = settings.TMP_PATH() + os.sep + "{}_depth.txt"
STATS_PATH = settings.TMP_PATH() + os.sep + "{}_stats.txt"
COVERAGE_PATH = settings.TMP_PATH() + os.sep + "{}_coverage.txt"
FLAGSTAT_PATH = settings.TMP_PATH() + os.sep + "{}_flagstat.txt"
UNMAPPED_BAM_PATH = settings.TMP_PATH() + os.sep + "{}_unmapped.bam"


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
    depth_path = DEPTH_PATH.format(os.path.basename(file))
    stats_path = STATS_PATH.format(os.path.basename(file))
    coverage_path = COVERAGE_PATH.format(os.path.basename(file))
    flagstat_path = FLAGSTAT_PATH.format(os.path.basename(file))
    unmapped_path = UNMAPPED_BAM_PATH.format(os.path.basename(file))

    if not cached:
        tools.run("samtools index " + file)
        tools.run("samtools depth " + file + " > " + depth_path)
        tools.run("samtools stats " + file + " > " + stats_path)
        tools.run("samtools coverage " + file + " > " + coverage_path)
        tools.run("samtools flagstat " + file + " > " + flagstat_path)
        tools.run("samtools view -b -f 4 " + file + " > " + unmapped_path)
        tools.run("samtools index " + unmapped_path)

    assert(os.path.exists(depth_path))
    assert(os.path.exists(stats_path))
    assert(os.path.exists(coverage_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(unmapped_path))

    stats = parse_stats(stats_path)
    pd_depth = pd.read_csv(depth_path, sep="\t")
    flag_stat = parse_flagstat(flagstat_path)
    pd_coverage = pd.read_csv(coverage_path, sep="\t")

    #ecoil = "data/GCF_000005845.2_ASM584v2_genomic.fna"
    #assert(os.path.exists(ecoil))

    # if not cached:
    #    tools.run("minimap2 -ax map-ont " + ecoil + " " + unmapped_path + " > $outputpath/$samplefile/$samplefile'_allpassedreads.sam'")

    return {"file": file,
            "stats": stats,
            "pd_depth": pd_depth,
            "flag_stat": flag_stat,
            "depth_path": depth_path,
            "stats_path": stats_path,
            "pd_coverage": pd_coverage,
            "coverage_path": coverage_path,
            "flagstat_path": flagstat_path,
            "unmapped_path": unmapped_path}
