import os
import report
import samtools
import bcftools
import bedtools
import settings
import pysamstats


def run(ref, ecoil, bam, mode, cached=True):
    assert(os.path.exists(ref))
    assert(os.path.exists(bam))
    assert(os.path.exists(ecoil))
    assert(mode == "mRNA" or mode == "plasmid")

    x1 = samtools.run(bam, cached)
    x2 = pysamstats.run(ref, bam, cached)
    x3 = bedtools.run(bam)
    x4 = bcftools.run(ref, bam)

    bed_path = x3["bed_path"]
    depth_path = x1["depth_path"]
    stats_path = x1["stats_path"]
    pysam_path = x2["pysam_path"]
    coverage_path = x1["coverage_path"]
    flagstat_path = x1["flagstat_path"]
    flagstat_path = x3["flagstat_path"]
    intersect_path = x3["intersect_path"]
    consensus_path = x4["consensus_path"]

    assert(os.path.exists(bed_path))
    assert(os.path.exists(pysam_path))
    assert(os.path.exists(depth_path))
    assert(os.path.exists(stats_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(coverage_path))
    assert(os.path.exists(flagstat_path))
    assert(os.path.exists(intersect_path))
    assert(os.path.exists(consensus_path))

    with open(settings.REPORT_FILE, "w") as w:
        w.write(report.run(x1, x2, x3, mode))

    return {"bed_path": bed_path,
            "depth_path": depth_path,
            "stats_path": stats_path,
            "pysam_path": pysam_path,
            "report_path": settings.REPORT_FILE,
            "flagstat_path": flagstat_path,
            "coverage_path": coverage_path,
            "flagstat_path": flagstat_path,
            "intersect_path": intersect_path,
            "consensus_path": consensus_path}
