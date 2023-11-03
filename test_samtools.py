import os
import samtools


def test_1():
    os.system("rm -f .mana/*.txt")
    x = samtools.run("data/simulated.bam")
    assert(x["depth_path"] == ".mana/simulated.bam_depth.txt")
    assert(x["stats_path"] == ".mana/simulated.bam_stats.txt")
    assert(x["coverage_path"] == ".mana/simulated.bam_coverage.txt")
    assert(x["flagstat_path"] == ".mana/simulated.bam_flagstat.txt")
