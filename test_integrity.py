import integrity


def test_1():
    x = integrity.run("data/simulated.bam")
    assert(x["bed_path"] == "data/mrna_target.bed")
    assert(x["intersect_path"] == ".mana/simulated.bam_integrity_intersected.bam")
    assert(x["flagstat_path"] == ".mana/simulated.bam_integrity_flagstat.txt")
