import os
import tools
import settings


BED_PATH = "data/mrna_target.bed"
FLAGSTAT_PATH = settings.TMP_PATH + os.sep + "{}_integrity_flagstat.txt"
INTERSECT_PATH = settings.TMP_PATH + os.sep + "{}_integrity_intersected.bam"


def run(file, bed_path=None, cached=False):
    assert(os.path.exists(file))

    bed_path = BED_PATH if bed_path is None else bed_path
    intersect_path = INTERSECT_PATH.format(os.path.basename(file))
    flagstat_path = FLAGSTAT_PATH.format(os.path.basename(file))

    if not cached:
        tools.run("bedtools intersect -a " + file + " -b " + BED_PATH + " > " + intersect_path)
        tools.run("samtools index " + intersect_path)
        tools.run("samtools flagstat " + intersect_path + " > " + flagstat_path)

    assert(os.path.exists(bed_path))
    assert(os.path.exists(intersect_path))
    assert(os.path.exists(flagstat_path))

    return {"bed_path": bed_path, "intersect_path": intersect_path, "flagstat_path": flagstat_path}
