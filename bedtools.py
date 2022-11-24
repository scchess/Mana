import os
import tools
import samtools
import settings


def run(file, bed_path=None, cached=False):
    assert(os.path.exists(file))
    bed_path = settings.BED_PATH()

    with open(bed_path) as r:
        toks = r.readline().split("\t")
        assert(len(toks) >= 3)
        chrom = toks[0]
        start = int(toks[1])
        end = int(toks[2])

    start_bed = "/tmp/start.bed"
    end_bed = "/tmp/end.bed"

    with open(start_bed, "w") as w:
        w.write(chrom + "\t" + str(start) + "\t" + str(start + 1))

    with open(end_bed, "w") as w:
        w.write(chrom + "\t" + str(end) + "\t" + str(end + 1))

    end = "/tmp/end.bam"
    start = "/tmp/start.bam"
    not_start = "/tmp/not_start.bam"
    start_end = "/tmp/start_end.bam"
    start_not_end = "/tmp/start_not_end.bam"
    end_not_start = "/tmp/end_not_start.bam"
    not_start_end = "/tmp/not_start_end.bam"
    not_start_not_end = "/tmp/not_start_not_end.bam"

    end_flagstat = "/tmp/end.bam_flagstat.txt"
    start_flagstat = "/tmp/start.bam_flagstat.txt"
    not_start_flagstat = "/tmp/not_start.bam_flagstat.txt"
    start_end_flagstat = "/tmp/start_end.bam_flagstat.txt"
    start_not_end_flagstat = "/tmp/start_not_end.bam_flagstat.txt"
    end_not_start_flagstat = "/tmp/end_not_start.bam_flagstat.txt"
    not_start_end_flagstat = "/tmp/not_start_end.bam_flagstat.txt"
    not_start_not_end_flagstat = "/tmp/not_start_not_end.bam_flagstat.txt"

    if not cached:
        tools.info("------------------------ BEDTOOLS START ------------------------")
        tools.run("bedtools intersect -a " + file + " -b " + start_bed + " > " + start)
        tools.run("samtools index " + start)
        tools.run("samtools flagstat " + start + " > " + start_flagstat)

        tools.run("bedtools intersect -a " + start + " -b " + end_bed + " > " + start_end)
        tools.run("samtools index " + start_end)
        tools.run("samtools flagstat " + start_end + " > " + start_end_flagstat)

        tools.run("bedtools intersect -a " + start + " -b " + end_bed + " -v > " + start_not_end)
        tools.run("samtools index " + start_not_end)
        tools.run("samtools flagstat " + start_not_end + " > " + start_not_end_flagstat)

        tools.run("bedtools intersect -a " + file + " -b " + end_bed + " > " + end)
        tools.run("samtools index " + end)
        tools.run("samtools flagstat " + end + " > " + end_flagstat)

        tools.run("bedtools intersect -a " + end + " -b " + start_bed + " -v > " + end_not_start)
        tools.run("samtools index " + end_not_start)
        tools.run("samtools flagstat " + end_not_start + " > " + end_not_start_flagstat)

        tools.run("bedtools intersect -a " + file + " -b " + start_bed + " -v > " + not_start)
        tools.run("samtools index " + not_start)
        tools.run("samtools flagstat " + not_start + " > " + not_start_flagstat)

        tools.run("bedtools intersect -a " + not_start + " -b " + end_bed + " > " + not_start_end)
        tools.run("samtools index " + not_start_end)
        tools.run("samtools flagstat " + not_start_end + " > " + not_start_end_flagstat)

        tools.run("bedtools intersect -a " + not_start + " -b " + end_bed + " -v > " + not_start_not_end)
        tools.run("samtools index " + not_start_not_end)
        tools.run("samtools flagstat " + not_start_not_end + " > " + not_start_not_end_flagstat)
        tools.info("------------------------ BEDTOOLS END ------------------------")

    # full_length_reads
    start_end_bam_flagstat = samtools.parse_flagstat(start_end_flagstat)

    # degraded_from5prime
    end_not_start_bam_flagstat = samtools.parse_flagstat(end_not_start_flagstat)

    # degraded_from3prime
    start_not_end_bam_flagstat = samtools.parse_flagstat(start_not_end_flagstat)

    # off-target
    no_start_no_end_bam_flagstat = samtools.parse_flagstat(not_start_not_end_flagstat)

    assert(os.path.exists(bed_path))
    return {"bed_path": bed_path,
            "start_end_bam_flagstat": start_end_bam_flagstat,
            "end_not_start_bam_flagstat": end_not_start_bam_flagstat,
            "start_not_end_bam_flagstat": start_not_end_bam_flagstat,
            "no_start_no_end_bam_flagstat": no_start_no_end_bam_flagstat}
