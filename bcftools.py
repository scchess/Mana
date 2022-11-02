import os
import tools
import settings


BCFTOOLS_VARIANT = settings.TMP_PATH() + os.sep + "{}_variants.vcf.gz"
BCFTOOLS_CONSENSUS = settings.TMP_PATH() + os.sep + "{}_consensus.fa"


def run(mode, ref, file, cached=False):
    assert(os.path.exists(ref))
    assert(os.path.exists(file))

    consensus_path = BCFTOOLS_CONSENSUS.format(os.path.basename(file))
    variants_path = BCFTOOLS_VARIANT.format(os.path.basename(file))
    bed_path = settings.BED_PATH()

    with open(bed_path) as r:
        toks = r.readline().split("\t")
        assert(len(toks) >= 3)
        chrom = toks[0]
        start = toks[1]
        end = toks[2]

    targets = chrom + ":" + start + "-" + end
    assert(mode == "mRNA" or mode == "plasmid")

    if not cached:
        if mode == "plasmid":
            cmd = "bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f {} {} | bcftools call -c -M --ploidy 1 -Oz -o {}"
        else:
            cmd = "bcftools mpileup --targets " + targets + " -d 300000000 --no-BAQ --min-BQ 0 -Ou -f {} {} | bcftools call -c -M --ploidy 1 -Oz -o {}"

        cmd = cmd.format(ref, file, variants_path)
        tools.run(cmd)

        cmd = "bcftools index {}"
        cmd = cmd.format(variants_path)
        tools.run(cmd)

        cmd = "bcftools consensus -f {} {} > {}"
        cmd = cmd.format(ref, variants_path, consensus_path)
        tools.run(cmd)

    assert(os.path.exists(consensus_path))
    return {"consensus_path": consensus_path, "variants_path": variants_path}
