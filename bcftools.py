import os
import tools
import settings


BCFTOOLS_VARIANT = settings.TMP_PATH + os.sep + "{}_variants.vcf"
BCFTOOLS_CONSENSUS = settings.TMP_PATH + os.sep + "{}_consensus.fa"


def run(ref, file, cached=False):
    assert(os.path.exists(ref))
    assert(os.path.exists(file))

    consensus_path = BCFTOOLS_CONSENSUS.format(file)
    variants_path = BCFTOOLS_VARIANT.format(file)

    if not cached:
        cmd = "bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f {} {} | bcftools call -c -M --ploidy 1 -Oz -o {}"
        cmd = cmd.format(ref, file, variants_path)
        tools.run(cmd)

        cmd = "bcftools index {}"
        cmd = cmd.format(variants_path)
        tools.run(cmd)

        cmd = "bcftools consensus -a - -f {} {} > {}"
        cmd = cmd.replace(ref, variants_path, consensus_path)

    assert(os.path.exists(consensus_path))
    return {"consensus_path": consensus_path, "variants_path": variants_path}
