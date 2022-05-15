import os


def parse(file):
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

    return {"total": total,
            "primary": primary,
            "secondary": secondary,
            "supplementary": supplementary,
            "duplicates": duplicates,
            "primary_duplicates": primary_duplicates,
            "mapped": mapped,
            "primary_mapped": primary_mapped}
