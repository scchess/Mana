import os


def parse(file):
    assert(os.path.exists(file))
    read_average = None
    with open(file) as r:
        for line in r:
            toks = line.strip().split("\t")
            if "average length" in line:
                read_average = int(toks[2])
    return {"read_average": read_average}
