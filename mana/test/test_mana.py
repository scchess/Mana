#!/usr/bin/env python
import os
import mana


def test_1():
    x = mana.create_bed("data/mrna_gfp_ref.fasta", "/tmp", 100, 200)
    assert os.path.exists(x)
    with open(x) as r:
        line = r.read()
        assert("mRNA_GFP_ref\t100\t200" in line)
