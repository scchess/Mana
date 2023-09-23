#!/usr/bin/env python
import math
import mana.pysamstats as pysamstats


#
# https://www.dropbox.com/personal/files%20for%20ted%20and%20mana/files%20and%20results%20plasmid%20testing%20mana%20vs%20manual/manual%20script%20files%20and%20results%20plasmid/vac%2014%20plasmid%20manual%20results%20and%20files?preview=mrna14_plasmid_passed_pysamreporttable.txt
#

def test_1():
    x = pysamstats.parse("data/mrna14_plasmid_allpassedreads_sorted.bam_pysam.txt")
    print(x)

    assert x["total_min"] == 0
    assert math.isclose(x["total_max"], 0.951970989975462)
    assert math.isclose(x["total_mean"], 0.03042945189367441)

    assert math.isclose(x["match_mean"], 0.9391410962126511)
    assert math.isclose(x["match_min"], 0.027519324662708216)
    assert x["match_max"] == 1

    assert math.isclose(x["insertions_mean"], 0.011543892755475035)
    assert math.isclose(x["insertions_min"], 0.0)
    assert math.isclose(x["insertions_max"], 0.22297016192108635)

    assert math.isclose(x["insertions_mean"], 0.011543892755475035)
    assert math.isclose(x["insertions_min"], 0.0)
    assert math.isclose(x["insertions_max"], 0.22297016192108635)

    assert math.isclose(x["deletions_mean"], 0.032999462220377335)
    assert math.isclose(x["deletions_min"], 0.0)
    assert math.isclose(x["deletions_max"], 0.8527974797319545)

    assert math.isclose(x["mismatches_mean"], 0.027859441566971482)
    assert math.isclose(x["mismatches_min"], 0.0)
    assert math.isclose(x["mismatches_max"], 0.951970989975462)

    assert math.isclose(x["mean_coverage"], 119458.85762711865)
    assert math.isclose(x["min_coverage"], 40873)
    assert math.isclose(x["max_coverage"], 126687)
