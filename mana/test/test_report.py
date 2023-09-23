#!/usr/bin/env python
import mana.report as report


def test_1():
    x = report.template()
    assert len(x) > 0
    assert "@@" in x

test_1()
