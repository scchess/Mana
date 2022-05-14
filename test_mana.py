import mana


def test_1():
    x = mana.usage()
    assert "USAGE EXAMPLE" in x
    assert "OUTPUTS" in x
    assert "Oxford" in x
