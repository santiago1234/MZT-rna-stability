from mztrnastabilitty.utils import find_stop_codon

def test_find_stop_codon():
    s1 = "ATGC"
    s2 = "acgTTTtag"
    s3 = "ggTAGc"
    s4 = "gggTGAtaa"
    assert find_stop_codon(s1) is None, "no stop codon"
    assert find_stop_codon(s2) == 6, "stop in position 6"
    assert find_stop_codon(s3) is None, "should look in 1st frame"
    assert find_stop_codon(s4) == 3, "get always 1st position"
