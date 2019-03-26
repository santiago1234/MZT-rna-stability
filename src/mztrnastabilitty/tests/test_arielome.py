import pandas as pd
from mztrnastabilitty.dataprocessing import arielome


def test_processarielome():
    """
    test the processing of the arielome
    """
    # define toy data
    arielome_input = pd.DataFrame({
        'seq_id': ['a', 'b'],
        'sequence': ['atgCCCgt', 'atgTAAccc']
    })
    # expected results
    expected_output = pd.DataFrame({
        'seq_id': ['a', 'b'],
        'sequence': ['atgCCCgt', 'atgTAAccc'],
        'cds': ['atgCCCgt', 'atg'],
        'utr3': ['', 'TAAccc']
    })

    results = arielome.processarielome(arielome_input).drop(
        'stop_position', axis=1)
    # correct 3' utr sequence
    assert all(results.utr3 == expected_output.utr3)
    # correct cds
    assert all(results.cds == expected_output.cds)
    # merging both should be equal to the sequence
    assert all((results.cds + results.utr3) == results.sequence)
