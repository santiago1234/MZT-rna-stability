"""Utility functions for data processing and other tasks

TODO:
    add a better description for this module
"""


def find_stop_codon(sequence):
    """
    finds the nucleotide position of the first stop codon
    Note: the search is done in the 1st fram of the sequence
    Args:
        sequence (str): dna sequence 4 letter code
    Returns:
        int: nucleotide position of stop codon
        None: sequence does not contain stop codon
    """
    stop_codons = ['TAG', 'TAA', 'TGA']
    sequence = sequence.upper()
    # iterate over the codons
    for i in range(0, len(sequence) - 2, 3):
        if sequence[i: i + 3] in stop_codons:
            return i
    return None  # no stop codon in sequence
