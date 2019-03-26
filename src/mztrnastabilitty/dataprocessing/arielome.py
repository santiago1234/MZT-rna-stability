"""Arielome


Define the path to different project data sets, the paths always start from
the top of the project top path
"""

import os
import fnmatch
import pandas as pd
from mztrnastabilitty.utils import find_stop_codon


def arielome_dpaths(relative_path_to_project):
    """
    Get A dic mapping sample ids to the Arielome raw files
    Args:
        relative_path_to_project (str): the relative path to the
        190108-mzt-rna-stability project
    Returns:
        dict: Sample ID -> relative file to raw arielome file
    """
    path_to_arielome = "data/19-03-23-Arielome/Arielome"
    relative_path_to_arielome = os.path.join(relative_path_to_project,
                                             path_to_arielome)
    if not os.path.isdir(relative_path_to_arielome):
        raise ValueError("wrong path to {}".format(path_to_arielome))
    rawfiles = {}
    for root, _, files in os.walk(relative_path_to_arielome):
        for filename in files:
            if fnmatch.fnmatch(filename, "*.txt"):
                rawfiles[os.path.basename(root)] = os.path.join(root, filename)

    return rawfiles


def getrawfile(filepath):
    """Loads raw Arielome data
    Args:
        filepath (str): the path to the raw arielome data
    Returns:
        Iterator (pd.io.parsers.TextFileReade): since
        this files are big an iterator object is returned to process
        data in chunks
    """
    chunksize = 50000
    colnames = ['seq_id', 'sequence']
    return pd.read_csv(
        filepath, chunksize=chunksize, sep="\t", header=None, names=colnames)


def processarielome(data_ari, process_output=lambda x: x):
    """
    Process the arielome data:
        + finds the position of the 1st premature stop codon
        + creates two new columns: cds (coding sequence), utr3 (3' utr sequence)
    Args:
        data_ari (pd.DataFrame): raw arielome file
        process_output (function): you can pass a function to process the output
            for example to deleter some columns or add something
    Returns:
        pd.DataFrame
    """
    assert 'sequence' in data_ari.columns
    # need to split sequence at the stop position
    data_ari['stop_position'] = data_ari.sequence.apply(find_stop_codon)
    # replace the missing values with the string length
    # missing values are sequences with no stop codon
    data_ari['stop_position'] = data_ari.stop_position.fillna(
        data_ari.sequence.str.len())

    # to spilit in cds and 3'utr sequences
    def split_at_i(seq, i):
        return [seq[:int(i)], seq[int(i):]]

    splits = data_ari.apply(
        lambda x: split_at_i(x.sequence, x.stop_position), axis=1)
    data_ari[['cds', 'utr3']] = pd.DataFrame(
        splits.values.tolist(), index=data_ari.index)

    # apply process_output function to data
    data_ari = process_output(data_ari)
    # output should be pd.DataFrame
    assert isinstance(data_ari, pd.DataFrame)
    return data_ari
