"""Prepare Data for Analysis
I use the following filters:
    cds length in nucleotides >= 150
    utr length in nucleotides >= 200
Then the 3' utrs are collapsed into a 
single utr separated by "*" 
"""
import os
import fnmatch
from functools import reduce
import numpy as np
import pandas as pd

path_to_files = "../../data/19-03-23-Arielome/processed_arielome_data/"
files = [
    os.path.join(path_to_files, x) for x in os.listdir(path_to_files)
    if fnmatch.fnmatch(x, '*premature_stop.csv')
]


def combine_utrs(data):
    join_utrs = lambda x, y: x + "*" + y
    return reduce(join_utrs, data.utr3.values)


def collapse_utr(dta, collapse_n=100):
    """
    filter the data:
        cds length in nucleotides >= 150
        utr length in nucleotides >= 200
    and collapses the utrs every n records,
    utrs are separated by *
    Args:
        data (pd.DataFrame): data to collapse
        n (int): number of sequence to collapse
    Returns:
        pd.DataFrame: collapsed 3' utrs
    """
    filters = (dta.utr3.str.len() >= 200) & (dta.cds.str.len() >= 150)
    dta = dta[filters]
    grouping_var = np.arange(dta.shape[0]) % int(dta.shape[0] / collapse_n)
    collapsed_utrs = dta.groupby(grouping_var).apply(combine_utrs)
    return pd.DataFrame(collapsed_utrs, columns=['collapsed_utrs'])


dsets = []
for filepath in files:
    ## get the file id to insert it as column
    print("getting {}".format(filepath))
    idname = os.path.basename(filepath).split('.')[0]
    # load data set
    dset = pd.read_csv(filepath)
    dset_utrs = collapse_utr(dset)
    dset_utrs['sample_name'] = idname
    dsets.append(dset_utrs)

dsets = pd.concat(dsets)
# load metadata
path_to_meta = '../../data/19-03-23-Arielome/Arielome/metadata.csv'
meta = pd.read_csv(path_to_meta)
pd.merge(
    meta, dsets, on='sample_name').to_csv(
        "./results_data/collapsed_3utrs_every100.csv", index=False)
