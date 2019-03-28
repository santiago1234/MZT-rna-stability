"""
get a random sample of the arielome
to test k-mer elements
"""
import os
import fnmatch
import pandas as pd

# get the metadata files
path_to_meta = '../../data/19-03-23-Arielome/Arielome/metadata.csv'
meta = pd.read_csv(path_to_meta)

path_to_files = "../../data/19-03-23-Arielome/processed_arielome_data/"
files = [
    os.path.join(path_to_files, x) for x in os.listdir(path_to_files)
    if fnmatch.fnmatch(x, '*premature_stop.csv')
]


def load_and_sample_datapoint(filepath, n=10000):
    dta = pd.read_csv(files[1])
    return dta.sample(n=int(n), random_state=42)


dsets = []
for filepath in files:
    ## get the file id to insert it as column
    print("getting {}".format(filepath))
    idname = os.path.basename(filepath).split('.')[0]
    dsample = load_and_sample_datapoint(filepath)
    dsets.append(dsample.assign(sample_name=idname))

dsets = pd.concat(dsets)
# merge
pd.merge(
    meta, dsets, on='sample_name').to_csv(
        "./results_data/sampledataset.csv", index=False)
