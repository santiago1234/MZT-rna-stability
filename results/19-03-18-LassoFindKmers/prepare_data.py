import pandas as pd
from sklearn.model_selection import train_test_split

# load data sets
# 3' utrs composition
utrs = pd.read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_3utr6mer_composition.csv")
utrs = utrs.rename(columns={'ensembl_gene_id': 'Gene_ID'}).drop('3utr', axis=1)
# optimality
pls = pd.read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")
# keep the columns we need
# drop non maternal genes
pls = pls[pls['is_maternal']][['Gene_ID', 'PLS1', 'PLS2', 'is_maternal']]

# join the data
dta = pd.merge(pls, utrs, how='inner', on='Gene_ID').drop('is_maternal',
                                                          axis=1)

# split into a train and test set for validation
train, test = train_test_split(dta, test_size=.2, random_state=42)

# save predictors
test.to_csv('results_data/kpreds_test.csv', index=False)
train.to_csv('results_data/kpreds_train.csv', index=False)
