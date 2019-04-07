# prepare data for analysis
# the output: analysis_data.csv: data frame with
# decay estimate for wt data
# regulatory pathways information:
import numpy as np
import pandas as pd

decay_estimate = pd.read_csv("../19-03-17-ObatainStabilityMeasureInWTConditions/results_data/wt_mRNA_stability_estimate.csv")
pathways = pd.read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")
my_m6a_targets = pd.read_csv("../19-02-26-M6A-targets/results_data/m6atargets_affected_in_mutant.csv")

# tidy the data for analysis
decay_estimate = (decay_estimate[['Gene_ID','sample_condition' ,'estimate']]
        .rename(columns={'estimate': 'decay_rate'})
        )
decay_estimate = decay_estimate.pivot(index='Gene_ID', columns='sample_condition', values='decay_rate').reset_index()

# subset pathwas only use maternal genes
pathways = pathways[pathways.is_maternal][['Gene_ID', 'm6A', 'MiR430', 'PLS1', 'PLS2']]
pathways['m6a_affected_in_mutant'] = pathways.Gene_ID.apply(lambda x: x in my_m6a_targets.gene_id.values)
analysis_data = pd.merge(decay_estimate, pathways, how='inner', on='Gene_ID')
analysis_data.to_csv('results_data/analysis_data.csv', index=False)
