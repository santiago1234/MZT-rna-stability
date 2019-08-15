import pandas as pd
import numpy as np
from joblib import load
from optimalcodon.projects.rnastability.predictuncertainty import predict_seq_with_uncertainty

# load model file
model = load(
    "../../data/19-08-08-PredictiveModelStability/predictivemodel.joblib")

# set seed to replicate results
np.random.seed(123)

# choose a random set of 3000 orfs from the test data
orfs = pd.read_csv(
    "../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/validation_set.csv")
unique_orfs = orfs.loc[:, ['coding']].drop_duplicates().sample(3000)


# predict each unique orf for each specie
fish = dict(specie="fish", datatype="aamanitin polya", cell_type="embryo mzt")
mouse = dict(specie="mouse", datatype="slam-seq", cell_type="mES cells")
xenopus = dict(specie="xenopus", datatype="aamanitin ribo",
               cell_type="embryo mzt")
human = dict(specie="human", datatype="slam-seq", cell_type="k562")

# predict the sequences
preds_f = unique_orfs.coding.apply(
    predict_seq_with_uncertainty, models=model, **fish)
preds_m = unique_orfs.coding.apply(
    predict_seq_with_uncertainty, models=model, **mouse)
preds_x = unique_orfs.coding.apply(
    predict_seq_with_uncertainty, models=model, **xenopus)
preds_h = unique_orfs.coding.apply(
    predict_seq_with_uncertainty, models=model, **human)


# add the median prediction to the table
unique_orfs['fish'] = preds_f.map(lambda x: x[1])
unique_orfs['mouse'] = preds_m.map(lambda x: x[1])
unique_orfs['xenopus'] = preds_x.map(lambda x: x[1])
unique_orfs['human'] = preds_h.map(lambda x: x[1])

# save results to file
unique_orfs.to_csv("results_data/orfs_prdictions.csv", index=False)
