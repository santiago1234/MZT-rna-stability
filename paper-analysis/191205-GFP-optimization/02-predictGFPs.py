import joblib
import pandas as pd
from optimalcodon.projects.rnastability.predictdecay import predict_sequence


mdl = joblib.load("../191005-EvaluateModelLearningCurve/results-data/final_model.joblib")
seqs = pd.read_csv("results-data/gfp_sequences.csv")

mdl_args = dict(specie="human", cell_type="293t", datatype='endogenous')

seqs['predicted_stability'] = seqs.seqs.apply(predict_sequence, mdl=mdl, **mdl_args)


seqs.to_csv('results-data/seq-predictions.csv', index=False)
