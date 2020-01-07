import pandas as pd
import joblib
from optimalcodon.projects.rnastability.predictdecay import predict_sequence

mdl = joblib.load("../191005-EvaluateModelLearningCurve/results-data/final_model.joblib")
mdl_args = dict(specie="human", cell_type="293t", datatype='endogenous')
silent_reporters = pd.read_csv('silent_reporters_seqs.csv').drop(['Unnamed: 0'], axis=1)

# model prediction

silent_reporters['prediction'] = silent_reporters.seqs.apply(predict_sequence, mdl=mdl, **mdl_args)
silent_reporters.to_csv('silent_predictions.csv', index=False)
