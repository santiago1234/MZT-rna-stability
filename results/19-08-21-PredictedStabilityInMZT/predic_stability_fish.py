import pandas as pd
import numpy as np
from joblib import load

from optimalcodon.projects.rnastability.dataprocessing import get_data, general_preprocesing_pipeline
from optimalcodon.projects.rnastability.predictuncertainty import predict_seq_with_uncertainty

(train_x, train_y), (test_x, test_y) = get_data(
    "../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/")

models = load(
    "../../data/19-08-08-PredictiveModelStability/predictivemodel.joblib")

train_x['y_observed'] = train_y
train_x['median_prediction'] = models['median_gbm'].predict(train_x)
test_x['y_observed'] = test_y
test_x['median_prediction'] = models['median_gbm'].predict(test_x)


preds_all = pd.concat([train_x, test_x])
(
    preds_all[(preds_all.specie == "fish") & (
        preds_all.datatype == "aamanitin polya")]
    .reset_index()
    .loc[:, ['gene_id', 'median_prediction', 'coding']]
    .rename(columns={"median_prediction": 'predicted_stability'})
    .to_csv("results_data/stability_predictions_fish.csv", index=False)
)
