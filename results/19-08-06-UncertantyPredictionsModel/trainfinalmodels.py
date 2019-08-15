"""
train and save predictive models
the output is a python dictionary that is saved with the joblib library
to the main data folder: data/19-08-08-PredictiveModelStability/predictivemodel.joblib
"""
from joblib import dump, load
from optimalcodon.projects.rnastability.predictuncertainty import train_gbm_quantiles

outputfile = "../../data/19-08-08-PredictiveModelStability/predictivemodel.joblib"
path_trainingdata = "../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/"

# run the training code
models = train_gbm_quantiles(path_trainingdata)
dump(models, outputfile)
