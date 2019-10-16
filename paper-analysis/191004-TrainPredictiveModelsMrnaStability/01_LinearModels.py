import numpy as np
import pandas as pd

# sklearn import
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression, Lasso, ElasticNet
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import VarianceThreshold
from sklearn.cross_decomposition import PLSRegression

# my module imports
from optimalcodon.projects.rnastability.dataprocessing import get_data, general_preprocesing_pipeline
from optimalcodon.projects.rnastability import modelevaluation



(train_x, train_y), (test_x, test_y) = get_data('../data/191004-TrainAndTestSets/')
print("{} points for training and {} for testing with {} features".format(
    train_x.shape[0], test_x.shape[0], test_x.shape[1]))


# DATA PRE-PROCESSING

# pre-process Pipeline

preprocessing = Pipeline([
    ('general', general_preprocesing_pipeline(train_x)), # see the code for general_preprocesing_pipeline
    ('polyfeaturs', PolynomialFeatures(degree=2)),
    ('zerovar', VarianceThreshold(threshold=0.0)),
    ('scaling', StandardScaler()) # I scale again not all polynomial features may be with scaled
])


preprocessing.fit(train_x)
train_x_transformed = preprocessing.transform(train_x)
# define the groups for the cross validation
groups = train_x.index.values

# LINEAR-MODEL
lm_reg = Pipeline([
    ('lm', LinearRegression())
])

lm_grid = dict()

# GRID SEARCH

lm_search = modelevaluation.gridsearch(lm_reg, lm_grid, train_x_transformed, train_y, groups)


# PLS
pls_reg = Pipeline([
    ('pls', PLSRegression())
])

pls_grid = dict(
    pls__n_components = np.arange(6, 15, 1)
)
pls_search = modelevaluation.gridsearch(pls_reg, pls_grid, train_x_transformed, train_y, groups)

# LASSO
lasso = Lasso()
alphas = np.logspace(-4, -0.5, 10)
lasso_grid = [{'alpha': alphas}]
lasso_search = modelevaluation.gridsearch(lasso, lasso_grid, train_x_transformed, train_y, groups)

# ENET
enet = ElasticNet()
alphas = np.logspace(-4, -0.5, 10)
enet_grid = [{'alpha': alphas, 'l1_ratio' : np.linspace(0, 1, 5)}]
enet_search = modelevaluation.gridsearch(enet, enet_grid, train_x_transformed, train_y, groups)


# EVALUATE IN VALIDATION TEST
mymodels = {
    'linear_reg': lm_search.best_estimator_,
    'PLS': pls_search.best_estimator_,
    'lasso': lasso_search.best_estimator_,
    'enet': enet_search.best_estimator_
}

# save the trained models

for model_name, model_trained in mymodels.items():
    modelevaluation.save_model(model_trained, model_name, "results_data/trained_models/")

modelevaluation.eval_models(mymodels, preprocessing, test_x, test_y).to_csv("results_data/val_linearmodels.csv", index=False)

