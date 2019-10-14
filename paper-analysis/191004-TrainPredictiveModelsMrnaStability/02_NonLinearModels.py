import pickle
import numpy as np
import pandas as pd

# sklearn imports
from sklearn.neighbors import KNeighborsRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR

# my module imports
from optimalcodon.projects.rnastability.dataprocessing import get_data, general_preprocesing_pipeline
from optimalcodon.projects.rnastability import modelevaluation


# load data and preprocessing
(train_x, train_y), (test_x, test_y) = get_data(
    '../data/191004-TrainAndTestSets/')
print("{} points for training and {} for testing with {} features".format(
    train_x.shape[0], test_x.shape[0], test_x.shape[1]))

# define the groups for the cross validation
# the group are the gene ids, with this we make sure that
# a gene id is not in training and testing at the same time
groups = train_x.index.values

preprocessing = general_preprocesing_pipeline(train_x)

preprocessing.fit(train_x)
train_x_transformed = preprocessing.transform(train_x)

print(train_x_transformed.shape)


####################################
# GRID-SEARCH (MODEL TUNNING)
####################################


#########
# KNN
#########
# THIS DISTANCE GIVES THE BEST RESULTS
knn_reg = KNeighborsRegressor(weights='distance')

knn_grid = dict(
    n_neighbors=np.arange(5, 10)
)

knn_search = modelevaluation.gridsearch(
    knn_reg, knn_grid, train_x_transformed, train_y, groups, cores=15)


#########
# Kernel Ridge
#########
kerRidge_reg = KernelRidge(kernel='rbf', gamma=0.1)
kerRidge_grid = {"alpha": [0.9, 0.5, 0.2, 0.1],
                 "gamma": np.array([1e-03, 1.e-02, 1.e-01, 1.e+0])}

kerRidge_search = modelevaluation.gridsearch(
    kerRidge_reg, kerRidge_grid, train_x_transformed, train_y, groups, cores=10)

#########
# SVM RBF kernel
#########
svr_reg = SVR()
svr_grid = {'C': [10, 100, 1000],
            'gamma': [0.01, 0.001],
            'kernel': ['rbf']},

svr_search = modelevaluation.gridsearch(
    svr_reg, svr_grid,
    train_x_transformed,
    train_y, groups,
    cores=30)


####################################
# EVALUATE MODELS ON TEST DATA
####################################

mymodels = {
    'knn': knn_search.best_estimator_,
    'kernel_ridge': kerRidge_search.best_estimator_,
    'svm': svr_search.best_estimator_
}

# save the trained models

for model_name, model_trained in mymodels.items():
    modelevaluation.save_model(model_trained, model_name, "results_data/trained_models/")


modelevaluation.eval_models(mymodels, preprocessing, test_x, test_y).to_csv(
    "results_data/val_non-linearmodels.csv", index=False)

####################################
# 10-FOLD CV IN TRAINING SET
####################################
results = modelevaluation.crossvalidation(
    mymodels, train_x_transformed, train_y, groups)
results.to_csv('results_data/cv_linearmodels.csv', index=False)
