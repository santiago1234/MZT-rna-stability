import pickle
import numpy as np
import pandas as pd

# sklearn imports
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor

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
# DECISION TREE REGRESSOR
#########
tree_reg = DecisionTreeRegressor()

tree_grid = {
    'min_samples_split': np.linspace(0.001, .03, 10),
    'max_features': [None],
    'splitter': ['best'],
    'max_depth': np.arange(10, 25)
}

tree_search = modelevaluation.gridsearch(tree_reg, tree_grid, train_x_transformed, train_y, groups, cores=15)


#########
# RANDOM FOREST
#########
rf_reg = RandomForestRegressor(max_depth=14, n_jobs=8)

rf_grid = {
    'n_estimators': np.arange(100, 2000, 150)
}

rf_search = modelevaluation.gridsearch(rf_reg, rf_grid, train_x_transformed, train_y, groups, cores=3)

#########
# ADA BOOST
#########

ada_reg = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4))

ada_grid = {
    'n_estimators': np.arange(10, 300, 100),
    'learning_rate' : [0.01, 0.1,  1],
    'loss' : ['linear', 'square', 'exponential']
}

ada_search = modelevaluation.gridsearch(ada_reg, ada_grid, train_x_transformed, train_y, groups, cores=25)


#########
# GRADIENT BOOSTING
#########
gbm_reg = GradientBoostingRegressor()
gbm_grid = {'loss':['huber'],
          'learning_rate':[0.1, 0.01, 0.001],
          'n_estimators':[700, 1200, 2000],
          'min_samples_split':[4, 8],
          'min_samples_leaf':[2],
          'max_depth':[10, 15, 20],
          'max_features':['log2']}

gbm_search = modelevaluation.gridsearch(gbm_reg, gbm_grid, train_x_transformed, train_y, groups, cores=25)

####################################
# EVALUATE MODELS ON TEST DATA
####################################

mymodels = {
    'decision tree': tree_search.best_estimator_,
    'AdaBoost': ada_search.best_estimator_,
    'gbm': gbm_search.best_estimator_,
    'random forest': rf_search.best_estimator_.set_params(n_jobs=2) # set params for cross validation
}



# save the trained models

for model_name, model_trained in mymodels.items():
    modelevaluation.save_model(model_trained, model_name, "results_data/trained_models/")

# evaluate test data
modelevaluation.eval_models(mymodels, preprocessing, test_x, test_y).to_csv(
    "results_data/val_trees-Ensebmls.csv", index=False)

