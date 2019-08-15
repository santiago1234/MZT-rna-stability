# Train GradientBoostingRegressor to include prediction intervals

### 1. Hyperparameter Search

** TrainGBMForPredictionIntervals.ipynb |** Performs the hyperparameter search, finds the parameters tha optimize the R2 score.

### 2. Train models

I train 3 models to get the median prediction and a 80% confidence interval. Training code is in the module: `optimalcodon.optimalcodon.projects.rnastability.predictuncertainty`


The 3 models were saved to:
