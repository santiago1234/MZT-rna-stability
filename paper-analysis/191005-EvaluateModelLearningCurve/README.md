# Model Evaluation and Model Selection

Here, I evaluated the set of trained models using a learning curve with the grouped cross fold validation.
Gene ids are used as the group variable.

Analysis:

+ **01-learning_curve.ipynb|** Generates the learning curve.
+ **02-PlotLearningCurve.R|** plots the learning curve across the models
+ **03-TrainFinalModelLasso.ipynb|** Retrain the model using a larger hyperparameter exploration.
+ **04-FinalModel.ipynb|** This trains the final model and saves the model in a sklearn Pipeline
that includes the data preprocessing. The final model is **results-data/final_model.joblib**. Also, 
generates predictions in the test data.
+ **05_R2scoreBoostrap.R|** Boostrap to compute a confidence interval for $R^2$ score.
+ **plots.R|** Visualization


To reproduce the analysis run the scripts in order.
