---
title: "Tunning Machine Learning Models (Results)"
author: "Santiago Medina"
date: "5/8/2019"
output: 
  html_document:
    keep_md: true
---


Here, I show the results of the machine learning tuning process.




## Cross Validation Profiles For Tunned Models

![](./figures/cross_val_models-1.png)<!-- -->

## Diagnostics plots for top 3 models

![](./figures/diagnosticsGBM-1.png)<!-- -->

Observation: It seems that is more difficult to predict unstable genes than stable genes. This is good in the case that we 
want to optimize proteins since the accuracy is better for more stable genes.

## Prediction Across Species and Data Sets

![](./figures/bydataset-1.png)<!-- -->

## Reporter Sequences Predictions

![](./figures/reporters_preds-1.png)<!-- -->


## Learning Curves

![](./figures/learningcurve-1.png)<!-- -->

