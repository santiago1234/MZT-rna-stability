---
title: "Evaluatio the Model in WT (MZT) conditions"
author: "Santiago Medina"
date: "5/5/2019"
output: 
  html_document:
    keep_md: true
---



This analysis is done for mainly maternal genes.

## The data sets
1. Model Predictions
2. log2FC data
3. Optimality (PLS1 or categories from Ariel's embo2016)
4. MiR-430 sites

## 1. Model Predictions Bassed on the Codon Compositio




## 2. Log Fold Change Data


## 3. Optimality AND MIR-430 data


## Fig 2

![](./figures/fig02-1.png)<!-- -->

![](./figures/fig02aAllGenes-1.png)<!-- -->

## Residual Plots To show MiR-430 are over-estimated

To get the residual I will fit a linear model:

$$
log2fc\text{(Xhrs / 3hrs)} = \alpha + \beta\;\hat{y}
$$
Where $\hat{y}$ is the predicted stability bassed on the codon composition. The linear model is to 
bassically put the predictions in the same scale as the log2FC.


![](./figures/residuals_M6A-1.png)<!-- -->


![](./figures/residuals_MiR430-1.png)<!-- -->

