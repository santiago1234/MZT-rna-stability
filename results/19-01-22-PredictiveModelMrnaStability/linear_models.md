---
title: "Linear Models"
author: "Santiago Medina"
date: "1/22/2019"
output: 
  html_document:
    keep_md: true
---





### Overview

In this notebook I will evaluate the folling linear models:

+ Linear Regression
+ Partial Least Squares
+ Ridge Regression
+ Lasso Regression
+ Elastic Net


First, I define a common recipe to preprocces the data.





Visualize tuning results.

![](./figures/linear_data_viz_models-1.png)<!-- -->![](./figures/linear_data_viz_models-2.png)<!-- -->![](./figures/linear_data_viz_models-3.png)<!-- -->


![](./figures/linear_unnamed-chunk-1-1.png)<!-- -->

![](./figures/linear_results_linear_models-1.png)<!-- -->


![](./figures/linear_residual_predicted-1.png)<!-- -->


**Conclusion**

The performance of these models is about 20% (R-squared). If I have to pick one, I will choose the PLS with the 3 components since these components can be used for the other modeling purposes (defining new elements).
