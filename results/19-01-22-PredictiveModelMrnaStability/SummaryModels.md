---
title: "Summary of Decay Models"
author: "Santiago Medina"
date: "1/30/2019"
output: 
  html_document:
    keep_md: true
---



Here, I will look at the results of the predictive models that were trained. 
**How do the models compare for these data and which one should be selected for the final model?**




mdl_name     performance
----------  ------------
rpart          0.0952195
treebag        0.1226168
mt             0.1983205
linearreg      0.1989921
enet           0.1996556
mars           0.2002860
pls            0.2015970
cb             0.2174089
nnet           0.2219086
rf             0.2226600
gbm            0.2247788
svmR           0.2452180

![](./figures/results-1.png)<!-- -->


### What about the test set?

I will plot the predictions in the test set

![](./figures/unnamed-chunk-1-1.png)<!-- -->


Here, the top performing models.

![](./figures/top_models-1.png)<!-- -->


In summary the best performing model in the Suport Vector Machine. This model explain ~1/4 of the data variablity.

![](./figures/test_preds-1.png)<!-- -->

