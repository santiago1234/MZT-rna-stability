---
title: "Lasso Analysis To Find Regulatory K-mers in MZT"
author: "Santiago Medina"
date: "3/22/2019"
output: 
  html_document:
    keep_md: true
---



## Overview

I run a lasso analysis to find regulatory k-mer with a response variable and a set of predictors.

Response: mRNA stability indicator (decay rate), this indicator was estimated with a linear model, TPM  ~ time.
Predictors: 2 PLSs optimality components and all the possible 6-mers.

This analysis is run 4 times:

+ all maternal genes (~4k) for the Ribo0 and PolyA data sets.
+ all maternal genes and no MiR-430 6-mer sits (~3k) for he Ribo0 and PolyA data sets.

In each of the cases above, we performed 6-Kfold validation to obtain reliable estimates.

The output of this analysis if a regularization path (with the LASSO) and the cross-validation for each fold to detect the best $\lambda$ parameter that maximizes the R2 score.


## Lasso Paths

![](./figures/lassopaths-1.png)<!-- -->

## Best $\lambda$

![](./figures/best_lambda-1.png)<!-- -->

For each $\lambda$ in the plot above we have six points coming from the 6s fold in CV.
Now let's get the best lambda value for each case.



sample_condition   which_genes     param_alpha   mean_test_score
-----------------  -------------  ------------  ----------------
wt_polya           all_genes         0.0085317         0.0529950
wt_ribo            all_genes         0.0064617         0.0840074
wt_polya           no_mir_genes      0.0085317         0.0192610
wt_ribo            no_mir_genes      0.0085317         0.0466361



![](./figures/candidates-1.png)<!-- -->

**NOTE**: If a k-mer is not significant for the condition in the heatmap above there is not color.



