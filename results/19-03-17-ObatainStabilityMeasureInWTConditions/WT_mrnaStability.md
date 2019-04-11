---
title: "Obtain Indicator of mRNA stability in WT conditions"
author: "Santiago Medina"
date: "3/17/2019"
output: 
  html_document:
    keep_md: true
---



## Overview

The goal here is to obtain an indicator of mRNA stability in the wild type conditions. Previously I've been using the log2 fold change between 6 hrs and 2hrs. But here I fit an exponential model (linear in the log2 scale) similar to the decay model, to obtain this indicator, I only use the time points 3, 4, 5, 6 , 7, 8.



for each gene I fit the following model.

$$
logRNA = \alpha + \beta\;time
$$

The $\beta$ parameter is our stability indicator. Note that this parameter will only be meaningful for purely maternal genes since
the zygotic expression will affect the parameter.

![](./figures/stability-1.png)<!-- -->

Test check

+ MiR-430 targets should go down
+ We should see optimality

At least for the maternal genes.

![](./figures/unnamed-chunk-2-1.png)<!-- -->

![](./figures/unnamed-chunk-3-1.png)<!-- -->

### Figure 2 (Codon composition explains maternal stability)

Here, I need to load the predictive model.




define 500 optimal and non optimal codons bassed on pls 1

![](./figures/fig02-1.png)<!-- -->

save the stability estimate in the wild type conditions



