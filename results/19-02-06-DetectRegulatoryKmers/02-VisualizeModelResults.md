---
title: "Visualize Models, Identify Potential k-mers"
author: "Santiago Medina"
date: "2/8/2019"
output: 
  html_document:
    keep_md: true
---





## Overview






## plot parameters of default model

![](./figures/02_default_model-1.png)<!-- -->


![](./figures/02_eda-1.png)<!-- -->



![](./figures/02_mdlsummary1-1.png)<!-- -->

NOTE: In general there seems to be a positive effect of each k-mer. This maybe because of the UTR length affects mRNA stability as suggested by Mishima 2016. i may need to add a second model where I include the UTR length.


correlate wild type effect with ribo0 effect.

![](./figures/02_corre_wt_polya-1.png)<!-- -->


## Candidates

What is a candidata k-mer?

+ effect in ribo 0
+ effect in poly A
+ WAIC better than the null model

By effect I mean that the 95% of the posterior estimate does not contain 0.

I will write a function to test each element.


```
## How many candidates? 69
```

![](./figures/02_showcandidates-1.png)<!-- -->

## Candidate sequences

Plot profiles effects

![](./figures/02_candidates_heatmpa-1.png)<!-- -->

