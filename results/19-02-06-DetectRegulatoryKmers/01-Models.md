---
title: "Modelin mRNA Dynamics during MZT"
author: "Santiago Medina"
date: "2/6/2019"
output: 
  html_document:
    keep_md: true
---




## Overview

In section **3.3 Identifying Regulatory K-mers** of my lab notebook I disscuss the main idea of modeling this data. The final goal is to identify candidate regulatory sequences that affect mRNA stability.

I start by testing some models in the wildtype ribo0 data.

The goal of this notebook is just to explore the models.



## Data Viz

We want to model the following distribution.

![](./figures/01_viz-1.png)<!-- -->

The plot above shows sevaral aspects of the data.

+ Time is the top predictor (accounts for the dynamics)
+ Variation increases with time
+ All the other predictor interact with time

The point is to discover which other element (6-mer) affects mRNA dynamics after
controling for this predictors.


## Models

Here, I test a several models:

model 1, time only model

$$
\begin{aligned}
log2FC &\sim Normal(\mu_i, \sigma) \\
\mu_i &= \alpha + \beta_1 t + \beta_2 t^2
\end{aligned}
$$

model 2, model variance with square term on time
$$
\begin{aligned}
log2FC &\sim Normal(\mu_i, \sigma_i) \\
\mu_i &= \alpha + \beta_1 t + \beta_2 t^2 \\
\sigma_i &= \alpha + \gamma_1 t + \gamma_2 t^2
\end{aligned}
$$
model 3, add optimality interaction with time

$$
\begin{aligned}
log2FC &\sim Normal(\mu_i, \sigma_i) \\
\mu_i &= \alpha + \beta_1 t + \beta_2 t^2 + O_i \\
O_i &= PLS_1(\beta_3 + \beta_4 t) + PLS_2(\beta_5 + \beta_6 t) \\
\sigma_i &= \alpha + \gamma_1 t + \gamma_2 t^2
\end{aligned}
$$

model 4, add Mir-430 time interaction

$$
\begin{aligned}
log2FC &\sim Normal(\mu_i, \sigma_i) \\
\mu_i &= \alpha + \beta_1 t + \beta_2 t^2 + O_i + MiR_i\\
O_i &= PLS_1(\beta_3 + \beta_4 t) + PLS_2(\beta_5 + \beta_6 t) \\
MiR_i &= M(\beta_7 + \beta_8 t)\\
\sigma_i &= \alpha + \gamma_1 t + \gamma_2 t^2
\end{aligned}
$$

Variables definition

+ $t$ time in hours
+ $O_i$ codon optimality for gene$_i$: ($PLS_1$ and $PLS_2$ are the PLS components)
+ $M$ number of Mir430 6-mers sites in the 3' UTR

I scale the predictors for all the models and also i substract the time the minimum
time point so starts at 0





```
## 
## SAMPLING FOR MODEL '28e2465dc825a44020f16f6cb98f09f1' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.001405 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 14.05 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 4.97126 seconds (Warm-up)
## Chain 1:                8.75927 seconds (Sampling)
## Chain 1:                13.7305 seconds (Total)
## Chain 1:
## 
## SAMPLING FOR MODEL '6cc4d449cda838d5fb84ec20c4973399' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.005241 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 52.41 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 33.2863 seconds (Warm-up)
## Chain 1:                34.342 seconds (Sampling)
## Chain 1:                67.6283 seconds (Total)
## Chain 1:
## 
## SAMPLING FOR MODEL '6cc4d449cda838d5fb84ec20c4973399' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.003635 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 36.35 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 49.9673 seconds (Warm-up)
## Chain 1:                39.975 seconds (Sampling)
## Chain 1:                89.9422 seconds (Total)
## Chain 1:
## 
## SAMPLING FOR MODEL '6cc4d449cda838d5fb84ec20c4973399' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.00289 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 28.9 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 54.7199 seconds (Warm-up)
## Chain 1:                42.246 seconds (Sampling)
## Chain 1:                96.9659 seconds (Total)
## Chain 1:
```


### Model Checking 


![](./figures/01_params-1.png)<!-- -->




![](./figures/01_mdl_check-1.png)<!-- -->

Ok, the goal here was just to stablish the final model:

$$
\begin{aligned}
log2FC &\sim Normal(\mu_i, \sigma_i) \\
\mu_i &= \alpha + \beta_1 t + \beta_2 t^2 + O_i + MiR_i + K_i\\
O_i &= PLS_1(\beta_3 + \beta_4 t) + PLS_2(\beta_5 + \beta_6 t) \\
MiR_i &= M(\beta_7 + \beta_8 t)\\
K_i &= K(\beta_9 + \beta_{10} t)\\
\sigma_i &= \alpha + \gamma_1 t + \gamma_2 t^2
\end{aligned}
$$

+ $K$ number of k-mer sites in the 3' UTR

The null hypothesis is $\beta_9 = 0$, $\beta_{10} = 0$, which kmers are significant?


I will run the above model for every combination of 6-mers.

See the corresponfing section in the lab notebook. There I describe the prior, but i end up using a lasso prior.
### session info


```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2  brms_2.6.0      Rcpp_1.0.0      gridExtra_2.3  
##  [5] forcats_0.3.0   stringr_1.3.1   dplyr_0.7.8     purrr_0.3.0    
##  [9] readr_1.3.1     tidyr_0.8.2     tibble_2.0.1    ggplot2_3.1.0  
## [13] tidyverse_1.2.1
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-137         matrixStats_0.54.0   xts_0.11-2          
##  [4] lubridate_1.7.4      threejs_0.3.1        httr_1.3.1          
##  [7] rprojroot_1.3-2      rstan_2.18.2         tools_3.5.1         
## [10] backports_1.1.2      R6_2.3.0             DT_0.5              
## [13] mgcv_1.8-24          lazyeval_0.2.1       colorspace_1.3-2    
## [16] withr_2.1.2          processx_3.2.0       tidyselect_0.2.5    
## [19] Brobdingnag_1.2-6    compiler_3.5.1       cli_1.0.1           
## [22] rvest_0.3.2          xml2_1.2.0           shinyjs_1.0         
## [25] labeling_0.3         colourpicker_1.0     scales_1.0.0        
## [28] dygraphs_1.1.1.6     mvtnorm_1.0-8        callr_3.0.0         
## [31] ggridges_0.5.1       StanHeaders_2.18.0   digest_0.6.18       
## [34] rmarkdown_1.10       base64enc_0.1-3      pkgconfig_2.0.2     
## [37] htmltools_0.3.6      htmlwidgets_1.3      rlang_0.3.1.9000    
## [40] readxl_1.1.0         rstudioapi_0.8       shiny_1.1.0         
## [43] bindr_0.1.1          zoo_1.8-4            jsonlite_1.5        
## [46] crosstalk_1.0.0      gtools_3.8.1         inline_0.3.15       
## [49] magrittr_1.5         loo_2.0.0            bayesplot_1.6.0     
## [52] Matrix_1.2-14        munsell_0.5.0        abind_1.4-5         
## [55] stringi_1.2.4        yaml_2.2.0           pkgbuild_1.0.1      
## [58] plyr_1.8.4           grid_3.5.1           parallel_3.5.1      
## [61] promises_1.0.1       crayon_1.3.4         miniUI_0.1.1.1      
## [64] lattice_0.20-35      haven_1.1.2          hms_0.4.2           
## [67] ps_1.1.0             knitr_1.20           pillar_1.3.1.9000   
## [70] igraph_1.2.2         markdown_0.8         shinystan_2.5.0     
## [73] reshape2_1.4.3       stats4_3.5.1         rstantools_1.5.1    
## [76] glue_1.3.0           evaluate_0.11        modelr_0.1.2        
## [79] httpuv_1.4.5         cellranger_1.1.0     gtable_0.2.0        
## [82] assertthat_0.2.0     mime_0.5             xtable_1.8-3        
## [85] broom_0.5.0          codonr_0.1           coda_0.19-2         
## [88] later_0.7.5          viridisLite_0.3.0    rsconnect_0.8.11    
## [91] shinythemes_1.1.2    bridgesampling_0.6-0
```

