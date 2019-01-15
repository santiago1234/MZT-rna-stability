---
title: "Decay Rate Estimation"
author: "Santiago Medina"
date: "1/11/2019"
output:
   html_document:
    keep_md: true
---

### Overview:

The goal of this analysis is to estimate the decay rate for the endogenous genes with the a-amanitin data.



### 1. Tidy data for analysis

I use the a-amanitin data set. This includes 12 time points.



### 2. Detect active genes

To estimate the decay rate I will use only active genes. I will use the 
[zFPKM](https://bioconductor.org/packages/release/bioc/html/zFPKM.html) package for this.

Plot of the distribution of gene expression.

![](./figures/expression_distribution-1.png)<!-- -->

as the authors suggest a cut-off of zFPKM > -3 is used.


![](./figures/active_genes-1.png)<!-- -->

### 3. Estimate Decay Rate

Now, I display the time course for a set of random genes.

![](./figures/unnamed-chunk-1-1.png)<!-- -->

#### Model

To estimate the decay rate I will use the following linear model

$$
\begin{aligned}
logTPM &\sim Normal(\mu_i, \sigma) \\
\mu_i & = \alpha + \beta\;time
\end{aligned}
$$
Note: I subtract -2 from the time variable, so the $\alpha$ parameter has the interpretation of the RNA-level at time 2 hrs.


![](./figures/get_estimates-1.png)<!-- -->

Note that to the left 0 (expression) the genes tend to have a positive $\beta$, this maybe just noise.


***
### Session Info


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
##  [1] bindrcpp_0.2.2  brms_2.6.0      Rcpp_1.0.0      broom_0.5.0    
##  [5] zFPKM_1.2.0     forcats_0.3.0   stringr_1.3.1   dplyr_0.7.8    
##  [9] purrr_0.2.5     readr_1.1.1     tidyr_0.8.2     tibble_1.4.2   
## [13] ggplot2_3.1.0   tidyverse_1.2.1
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.3-2            ggridges_0.5.1             
##   [3] rsconnect_0.8.11            rprojroot_1.3-2            
##   [5] markdown_0.8                XVector_0.20.0             
##   [7] GenomicRanges_1.32.7        base64enc_0.1-3            
##   [9] rstudioapi_0.8              rstan_2.18.2               
##  [11] DT_0.5                      mvtnorm_1.0-8              
##  [13] lubridate_1.7.4             xml2_1.2.0                 
##  [15] bridgesampling_0.6-0        knitr_1.20                 
##  [17] shinythemes_1.1.2           bayesplot_1.6.0            
##  [19] jsonlite_1.5                shiny_1.1.0                
##  [21] compiler_3.5.1              httr_1.3.1                 
##  [23] backports_1.1.2             assertthat_0.2.0           
##  [25] Matrix_1.2-14               lazyeval_0.2.1             
##  [27] cli_1.0.1                   later_0.7.5                
##  [29] htmltools_0.3.6             tools_3.5.1                
##  [31] igraph_1.2.2                coda_0.19-2                
##  [33] gtable_0.2.0                glue_1.3.0                 
##  [35] GenomeInfoDbData_1.1.0      reshape2_1.4.3             
##  [37] Biobase_2.40.0              cellranger_1.1.0           
##  [39] nlme_3.1-137                crosstalk_1.0.0            
##  [41] ps_1.1.0                    rvest_0.3.2                
##  [43] mime_0.5                    miniUI_0.1.1.1             
##  [45] gtools_3.8.1                zlibbioc_1.26.0            
##  [47] zoo_1.8-4                   scales_1.0.0               
##  [49] colourpicker_1.0            hms_0.4.2                  
##  [51] promises_1.0.1              Brobdingnag_1.2-6          
##  [53] parallel_3.5.1              SummarizedExperiment_1.10.1
##  [55] inline_0.3.15               shinystan_2.5.0            
##  [57] yaml_2.2.0                  gridExtra_2.3              
##  [59] StanHeaders_2.18.0          loo_2.0.0                  
##  [61] stringi_1.2.4               S4Vectors_0.18.3           
##  [63] dygraphs_1.1.1.6            checkmate_1.8.5            
##  [65] BiocGenerics_0.26.0         pkgbuild_1.0.1             
##  [67] BiocParallel_1.14.2         GenomeInfoDb_1.16.0        
##  [69] rlang_0.3.0.9001            pkgconfig_2.0.2            
##  [71] matrixStats_0.54.0          bitops_1.0-6               
##  [73] evaluate_0.11               lattice_0.20-35            
##  [75] bindr_0.1.1                 labeling_0.3               
##  [77] rstantools_1.5.1            htmlwidgets_1.3            
##  [79] processx_3.2.0              tidyselect_0.2.5           
##  [81] plyr_1.8.4                  magrittr_1.5               
##  [83] R6_2.3.0                    IRanges_2.14.12            
##  [85] DelayedArray_0.6.6          mgcv_1.8-24                
##  [87] pillar_1.3.0.9001           haven_1.1.2                
##  [89] withr_2.1.2                 xts_0.11-2                 
##  [91] abind_1.4-5                 RCurl_1.95-4.11            
##  [93] modelr_0.1.2                crayon_1.3.4               
##  [95] rmarkdown_1.10              grid_3.5.1                 
##  [97] readxl_1.1.0                callr_3.0.0                
##  [99] threejs_0.3.1               digest_0.6.18              
## [101] xtable_1.8-3                httpuv_1.4.5               
## [103] stats4_3.5.1                munsell_0.5.0              
## [105] viridisLite_0.3.0           shinyjs_1.0
```

