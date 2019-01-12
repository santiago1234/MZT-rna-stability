Decay Rate Estimation
================
Santiago Medina
1/11/2019

### Overview:

The goal of this analysis is to estimate the decay rate for the endogenous genes with the a-amanitin data.

### 1. Tidy data for analysis

I use the a-amanitin data set. This includes 12 time points.

### 2. Detect active genes

To estimate the decay rate I will use only active genes. I will use the [zFPKM](https://bioconductor.org/packages/release/bioc/html/zFPKM.html) package for this.

Plot of the distribution of gene expression.

![](./figures/expression_distribution-1.png)

as the authors suggest a cut-off of zFPKM &gt; -3 is used.

![](./figures/active_genes-1.png)

### 3. Estimate Decay Rate

Can you go full Bayes?

![](./figures/unnamed-chunk-1-1.png)

------------------------------------------------------------------------

### Session Info

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
    ##  [1] bindrcpp_0.2.2  zFPKM_1.2.0     forcats_0.3.0   stringr_1.3.1  
    ##  [5] dplyr_0.7.8     purrr_0.2.5     readr_1.1.1     tidyr_0.8.2    
    ##  [9] tibble_1.4.2    ggplot2_3.1.0   tidyverse_1.2.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.0                  lubridate_1.7.4            
    ##  [3] lattice_0.20-35             assertthat_0.2.0           
    ##  [5] rprojroot_1.3-2             digest_0.6.18              
    ##  [7] R6_2.3.0                    GenomeInfoDb_1.16.0        
    ##  [9] cellranger_1.1.0            plyr_1.8.4                 
    ## [11] backports_1.1.2             stats4_3.5.1               
    ## [13] evaluate_0.11               httr_1.3.1                 
    ## [15] pillar_1.3.0.9001           zlibbioc_1.26.0            
    ## [17] rlang_0.3.0.9001            lazyeval_0.2.1             
    ## [19] readxl_1.1.0                rstudioapi_0.8             
    ## [21] S4Vectors_0.18.3            Matrix_1.2-14              
    ## [23] checkmate_1.8.5             rmarkdown_1.10             
    ## [25] labeling_0.3                BiocParallel_1.14.2        
    ## [27] RCurl_1.95-4.11             munsell_0.5.0              
    ## [29] DelayedArray_0.6.6          broom_0.5.0                
    ## [31] compiler_3.5.1              modelr_0.1.2               
    ## [33] pkgconfig_2.0.2             BiocGenerics_0.26.0        
    ## [35] htmltools_0.3.6             tidyselect_0.2.5           
    ## [37] SummarizedExperiment_1.10.1 GenomeInfoDbData_1.1.0     
    ## [39] matrixStats_0.54.0          IRanges_2.14.12            
    ## [41] viridisLite_0.3.0           crayon_1.3.4               
    ## [43] withr_2.1.2                 bitops_1.0-6               
    ## [45] grid_3.5.1                  nlme_3.1-137               
    ## [47] jsonlite_1.5                gtable_0.2.0               
    ## [49] magrittr_1.5                scales_1.0.0               
    ## [51] cli_1.0.1                   stringi_1.2.4              
    ## [53] XVector_0.20.0              xml2_1.2.0                 
    ## [55] tools_3.5.1                 Biobase_2.40.0             
    ## [57] glue_1.3.0                  hms_0.4.2                  
    ## [59] parallel_3.5.1              yaml_2.2.0                 
    ## [61] colorspace_1.3-2            GenomicRanges_1.32.7       
    ## [63] rvest_0.3.2                 knitr_1.20                 
    ## [65] bindr_0.1.1                 haven_1.1.2
