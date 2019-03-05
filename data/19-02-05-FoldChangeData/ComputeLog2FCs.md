---
title: "Generate log2 Fold Change Data"
author: "Santiago Medina"
date: "2/5/2019"
output: 
  html_document:
    keep_md: true
---



The goal of this notebook is to generate the log2 Fold Change data, this data will be usefull in other analysis.





![](./figures/computelog2FC-1.png)<!-- -->


I will ad an indicator column if the gene is maternal. I have a list of genes that are supposly maternal (*../19-01-17-EMBO2016DATA/datasets/Half_life_Zebrafish.txt*) although I am not 100% sure.

![](./figures/maternals-1.png)<!-- -->

Save the data, I will also save a tidy version of the time course data in case i want to see effect of regulatory sequences at the RNA level.




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
##  [1] bindrcpp_0.2.2  forcats_0.3.0   stringr_1.4.0   dplyr_0.7.8    
##  [5] purrr_0.3.0     readr_1.3.1     tidyr_0.8.2     tibble_2.0.1   
##  [9] ggplot2_3.1.0   tidyverse_1.2.1
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.0        cellranger_1.1.0  pillar_1.3.1.9000
##  [4] compiler_3.5.1    plyr_1.8.4        bindr_0.1.1      
##  [7] tools_3.5.1       digest_0.6.18     lubridate_1.7.4  
## [10] jsonlite_1.6      evaluate_0.11     nlme_3.1-137     
## [13] gtable_0.2.0      lattice_0.20-35   pkgconfig_2.0.2  
## [16] rlang_0.3.1.9000  cli_1.0.1         rstudioapi_0.9.0 
## [19] yaml_2.2.0        haven_2.1.0       withr_2.1.2      
## [22] xml2_1.2.0        httr_1.3.1        knitr_1.20       
## [25] hms_0.4.2         rprojroot_1.3-2   grid_3.5.1       
## [28] tidyselect_0.2.5  glue_1.3.0        R6_2.3.0         
## [31] readxl_1.1.0      rmarkdown_1.10    reshape2_1.4.3   
## [34] modelr_0.1.2      magrittr_1.5      ggthemes_4.0.1   
## [37] backports_1.1.3   scales_1.0.0      htmltools_0.3.6  
## [40] rvest_0.3.2       assertthat_0.2.0  colorspace_1.4-0 
## [43] labeling_0.3      stringi_1.2.4     lazyeval_0.2.1   
## [46] munsell_0.5.0     broom_0.5.0       crayon_1.3.4
```





