---
title: "Compute CSC"
author: "Santiago Medina"
date: "1/17/2019"
output: 
  html_document:
    keep_md: true
---

## CSC analysis

Here, I estimate the CSC as defined in [Prensyak et al 2016](https://www.cell.com/cell/fulltext/S0092-8674(15)00195-6) for the a-amanitin time course.




I filter the data $\alpha > 0$, bassically drop low expressed genes.




### Data Exploration

![](./figures/vis-1.png)<!-- -->

There is an unpublished a-amanitin data set from Ariel's previous lab. I will load this data and see if they correlate.


![](./figures/cor_with_old_data-1.png)<!-- -->

GOOD!, the correlation is not bad.

## Estimate CSC values

![](./figures/csc-1.png)<!-- -->

#### Do these estimates correlate with the EMBO 2016 paper estimates?


![](./figures/csc_cors-1.png)<!-- -->

