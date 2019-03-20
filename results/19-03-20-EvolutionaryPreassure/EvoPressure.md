---
title: "Evolutionary Pressure Optimality Analysis"
author: "Santiago Medina"
date: "3/20/2019"
output: 
  html_document:
    keep_md: true
---



**NOTE** in this case I am using the 700 list of m6A targets.




Now we can visualize the data



### The groups

I define the following 3 metagroups.

1. All genes: these genes were detected to be expressed at 2 hrs in the alpha-amanitin data
2. Maternal Genes: 
3. Maternal Genes lower unstable quartile: the 25% most unstable genes of each group

![](./figures/4aOptimalityBoxplot-1.png)<!-- -->

### Herarchical Model to estimate group means


```
## 
## SAMPLING FOR MODEL '1068d9726dcd74e45af0788b195e7b68' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.000605 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 6.05 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 10000 [  0%]  (Warmup)
## Chain 1: Iteration: 1000 / 10000 [ 10%]  (Warmup)
## Chain 1: Iteration: 2000 / 10000 [ 20%]  (Warmup)
## Chain 1: Iteration: 3000 / 10000 [ 30%]  (Warmup)
## Chain 1: Iteration: 3001 / 10000 [ 30%]  (Sampling)
## Chain 1: Iteration: 4000 / 10000 [ 40%]  (Sampling)
## Chain 1: Iteration: 5000 / 10000 [ 50%]  (Sampling)
## Chain 1: Iteration: 6000 / 10000 [ 60%]  (Sampling)
## Chain 1: Iteration: 7000 / 10000 [ 70%]  (Sampling)
## Chain 1: Iteration: 8000 / 10000 [ 80%]  (Sampling)
## Chain 1: Iteration: 9000 / 10000 [ 90%]  (Sampling)
## Chain 1: Iteration: 10000 / 10000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 53.9022 seconds (Warm-up)
## Chain 1:                146.002 seconds (Sampling)
## Chain 1:                199.904 seconds (Total)
## Chain 1: 
## 
## SAMPLING FOR MODEL '1068d9726dcd74e45af0788b195e7b68' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.000185 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.85 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 10000 [  0%]  (Warmup)
## Chain 1: Iteration: 1000 / 10000 [ 10%]  (Warmup)
## Chain 1: Iteration: 2000 / 10000 [ 20%]  (Warmup)
## Chain 1: Iteration: 3000 / 10000 [ 30%]  (Warmup)
## Chain 1: Iteration: 3001 / 10000 [ 30%]  (Sampling)
## Chain 1: Iteration: 4000 / 10000 [ 40%]  (Sampling)
## Chain 1: Iteration: 5000 / 10000 [ 50%]  (Sampling)
## Chain 1: Iteration: 6000 / 10000 [ 60%]  (Sampling)
## Chain 1: Iteration: 7000 / 10000 [ 70%]  (Sampling)
## Chain 1: Iteration: 8000 / 10000 [ 80%]  (Sampling)
## Chain 1: Iteration: 9000 / 10000 [ 90%]  (Sampling)
## Chain 1: Iteration: 10000 / 10000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 11.9766 seconds (Warm-up)
## Chain 1:                23.7998 seconds (Sampling)
## Chain 1:                35.7764 seconds (Total)
## Chain 1: 
## 
## SAMPLING FOR MODEL '1068d9726dcd74e45af0788b195e7b68' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.000101 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.01 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:    1 / 10000 [  0%]  (Warmup)
## Chain 1: Iteration: 1000 / 10000 [ 10%]  (Warmup)
## Chain 1: Iteration: 2000 / 10000 [ 20%]  (Warmup)
## Chain 1: Iteration: 3000 / 10000 [ 30%]  (Warmup)
## Chain 1: Iteration: 3001 / 10000 [ 30%]  (Sampling)
## Chain 1: Iteration: 4000 / 10000 [ 40%]  (Sampling)
## Chain 1: Iteration: 5000 / 10000 [ 50%]  (Sampling)
## Chain 1: Iteration: 6000 / 10000 [ 60%]  (Sampling)
## Chain 1: Iteration: 7000 / 10000 [ 70%]  (Sampling)
## Chain 1: Iteration: 8000 / 10000 [ 80%]  (Sampling)
## Chain 1: Iteration: 9000 / 10000 [ 90%]  (Sampling)
## Chain 1: Iteration: 10000 / 10000 [100%]  (Sampling)
## Chain 1: 
## Chain 1:  Elapsed Time: 5.67766 seconds (Warm-up)
## Chain 1:                13.1993 seconds (Sampling)
## Chain 1:                18.877 seconds (Total)
## Chain 1:
```

![](./figures/4aOptimalityStatistics-1.png)<!-- -->
