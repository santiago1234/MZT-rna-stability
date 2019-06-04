---
title: "Exploratory Analysis"
author: "Santiago Medina"
date: "4/30/2019"
output: 
  html_document:
    keep_md: true
---

## Overview

Here, Before evaluation machine learning models, I will explore the data set. That
we are using. The decay rate may be in different scales hence I will scale it according to the data it was originated. Also
here I will add some new features (3' utr length, cds length)



***

## remove outliers

I will drop the outliers since they can affect the predictive models.

I wont use the orf-ome for model training. It will be drop it.

![](./figures/datawithoutliers-1.png)<!-- -->

## Does decay rate shows correlation bewtween cell types?

![](./figures/unnamed-chunk-1-1.png)<!-- -->

![](./figures/datanooutliers-1.png)<!-- -->

***

### Scale Decay Rate

![](./figures/decay-1.png)<!-- -->

It is easy to check that the data comes in different scales, also there are outliers present in the data.

I will save a data with the mean and variance of each data set, this can be used later in case I need
to convert the data back to the original scale.

![](./figures/decayscaled-1.png)<!-- -->

### number of data points

![](./figures/npoints-1.png)<!-- -->


***

### Add features: cds length and 3utr length

I will add the length of the cds and also the length of the 3' UTR. After this, I will drop the 3utr.





### Validation Set

Now, I will define a validation set that will be only used for evaluating the predictive models.




