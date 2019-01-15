---
title: "Notebook"
author: "Santiago Medina"
date: "1/9/2019"
output:
   html_document:
    keep_md: true
---

In this notebook I record my daily progress in detail. This notebook should record
my observations.

***
### 19-01-09

+ Created directory organization based on [Noble 2009](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)
+ Started version control with git. I host a private [github repo](https://github.com/santiago1234/MZT-rna-stability)
+ Downloaded mRNA seq time courses generated by the lab.

I have to re-run previews analysis and do some quality checks on this data. Ariel suggested to do: PCA, poly-A vs ribo (also the ratio) and then the model.

***
### 19-01-10

+ Run Exploratory analysis for RNA-seq time course (PCA), the PCA analysis looks find.

I also generated the following bi-variate distribution: ![](../../results/19-01-10-EDA-RNA-time-course/figures/bivariate_polyA-1.png)

- For the alpha-amanitin-time course I have to see what is the best way to perform the SPIKE-in normalization, I can
  ask Ariel or Comp-Bio. Also the following paper by [Risso](https://www.nature.com/articles/nbt.2931) can be useful to read.
