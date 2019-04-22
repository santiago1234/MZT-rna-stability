---
title: "Overlap MiR-430 by Number of 6-mer sites"
author: "Santiago Medina"
date: "4/11/2019"
output:
   html_document:
    keep_md: true
---




Here, I check the effect of codon optimality (PLS1) in miR-430 targets.

The question is: Do we see an optimality effect in MiR-430 targets genes?

Approach:

I will split MiR-430 targets in 4 categories: 0 sites, 1 site, 2 sites and 3 or more sites.

The response variable is the log2 fold change (X / 3) hrs. 

I have the conditions a-amanitin, polyA and ribo 0.

To access the effect of codon optimality in each of the groups. I fitted the following Hierarchical Model.

Optimality is PLS1 (partial least squares compenent 1).

$$
\begin{aligned}
y &\sim Normal(\mu_i, \sigma_i) \\
\mu_i &= A_i +  B_i * \text{PLS1}\\
A_i &= \alpha + \alpha_{time} + \alpha_{mir} + \alpha_{condition} \\
B_i &= \beta + \beta_{time} + \beta_{mir} + \beta_{condition}
\end{aligned}
$$

### Posterior summary




![](./figures/mir430optimalityEffect-1.png)<!-- -->


