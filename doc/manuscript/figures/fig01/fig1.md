---
title: "Codon composition strongly predicts mRNA stability"
author: "Santiago Medina"
date: "March 1, 2019"
output: 
  html_document:
    keep_md: true
---



##	Half-life or decay distribution (Santi)

![](./figures/f01_decaydistribution-1.png)<!-- -->

## C)	Model to predict mRNA decay in alpha (scatter plot, R and P value), How much can you explain (Santi)


Best tuning parameters

```
##        sigma    C
## 1 0.01000694 0.25
```

![](./figures/f01_mdl_test_data_predictions-1.png)<!-- -->


D)	Calculate the mRNA decay of the 1nt out of frame (codon vs nt) and Gopal can inject those. (Santi and Gopal)

TODO: tibble with reporter sequences
I do not think that the reporter sequences I am using from the EMBO paper are correct.


![](./figures/f01_reporters_predictions-1.png)<!-- -->
