# Fig3. Codon optimality regulation overlaps with the microRNA and m6A mRNA destabilizing pathways

*Define the truly m6A targets, rescue or partially rescue in m6A mutant.

Rationale: miR430 and m6A have been proposed as strong regulatory pathways, therefore we want to see if codon optimality overlaps with them and compare the strength.

A)	Boxplot decay for each time point (take top optimal, top non-optimal based on my EMBO code) for All, m6A, miR430 and none miR430 nor m6a (Santi)
B)	Scatter plot for decay vs your optimality for all, m6A, miR430 and none miR430 nor m6a (Santi)
C)	Cumulative plots FOLD change X-2 hin wild, m6A, Dicer, alpha: for all, miR430, m6A and miR430+m6A (Santi)
D)	Take top 500 most unstable mRNAs for each group and then plot the optimality (Santi)
E)	Think in another way to show D.
F)	Predict which one m6a, miR430 or optimality is the stronger. (Santi)
  
Supplementary 2

A)	Cumulative plots FOLD change X-2 hin wild, m6A, Dicer, alpha: for all, miR430, m6A and miR430+m6A (Santi). X, 3,4,5,6,7,8 etc
B)	Potential interesting plot, Xenopus fold change, top 500 and dive in miR430 and not. And plot optimality.

Main conclusion, codon optimality overlaps the other pathways, so we need to think about combinatory effect. Optimal mRNA that needs to go to decay might need miR430 to go down.


## Files and data

**build_data.R |** This script generates the file: *results_data/regulatory_pathways_matrix.csv*. This file contains the following information:

+ *Gene_ID* gene id for the genes I was able to estimate the decay rate.
+ *m6A* is the gene a m6A target?
+ *MiR430* number of MiR430 6-mer sites in 3' UTR.
+ *stable/unstable_mer* number of 6-mer sites in 3' UTR for candidate regulatory elements.
+ *PLS1/2* indicator of codon optimality
+ *is_maternal* is the gene maternal?
