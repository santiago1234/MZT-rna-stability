# Identify Regulatory k-mers affecting mRNA stabilit


The goal of this analysis is to develop a model on the residuals after controling
for the optimality and mir-430 pathways.

The analysis **01-Models.Rmd** explores models **02-VisualizeModelResults.Rmd** summarizes the fitted model
to each of the 6-mers. The output is a table with candidate regulatory k-mers.


## Run the test for all 6-mers


The script **test_kmer.R** runs the bayesian model to evaluate each k-mer.

ussage: `Rscript test_kmer.R {}`

This script will output the parameters posterior summary and the model deviance (WAIC). The results are written
to the dir: *results_data/test_kmers/mdl_deviance/*

```bash
head -n 1 ../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_3utr6mer_composition.csv  |\ # get the list of six mers
	tr ',' '\n' |\ # split csv line to line breaks
	tail -n +4 |\ # don show the first 3 lines (column names that are not k-mers)
	parallel -j 32 Rscript test_kmer.R {} # run the script in parallel
```

The tables containing the results for all 6-mer are: *estimated_parameters_6mers.csv* and *waic_6mers.csv*.

Finally the candidates mers are in the table *candidates_6mers.csv*. This table contains the columns:
kmer (all 6-mers) and effect (no effect, destabilazing or stabilizing).
