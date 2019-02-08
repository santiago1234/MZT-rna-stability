# Identify Regulatory k-mers affecting mRNA stabilit


## Run the test for all 6-mers

```bash
head -n 1 ../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_3utr6mer_composition.csv  |\ # get the list of six mers
	tr ',' '\n' |\ # split csv line to line breaks
	tail -n +4 |\ # don show the first 3 lines (column names that are not k-mers)
	parallel -j 32 Rscript test_kmer.R {} # run the script in parallel
```
