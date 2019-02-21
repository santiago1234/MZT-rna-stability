to get marginal effects in a cross validation seeting (other data that we did not use for training, to get marginal effects
of k-mers seqs I ran the following command in the following directory)

```bash
cut -f3 -d','  ../../../../190108-mzt-rna-stability/results/19-02-12-TopKmersFeatureSelection/results_data/top_candidates.csv  | sort  | uniq | grep -v -E "(PLS|kmer)"  | parallel -j 4 "python testmotif.py -in  ../results/rnasamples_and_utrs.csv -m {} -out ../../../../190108-mzt-rna-stability/results/19-02-12-TopKmersFeatureSelection/results_data/marginal-effec-crossval-{}"
#in dir:
# /home/smedina/my-projects/rna-stability-model/experiments/18-03-29-TestKmersInRNAseqData/src
```


