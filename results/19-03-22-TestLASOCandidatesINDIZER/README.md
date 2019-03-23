# Test Marginal Effects in Unseen data

I will use code from a bayesian herarchical model from old analysis.

Create a link to the codeo

```bash
ln -s ../../../rna-stability-model/experiments/18-03-29-TestKmersInRNAseqData/src/testmotifeffect/ testmotifeffect
ln -s ../../../rna-stability-model/experiments/18-03-29-TestKmersInRNAseqData/src/testmotif.py testmotif.py
# add link to the utr and log2 fold change data from old analysis
mkdir results_data
cd results_data
ln -s ../../../../rna-stability-model/experiments/18-03-29-TestKmersInRNAseqData/results/rnasamples_and_utrs.csv rnasamples_and_utrs.csv
cd ..

# Now get a list of the candidate k-mer
cat ../19-03-18-LassoFindKmers/results_data/candidate_elements_laso_analysis.csv  | cut -d',' -f3 | grep -E -v '(PL|vars)' | sort | uniq >candidates_kmers.txt

# run the models in parallel
cat candidates_kmers.txt | parallel -j 31 "python testmotif.py -in results_data/rnasamples_and_utrs.csv -m {}  -out {}-marginal-effects"

mv *marginal-effects results_data/
# remove tmp file with candidates
rm candidates_kmers.txt
```

## Analysis results

The script **plot_in_log2fc_data.R** generates the following plot with the marginal effects:

![Caption](figures/log2fc_val_data_candidates.png)

I save a frame with the effects: *results_data/log2fc_validation_test.csv*
