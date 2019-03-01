## Test k-mer sequences in human decay data


### Can we detect a k-mer effect in endogenous human data?

Here I took the following approach:
1. Get the decay rate profiles for human genes (slam-seq and endogenous data)
2. Get the 3' UTR sequences for these human genes
3. For the candidates, k-mers detect if the k-mer sequence is present in the 3' UTR (two groups genes with k-mers and gene with no k-mers)
4. Is there a difference in the decay rate of the genes with k-mer vs. genes without k-mer?

### Conclusion

no significant difference, therefore candidates do not affect human genes in the cell types I tested.
