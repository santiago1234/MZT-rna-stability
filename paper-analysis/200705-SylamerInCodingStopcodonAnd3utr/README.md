# Run sylamer in the coding sequence, a window around the stop codon, and  the 3' UTR

This is a generalization of the analysis presented in: ../200702-SylamerAddTimePoints/ but here I
also include the coding region and a window around the stop codon.

The point of this analysis is to see if we can detect enrichment of the m6A motif. We did
not find significant enrichment in the 3' UTR. The m6a motif according to the literatur should be
depleated in the 3' UTR. Hence my goal here is to see whether i can detect this enrichment in other regions of the transcript.

Conclusions:
I failed to detect the m6A motif.


Run the pipeline.

```bash
snakemake --cores 3 figures/syl-coding-stop-utr.pdf
```

Next I decided to plot the residual as a function of the number of m6A sites (GGACT) in the 3' UTR and the coding sequence. I got the expected result. 

- We can see the repressive effect of m6A
- We can see the effect is stronger in the coding compared to the 3' UTR sequence.

The output plot the number 4 represents 4 or more sites.
