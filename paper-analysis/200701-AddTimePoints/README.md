## Add time points

_1. The authors utilized a data set by Owens et al., which has samples in 30-minute intervals during Xenopus early embryogenesis. The authors chose 1 and 9 hpf to calculate mRNA stability. It would be helpful to get a sense of the robustness of these findings by using more time points._

Here, I add more time points for this reviewer comment.

The first script *01-PrefictFoldChanges.R* estimates the log2fold change based on the predicted stability (linear transformation).


To run the analysis:

```bash
snakemake figures/prediction-by-timepoint.pdf
```

Which time points I added?

- Time points in fish: 4 hrs to 8 hrs
- Time points in xenopus: 3 hrs to 12 hrs
