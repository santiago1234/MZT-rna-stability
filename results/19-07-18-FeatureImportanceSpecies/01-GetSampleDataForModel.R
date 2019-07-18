library(tidyverse)


# get the data and create sample data for analysis ------------------------
# for each specie I take a sample of 10K genes or in the number of genes
# is less than 10K then I take all of them

dset <- read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/training_set.csv")

# replace utr leng log NA with median

dset <- replace_na(dset, list(utrlenlog = mean(dset$utrlenlog, na.rm = T)))

sample_genes <- function(grp) {
  
  if (nrow(grp) < 10000) {
    grp
  }
  else {
    sample_n(grp, 10000)
  }
}


set.seed(134)
dset_sample <- 
  dset %>% 
  group_by(specie) %>% 
  nest() %>% 
  mutate(sampledata = map(data, sample_genes)) %>% 
  unnest(sampledata)

## compute the codon composition

orfs <- select(dset_sample, gene_id, coding) %>% distinct()

codon_composition <- codonr::codon_composition(orfs, orf_col = "coding", id_col = "gene_id")
codon_composition <- rename(codon_composition, gene_id = id_col)

inner_join(dset_sample, codon_composition) %>% 
  write_csv("results_data/sample_data_for_feature_importance_analysis.csv")
