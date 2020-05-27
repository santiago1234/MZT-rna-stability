library(tidyverse)

# built predictor matrix --------------------------------------------------

ids_mappings <- read_tsv("../19-02-26-M6A-targets/m6A-data/mart_export.txt") %>% 
  rename(Gene_ID = `Gene stable ID`, `tracking ID` = `Gene name`) %>% 
  select(Gene_ID, `tracking ID`) %>% 
  unique()

gene_grps <- read_csv("../19-02-27-m6aPathway/results_data/m6a_paper_data/gene_classification_accordinm6a_paper.csv")

mappings_grps <- tribble(
  ~gene_grp, ~class_grp,
  0, "semi stable 1",
  1, "semi stable 2",
  5, "maternal late",
  2, "maternal early",
  3, "zygotic late",
  4, "zygotic early"
)

gene_grps <- inner_join(gene_grps, ids_mappings) %>% 
  select(Gene_ID, `gene group`) %>% 
  rename(gene_grp = `gene group`) %>% 
  inner_join(mappings_grps) %>% 
  select(-gene_grp)

predictors <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  select(Gene_ID, m6A, MiR430, PLS1, PLS2) %>% 
  inner_join(gene_grps)


write_csv(predictors, "results_data/predictors.csv")
