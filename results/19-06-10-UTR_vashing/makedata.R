library(tidyverse)
library(brms)

predictors <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, m6A, MiR430, PLS1, PLS2)

# add utrs to include ARE elements ----------------------------------------

utrs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv")

predictors <- utrs %>% 
  mutate(
    ARE = str_count(`3utr`, pattern = "TATTTA"),
    MiR_4302nd6mer = str_count(`3utr`, pattern = "AGCACT"),
    utr_log = log(str_length(`3utr`) + 1)
  ) %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  select(Gene_ID, ARE, MiR_4302nd6mer, utr_log) %>% 
  inner_join(predictors)

log2fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv")


mdldata <- 
  log2fc %>% 
  filter(time == 6, sample_condition != "aamanitin_polya") %>% 
  inner_join(predictors)

mdldata2 <- mdldata %>% 
  spread(key = sample_condition, value = log2FC)

write_csv(mdldata2, "results_data/mdl_data.csv")
