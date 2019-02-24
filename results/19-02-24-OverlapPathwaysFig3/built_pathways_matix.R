library(tidyverse)

# this script produces a matrix with information for each gene about the 
# regulatory pathways
# load the m6a data (m6A targets) -----------------------------------------

m6A <- 
  read_csv("../../../../projectos/rna-stability-model/paper/figures/f_02/data/m6A_targets.csv")
idsmapping <- 
  read_tsv("../../../../projectos/rna-stability-model/paper/figures/f_02/data/mart_export.txt") %>% 
  rename(gene_id = `Gene stable ID`, tracking.ID = `Gene name`)
m6A <- 
  inner_join(m6A, idsmapping, by = "tracking.ID") %>% 
  select(gene_id) %>% 
  unique() %>% 
  mutate(m6A = T) %>% 
  rename(Gene_ID = gene_id)

utrs <- 
  read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  select(-coding) %>% 
  mutate(
    MiR430 = str_count(`3utr`, "GCACTT"),
    stable_TCGGCG = str_count(`3utr`, "TCGGCG"),
    stable_CTCCCG = str_count(`3utr`, "CTCCCG"),
    stable_CCTGGG = str_count(`3utr`, "CCTGGG"),
    unstable_AGACGG = str_count(`3utr`, "AGACGG"),
    unstable_TCCGTA = str_count(`3utr`, "TCCGTA"),
  ) %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  select(-`3utr`)

optimality <- 
  read_csv("../19-01-18-ObtainEstimateOfCodonOptimalityPLS/results_data/pls_components_fish_genes.csv")

# I will add an indicator in case i only want maternal genes in the analysis
maternal <- 
  read_tsv("../../data/19-01-17-EMBO2016DATA/datasets/Half_life_Zebrafish.txt") %>% 
  mutate(is_maternal = T) %>% 
  rename(Gene_ID = gene) %>% 
  select(-HL)


# make matrix with regulatory pathwyas ------------------------------------

reg_pathways <- 
  list(m6A, utrs, optimality, maternal) %>% 
  reduce(function(x, y) full_join(x, y, by = "Gene_ID")) %>% 
  filter(!is.na(decay_rate)) %>%
  replace_na(list(m6A = F, is_maternal = F))

write_csv(reg_pathways, "results_data/regulatory_pathways_matrix.csv")
