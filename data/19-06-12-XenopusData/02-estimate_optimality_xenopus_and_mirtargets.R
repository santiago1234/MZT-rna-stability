# estimate the optimality of xenopus bassed on the fish optimality
# PLS projection on xenopus genes
library(tidyverse)
library(recipes)

fish_preprocessing <- read_rds("../../results/19-01-18-ObtainEstimateOfCodonOptimalityPLS/results_data/extrac_pls_components.rds")
  

# get the codon composition of xenopus ------------------------------------

xenopus_data <- read_csv("xenopus_seq_data_cds_3utr.csv") 
# get the codon composition
cc_xen <- codonr::codon_composition(xenopus_data, "coding", "ensembl_gene_id") %>% 
  rename(Gene_ID = id_col)

optimality <- bake(fish_preprocessing, new_data = mutate(cc_xen, Gene_ID = NA)) %>% 
  select(-Gene_ID) %>% 
  bind_cols(select(cc_xen, Gene_ID))
  
# count the number of mir-430 six mer sites 3’ UTR ------------------------

reg_pathways <- 
  xenopus_data %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  mutate(
    MiR430 = str_count(`3utr`, "GCACTT"),
  ) %>% 
  inner_join(optimality) %>% 
  select(-coding, -`3utr`)

# define the top 4K genes expressed at time 0hrs --------------------------


expression <- read_csv("time_course_xenopus.csv") %>% 
  filter(time == 0)

top_4k_expressed <- expression %>% 
  group_by(ensembl_gene_id) %>% 
  summarise(expression = mean(expression_quantification)) %>% 
  arrange(-expression) %>% 
  pull(ensembl_gene_id) %>% 
  .[1:4000]

reg_pathways %>% 
  mutate(
    in_top_4k_expressed_at_0hrs = Gene_ID %in% top_4k_expressed
  ) %>% 
  write_csv("regulatory_pathways_xenopus.csv")
