library(tidyverse)
library(recipes)


# estimate pls components for human data ----------------------------------

datasets <- read_csv("../../data/19-04-29-AllDecayProfiles/decay_profiles_and_seqdata.csv")

humandata <- 
  datasets %>% 
  filter(specie == "human", datatype != "orf-ome") %>%
  select(-specie) %>% 
  select(-`3utr`)

# get the codon composition -----------------------------------------------

orfs <- humandata %>% 
  select(gene_id, coding) %>% 
  distinct()

cchuman <- codonr::codon_composition(orfs, 'coding', 'gene_id') %>%
  select(-contains("N")) %>% 
  rename(gene_id = id_col)


humandata <- humandata %>% 
  select(-coding) %>% 
  inner_join(cchuman)


# get PLS components ------------------------------------------------------


get_pls <- function(dgrp) {
  # extract the pls components for the given data set
  
  feature_extractor <-
    recipe(decay_rate ~ ., data = dgrp)  %>%
    update_role(gene_id, new_role = "id variable") %>%
    step_spatialsign(all_predictors()) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_pls(all_predictors(), outcome = "decay_rate") %>%
    prep(retain = T)
  
  juice(feature_extractor)
  
}

humandata %>% 
  group_by(datatype, cell_type) %>% 
  nest() %>% 
  mutate(plscomp = map(data, get_pls)) %>% 
  unnest(plscomp) %>% 
  write_csv("results_data/pls_components_human.csv")
