library(tidyverse)
library(recipes)
require(codonr) # my-package

# helper functions --------------------------------------------------------

get_cc <- function(d_train) {
  # gets the codon composition for the orf in d_trian
  d_train %>% 
    select(gene_id, coding) %>% 
    distinct() %>% 
    codonr::codon_composition("coding", "gene_id") %>% 
    rename(gene_id = id_col)
}


pls_decomposition <- function(d_train) {
  # apply the PLS decomposition to the codon
  # space, using the decay rate as the response
  # valdidate
  stopifnot(
    "gene_id" %in% colnames(d_train),
    "decay_rate" %in% colnames(d_train),
    length(colnames(d_train)) == 66
  )
  
  feature_extractor <-
    recipe(decay_rate ~ ., data = d_train)  %>%
    update_role(gene_id, new_role = "id variable") %>%
    step_spatialsign(all_predictors()) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_pls(all_predictors(), outcome = "decay_rate") %>%
    prep(retain = T)
  
  return(feature_extractor)
}


# get the training pipeline data ------------------------------------------
train <- read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/training_set.csv")
test <- read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/validation_set.csv")
train <- bind_rows(train, test)

# apply PLS to MZT --------------------------------------------------------
# here together I apply the preprocessing pipeline to the open readin
# frames for MZT fish and xenopus together a total of 10 profiles
# load the data -----------------------------------------------------------


mztdata <- 
  train %>% 
  filter(cell_type == "embryo mzt") %>% 
  select(-utrlenlog, -cdslenlog, -specie, -cell_type, -datatype)

# compute the codon compostion --------------------------------------------

cc_mzt <- get_cc(mztdata)
mztdata <- 
  mztdata %>% 
  inner_join(cc_mzt) %>% 
  select(-coding)

pls_mzt <- pls_decomposition(mztdata)


# get the orfs ------------------------------------------------------------

fish_orfs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  select(ensembl_gene_id, coding) %>% 
  rename(gene_id = ensembl_gene_id) %>% 
  get_cc()

fish_pls <- bake(pls_mzt, new_data = fish_orfs)

fish_pls <- 
  bind_cols(select(fish_orfs, gene_id), fish_pls) %>% 
  select(-gene_id1) %>% 
  mutate(specie = "fish")

xen_orfs <- read_csv("../../data/19-06-12-XenopusData/xenopus_seq_data_cds_3utr.csv") %>% 
  select(ensembl_gene_id, coding) %>% 
  rename(gene_id = ensembl_gene_id) %>% 
  get_cc()

xen_pls <- bake(pls_mzt, new_data = xen_orfs)

xen_pls <- 
  bind_cols(select(xen_orfs, gene_id), xen_pls) %>% 
  select(-gene_id1) %>% 
  mutate(specie = "xenopus")


# apply pls pipeline to human ---------------------------------------------

humandata <- 
  train %>% 
  filter(specie == "human") %>% 
  select(-utrlenlog, -cdslenlog, -specie, -cell_type, -datatype)

cc_human <- get_cc(humandata)

humandata <- 
  humandata %>% 
  inner_join(cc_human) %>% 
  select(-coding)

pls_human <- pls_decomposition(humandata)

human_orfs <- select(humandata, -decay_rate) %>% 
  distinct()

human_pls <- bake(pls_human, new_data = human_orfs) %>% 
  mutate(specie = "human")

# apply to mouse ----------------------------------------------------------

mousedata <- 
  train %>% 
  filter(specie == "mouse") %>% 
  select(-utrlenlog, -cdslenlog, -specie, -cell_type, -datatype)

cc_mouse <- get_cc(mousedata)

mousedata <- 
  mousedata %>% 
  inner_join(cc_mouse) %>% 
  select(-coding)

pls_mouse <- pls_decomposition(mousedata)

mouse_orfs <- select(mousedata, -decay_rate) %>% 
  distinct()

mouse_pls <- bake(pls_mouse, new_data = mouse_orfs) %>% 
  mutate(specie = "mouse")


# save the pls components -------------------------------------------------

bind_rows(fish_pls, xen_pls, human_pls, mouse_pls) %>% 
  write_csv("results_data/pls_species.csv")
