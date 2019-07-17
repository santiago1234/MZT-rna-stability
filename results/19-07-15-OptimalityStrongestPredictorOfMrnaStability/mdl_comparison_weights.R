library(tidyverse)
library(brms)

dset <- read_csv("results_data/pathways_mRNAstability.csv")
dset <- dset[!apply(dset, 1, function(x) any(is.na(x))), ,]


mdl_weights <- function(dset) {
  # bayesian model comparison analysis
  fml_opt <- bf(stability ~ PLS_1 + PLS_2 + PLS_3 + PLS_4 + PLS_5 + PLS_6 + PLS_7 + PLS_8 + PLS_9)
  fml_m6a <- bf(stability ~ m6A)
  fml_mir <- bf(stability ~ microRNAsites)
  
  fit_opt <- brm(fml_opt, family = gaussian(), data = dset, chains = 2, cores = 2)
  fit_m6a <- brm(fml_m6a, family = gaussian(), data = dset, chains = 2, cores = 2)
  fit_mir <- brm(fml_mir, family = gaussian(), data = dset, chains = 2, cores = 2)
  
  loo_opt <- loo(fit_opt)
  loo_m6a <- loo(fit_m6a)
  loo_mir <- loo(fit_mir)
  
  loo_list <- list(optimality=loo_opt, m6A=loo_m6a, microRNA=loo_mir)
  
  results <- loo_model_weights(loo_list)
  
  tibble(
    model = names(results),
    weights = as.vector(results)
  )
  
}

mdlwts <- dset %>% 
  group_by(specie, cell_type) %>% 
  nest() %>% 
  mutate(
    model_weights = map(data, mdl_weights)
  )


unnest(mdlwts, model_weights) %>% 
  write_csv("results_data/mdl_weights.csv")
