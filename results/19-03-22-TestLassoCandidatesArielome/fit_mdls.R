#!/usr/bin/env Rscript

# fit model to kmer -------------------------------------------------------
library(tidyverse)
library(brms)
library(parallel)


# load data ---------------------------------------------------------------

utrs <- read_csv("results_data/collapsed_3utrs_every100.csv")
# for now remove the +MO condition
utrs <- 
  utrs %>% 
  filter(condition == "wt") %>% # we drop the treated (+MO) samples, we are no interested in
  mutate(
    latetimepoint = time == 8
    )

# fit mdl -----------------------------------------------------------------

fit_kmer <- function(kmer) {
  utrs <- mutate(utrs, kmer = str_count(collapsed_utrs, kmer))
  # fit bayesian mdl
  fit <- brm(formula = kmer ~ latetimepoint, data = utrs, family = poisson(), chains = 2, cores = 2)
  
  # save posterior summary
  outname <- str_c("results_data/", 'mdlfit_', kmer, ".csv")
 
   posterior_summary(fit) %>% 
    as_tibble(rownames = "parameter") %>% 
    filter(parameter != "lp__") %>% 
    mutate(kmer = kmer) %>% 
    write_csv(
      outname
    )
}


# get k-mers to test ------------------------------------------------------
mir430mers <- c(
  "AGCACT",
  "GCACTT",
  "GCACTTA",
  "AGCACTT",
  "AGCACTTA"
)

candidates <- read_csv("../19-03-18-LassoFindKmers/results_data/candidate_elements_laso_analysis.csv") %>% 
  filter(!str_detect(vars, "PLS")) %>% 
  pull(vars) %>% 
  unique() %>% 
  c(.,mir430mers) # append mir 6,7,8,9 mers for validation


# run mdls in paralell ----------------------------------------------------

mclapply(
  candidates,
  FUN = fit_kmer,
  mc.cores = 1
)