#!/usr/bin/env Rscript
# usage: Rsciprt script.R <KMER>

library(tidyverse)
library(brms)

scaler <- function(x) (x - mean(x)) / sd(x)

fc_tc <- "../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv"
pls <- "../../results/19-01-18-ObtainEstimateOfCodonOptimalityPLS/results_data/pls_components_fish_genes.csv"
utrs <- "../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv"

# helper functions --------------------------------------------------------

data_loader <- function(sample_cond, kmer) {
  codonr::make_mdl_frame(
    log2FCtimecourse_dp = fc_tc,
    pls_optimality_dp = pls,
    utrs_dp = utrs,
    kmers_to_keep = c("GCACTT", kmer), 
    .sample_condition = sample_cond,
    minimum_time_point = 3,
    .maternal = T
  ) %>% 
    mutate_if(is.numeric, scaler)
}
 
## model formula

make_fml <- function(kmer) {
  if (kmer == "GCACTT") {
    bform <- bf(
      log2FC ~ I(time^2) + time*PLS1 + time*PLS2 + time*GCACTT,
      sigma ~ poly(time, 2)
    )
  } else {
    fml_mean <- str_c("log2FC ~ I(time^2) + time*PLS1 + time*PLS2 + time*GCACTT + time*", kmer)
    bform <- bf(
      fml_mean,
      sigma ~ poly(time, 2)
    )
  }
  return(bform)
}

run_mdl <- function(fml, data) {
  brm(
    formula = fml,
    data = data,
    family = gaussian(),
    chains = 1,
    cores = 1,
    prior = prior(lasso(), class = "b")
  )
}

# model metrics
get_waic <- function(mdl) {
  waic(mdl) %>% 
    (function(x) x$estimates)() %>% 
    as_tibble(rownames = "val")
}


get_posterior_summary <- function(mdl) {
  posterior_summary(mdl) %>% 
    as_tibble(rownames = "parameter")
}


fit_model_to_data <- function(sample_name, kmer) {
  ## get the data
  datam <- data_loader(sample_name, kmer)
  model <- run_mdl(make_fml(kmer), datam)
  ## summarize
  waic_res <- get_waic(model)
  params <- get_posterior_summary(model)
  list(
    waic = mutate(waic_res, sample_name = sample_name, kmer = kmer),
    params = mutate(params,  sample_name = sample_name, kmer = kmer)
    )
}


run_analysis <- function(kmer) {
  wt_poly <- fit_model_to_data("wt_polya", kmer)
  wt_ribo <- fit_model_to_data("wt_ribo", kmer)
  aamanitin_polya <- fit_model_to_data("aamanitin_polya", kmer)
  
  ## gather results in a tible
  params_res <- bind_rows(
    wt_poly$params,
    wt_ribo$params,
    aamanitin_polya$params
  )
  
  waic_res <- bind_rows(
    wt_poly$waic,
    wt_ribo$waic,
    aamanitin_polya$waic
  )
  write_csv(waic_res, str_c("results_data/test_kmers/mdl_deviance/", kmer, "_waic.csv"))
  write_csv(params_res, str_c("results_data/test_kmers/mdl_deviance/", kmer, "_params.csv"))
}



# run analysis for use args -----------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
kmer <- args[1]
cat("runing analysis ...", kmer, "\n")
run_analysis(kmer)