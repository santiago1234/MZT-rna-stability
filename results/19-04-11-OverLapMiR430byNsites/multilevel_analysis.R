library(tidyverse)
library(brms)
library(ggridges)

# load the data -----------------------------------------------------------

predictors <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv")

# Ariel suggested to control for the UTR length let's load the data
# I will take the log of the UTR length

utrlen <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  mutate(utrlenlog = log(str_length(`3utr`) + 1)) %>% # a pesudo count is added to the log
  select(ensembl_gene_id, utrlenlog) %>% 
  rename(Gene_ID = ensembl_gene_id)

# merge the data sets
predictors <- inner_join(utrlen, predictors)

# load the response data (log2FC)
log2fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(time > 3) %>% 
  filter(is_maternal, time %% 1 == 0) %>% 
  select(-is_maternal) %>% 
  filter(!is.na(log2FC), !is.infinite(log2FC))

# tidy the data for modeling ----------------------------------------------

dta <- 
  predictors %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, MiR430, PLS1, utrlenlog) %>% 
  inner_join(log2fc)

# recode the MiR430 targets

code_mir_targets <- function(x) {
  if (x >= 3) return("3 or more sites")
  else paste0(x, " sites")
}

dta <- dta %>% 
  mutate(
    MiR430 = map_chr(MiR430, code_mir_targets)
  )

# save a copy in the data for visualization

write_csv(dta, "results_data/mdl_dta.csv")

# Scale log2fc
# the variance of the data increases with time, this is because as time goes
# on the TPMs look more different

dta <- dta %>% 
  group_by(time, sample_condition) %>% 
  mutate(log2FC = scale(log2FC) %>% as.numeric()) %>% 
  ungroup()

# make time factor to use as a level

dta <- mutate(dta, time = as.factor(time))

## THE MODEL FORMULA
## we assume the utr effect is constant over time
fml <- bf(log2FC ~ 1 + PLS1 + utrlenlog +  (1 + PLS1 | MiR430) + (1 + PLS1 | sample_condition) + (1 + PLS1 | time))

fit <- brm(fml, family = gaussian(), data = dta, chains = 1, cores = 1)

## save multivariate posterior samples
## and mdl and posterior samples
write_rds(fit, "results_data/multilevel_mdl.rds")
posterior_samples(fit) %>% 
  as_tibble() %>% 
  write_csv("results_data/posterior_samples.csv")
