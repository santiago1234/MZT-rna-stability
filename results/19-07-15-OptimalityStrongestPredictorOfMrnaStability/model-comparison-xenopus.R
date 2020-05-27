library(tidyverse)
library(brms)

theme_set(theme_tufte(base_family = "Helvetica"))

# gather the data ---------------------------------------------------------
# xenopus data
expression_xen <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

pathways_xen <- read_csv("../../data/19-06-12-XenopusData/regulatory_pathways_xenopus.csv")

dataset_xen <- expression_xen %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  # estimate the log2 fold change
  filter(
    time %in% c(2, 9)
  ) %>% 
  spread(key = time, value = expression_quantification) %>% 
  mutate(
    log2FC = log2(`9` / (`2` + 0.00001))
  ) %>% 
  left_join(pathways_xen) %>% 
  filter(
    !is.na(MiR430),
    in_top_4k_expressed_at_0hrs,
    sample_condition == "wt_ribo"
  ) %>% 
  select(Gene_ID, log2FC, MiR430, PLS1, PLS2)


# model comparison --------------------------------------------------------

fml_opt <- bf(log2FC ~ PLS1 + PLS2)
fml_mir <- bf(log2FC ~ MiR430)

fit_opt <- brm(fml_opt, family = gaussian(), data = dataset_xen, chains = 2, cores = 2)
fit_mir <- brm(fml_mir, family = gaussian(), data = dataset_xen, chains = 2, cores = 2)

loo_opt <- loo(fit_opt)
loo_mir <- loo(fit_mir)

loo_list <- list(optimality=loo_opt, microRNA=loo_mir)

results <- loo_model_weights(loo_list)

tibble(
  model = names(results),
  weights = as.vector(results)
) %>% 
  write_csv("results_data/mdl_weights_xenopus.csv")
