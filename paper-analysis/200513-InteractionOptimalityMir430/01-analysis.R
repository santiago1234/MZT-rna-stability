library(tidyverse)

# xenopus data
expression_xen <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")
pathways_xen <- read_csv("../../data/19-06-12-XenopusData/regulatory_pathways_xenopus.csv") %>% 
  select(-PLS1, -PLS2)

xen_data <- expression_xen %>% 
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
  select(Gene_ID, log2FC, MiR430)


# fish data

fish_fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(
    sample_condition == "wt_ribo",
    is_maternal,
    time == 6
  ) %>% 
  select(-is_maternal, -sample_condition, -time)

mir430 <- 
  read_csv("../../results/19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, MiR430)

fish_data <- inner_join(fish_fc, mir430)

# load optimality (PLS) ---------------------------------------------------

optimality <- read_csv("../../results/19-08-19-OverlapFinal/results_data/pls_species.csv")

datsets <- bind_rows(fish_data, xen_data) %>% 
  rename(gene_id = Gene_ID)
datum <- inner_join(datsets, optimality)

write_csv(datum, "results-data/mir_optimality_fc.csv")
