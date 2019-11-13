library(tidyverse)
library(parallel)

# get wild-type stability -------------------------------------------------
# stability is defined as the log2 fold change 6 / 3 hrs
# important parameters: for fish, time point and sample condition (ribo orpolyA)
# which ones gives better results
stability_mzt_fish <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(sample_condition == "wt_ribo", time == 6, is_maternal) %>% 
  select(-time, -sample_condition, -is_maternal) %>% 
  rename(gene_id = Gene_ID) %>% 
  mutate(
    specie = "fish"
  )

# xenous stability

xen_rna <- read_csv("../../data/19-06-12-XenopusData/time_course_xenopus.csv")

xen_rna <-
  xen_rna %>%
  rename(Gene_ID = ensembl_gene_id) %>%
  filter(
    time %in% c(2, 9),
    sample_condition == "wt_polya"
  ) %>%
  spread(key = time, value = expression_quantification) %>%
  arrange(-`2`) %>%
  slice(1:6000) %>%
  mutate(
    log2FC = log2(`9` / (`2` + 0.00001)),
    specie = "xenopus"
  ) %>%
  select(Gene_ID, log2FC, specie) %>% 
  rename(gene_id = Gene_ID) %>% 
  mutate(specie = "xenopus")


predicted_stability <- read_csv("results-data/predicted_stabilit_fish_and_xenopus_all_genes.csv") %>% 
  select(-coding, -cell_type, -datatype, -utrlenlog)


# add the 3’ UTR to color code the residual plot --------------------------


stability_mzt_fish <- 
  read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  select(ensembl_gene_id, `3utr`, coding) %>% 
  rename(gene_id = ensembl_gene_id) %>% 
  right_join(stability_mzt_fish)

xen_rna <- 
  read_csv("../../data/19-06-12-XenopusData/xenopus_seq_data_cds_3utr.csv") %>% 
  select(ensembl_gene_id, `3utr`, coding) %>% 
  rename(gene_id = ensembl_gene_id) %>% 
  inner_join(xen_rna)

# compute the % of optimal codons -----------------------------------------

opt_ratio <- function(orf) {
  # get the percentage of optimal codons according to bazzini embo2016
  codonr::count_codons(orf) %>%
    full_join(codonr::optimality_code_embo2016) %>%
    replace_na(list(n = 0)) %>%
    group_by(optimality) %>%
    summarise(suma = sum(n)) %>%
    spread(key = optimality, value = suma) %>%
    mutate(
      ratio_log = log2(optimal / (`non-optimal` + 1))
    ) %>%
    pull(ratio_log)
}

mzt_stability <- bind_rows(stability_mzt_fish, xen_rna) %>% 
  inner_join(predicted_stability)

mzt_stability$optimality_ratio <- mclapply(mzt_stability$coding, FUN = opt_ratio, mc.cores = 4) %>% unlist()

write_csv(mzt_stability, "results-data/mzt_data_predictionsAndObservedStability.csv")
