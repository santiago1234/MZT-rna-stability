# Outliers List
# The outliers genes were defined based on a linear model to predict log2 fold 
# change (6 / 2 hrs), this linear model included as predictors the m6A targets, 
# miR430 sites, optimality (PLS1 and PLS2) and ARE sites. The outliers were 
# defined based on the residual value (IQR test). I select the final outliers 
# to be outliers in WT-ribo0, WT-polyA, and DICER data. Also, I filter for 
#genes that are not m6a targets and do not have miR-430 sites and not ARE sites.
library(tidyverse)
library(brms)


# get the predictors ------------------------------------------------------

predictors <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(is_maternal) %>% 
  select(Gene_ID, m6A, MiR430, PLS1, PLS2)

# add utrs to include ARE elements ----------------------------------------

utrs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv")

predictors <- utrs %>% 
  mutate(
    ARE = str_count(`3utr`, pattern = "TATTTA"),
    MiR_4302nd6mer = str_count(`3utr`, pattern = "AGCACT"),
    utr_log = log(str_length(`3utr`) + 1)
  ) %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  select(Gene_ID, ARE, MiR_4302nd6mer, utr_log) %>% 
  inner_join(predictors)




log2fc <- read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv")

# get dicer fold change ---------------------------------------------------

# the data is poly a
# get the mrna level at 2 hrs and take the fold change
rna_2hrs_polt <- read_csv("../../data/19-01-09-RNAseqProfilesFish/rna-seq-profiles/polyA-profile.csv") %>% 
  select(Gene_ID, s_wt_2h)

dicer <- read_csv("../../../181029-Paper/data/rna-time-course/dicer-data/cufflinks_fpkm_all.csv") %>% 
  select(GeneID, s_MZdicer_6h) %>% 
  rename(Gene_ID = GeneID)

fc_dicer <- inner_join(rna_2hrs_polt, dicer) %>% 
  mutate(log2FC = log2(s_MZdicer_6h / s_wt_2h)) %>% 
  select(-s_wt_2h, -s_MZdicer_6h) %>% 
  mutate(sample_condition = "dizer", time = 6)



mdldata <- 
  log2fc %>% 
  filter(time == 6, sample_condition != "aamanitin_polya") %>% 
  bind_rows(fc_dicer) %>% 
  inner_join(predictors)



# fit linear model --------------------------------------------------------

results <- mdldata %>% 
  group_by(sample_condition) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log2FC ~ PLS1 + PLS2 + utr_log + m6A + ARE + MiR_4302nd6mer + MiR430, data = .)),
    augmented = map(fit, broom::augment)
  ) %>% 
  unnest(augmented) %>% 
  select(sample_condition, log2FC, PLS1, .fitted:.std.resid) %>% 
  inner_join(mdldata, by = c("sample_condition", "log2FC", "PLS1"))


results <- results %>% 
  group_by(sample_condition) %>% 
  mutate(
    outlier = (.resid < (quantile(.resid)[2]  - 1.5 * IQR(.resid))) | (.resid > (quantile(.resid)[4]  + 1.5 * IQR(.resid)))
  ) %>% 
  ungroup()


results %>% 
  ggplot(aes(x=.fitted, y=log2FC, color=outlier)) +
  geom_point() +
  facet_grid(~sample_condition) +
  scale_color_manual(values = c("black", "purple"))



# select elements to clone  -----------------------------------------------
# a outlier shold be consisten in the 3 conditions
outliers_genes <- results %>% 
  filter(m6A == 0, MiR430 == 0, MiR_4302nd6mer == 0, ARE < 3) %>% 
  filter(outlier) %>% 
  group_by(Gene_ID) %>% 
  summarise(n_outliers = sum(outlier)) %>% 
  filter(n_outliers == 3) %>% 
  pull(Gene_ID)


results %>% 
  filter(Gene_ID %in% outliers_genes) %>% 
  select(sample_condition, log2FC, PLS1, .resid, outlier, Gene_ID) %>% 
  arrange(Gene_ID) %>% 
  group_by(Gene_ID) %>% 
  summarise(
    PLS1 = median(PLS1), # this does not matter since is the same
    median_log2FC = median(log2FC),
    median_residual = median(.resid)
  ) %>% 
  inner_join(rename(utrs, Gene_ID = ensembl_gene_id)) %>% 
  write_csv("outlier_no_mir430_no_m6a.csv")
