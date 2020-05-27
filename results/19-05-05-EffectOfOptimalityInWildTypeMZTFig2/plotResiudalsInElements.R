library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

# define some ele

utrs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv") %>% 
  select(-coding)

m6a <- read_csv("../19-02-24-OverlapPathwaysFig3/results_data/regulatory_pathways_matrix.csv") %>% 
  filter(m6A) %>% 
  pull(Gene_ID)


# load log2-fold change and predictions -----------------------------------

genes_to_predict <- bind_rows(
  read_csv("../19-01-22-PredictiveModelMrnaStability/results_data/test_data.csv"),
  read_csv("../19-01-22-PredictiveModelMrnaStability/results_data/train_data.csv")
)

mdl <- read_rds("../19-01-22-PredictiveModelMrnaStability/results_data/trained_models/svmRModel.rds")

genes_to_predict$predictedDecay <- predict(mdl, newdata = genes_to_predict)

# subset the important columns
predicted_stability <- 
  genes_to_predict %>% 
  select(Gene_ID, predictedDecay)

# load the m6a-genes

programs <- 
  utrs %>% 
  mutate(
    `all maternal genes` = TRUE,
    ARE = str_detect(`3utr`, "TATTTATT"),
    `stable polyU` = str_detect(`3utr`, "TTTTTTT"),
    `miR-430 6-mer` = str_detect(`3utr`, "(GCACTT|AGCACT)"),
    `miR-430 7-mer` = str_detect(`3utr`, "(GCACTTA|AGCACTT)"),
    `miR-430 8-mer` = str_detect(`3utr`, "AGCACTTA"),
    m6A = ensembl_gene_id %in% m6a
  ) %>% 
  select(-`3utr`)

programs <- 
  programs %>% 
  rename(Gene_ID = ensembl_gene_id) %>% 
  gather(key = program, value = gene_in_program, -Gene_ID) %>% 
  filter(gene_in_program) %>% 
  select(-gene_in_program)

# load log fold change ----------------------------------------------------

log2fc <-
  read_csv("../../data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv") %>% 
  filter(is_maternal, time %% 1 == 0, time > 4, time < 8) %>% 
  select(-is_maternal)
  
dta <- inner_join(predicted_stability, log2fc)


# estimate the residuals of the model -------------------------------------
# each residual is estimated for each group: time, sample_condition

residuals <- 
  dta %>% 
  group_by(time, sample_condition) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(log2FC ~ predictedDecay, data = .)),
    augemented = map(fit, broom::augment)
  ) %>% 
  unnest(augemented) %>% 
  inner_join(dta, by=c("time", "sample_condition", "log2FC", "predictedDecay"))


resdata <- inner_join(programs, residuals) %>% 
  mutate(
    program = factor(
      program,
      levels = c("all maternal genes", "miR-430 8-mer", "miR-430 7-mer", "miR-430 6-mer",
                 "ARE", "stable polyU", "m6A" )
    )
  )
  

resdata %>% 
  mutate(time = paste0(time, " hrs")) %>% 
  ggplot(aes(
    x = program, y = .resid, fill = program, alpha = program
  )) +
  geom_boxplot(outlier.shape = NA, size = 1/3) +
  geom_rangeframe(color = "black") +
  scale_fill_manual(
    values = c("grey", "#56B4E9",
               "#56B4E9", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")
  ) +
  scale_alpha_manual(
    values = c(1, 1, .8, .6, 1, 1, 1)
  ) +
  coord_cartesian(ylim = c(-5.5, 5)) +
  facet_grid(time~sample_condition) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "residual\n(observed - predicted)"
  )

ggsave("figures/residuals_all_programs.pdf", height = 6, width = 5)
