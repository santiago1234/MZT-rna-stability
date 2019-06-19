library(tidyverse)
library(ggthemes)
theme_set(theme_tufte(base_family = "Helvetica"))

# get the expression data -------------------------------------------------
rpe <- read_csv("../../../180815-orfome/data/Q-experimental-data/RPE_RPKM.csv")
hela <- read_csv("../../../180815-orfome/data/Q-experimental-data/hela_rpkm.csv")
k562 <- read_csv("../../../180815-orfome/data/Q-experimental-data/unormalk562ad.csv")
c293t <- read_csv("../../../180815-orfome/data/Q-experimental-data/unormal293tad.csv")

# tidy --------------------------------------------------------------------
# I only take the time zero and take the mean expression level of the replicates
# then I take the log2 of the expression

rpe <- rpe %>% 
  select(sample:s_rep.0.rep) %>% 
  mutate(expression =  (s_rep.0 + s_rep.0.rep) / 2) %>% 
  rename(gene_id = sample) %>% 
  select(gene_id, expression) %>% 
  mutate(
    expression = log2(expression),
    cell_type = "RPE"
  )

hela <- hela %>% 
  select(sample: s_hela.0.rep) %>% 
  mutate(expression =  (s_hela.0 + s_hela.0.rep) / 2) %>% 
  rename(gene_id = sample) %>% 
  select(gene_id, expression) %>% 
  mutate(
    expression = log2(expression),
    cell_type = "hela"
  )

k562 <- k562 %>% 
  select(GeneID:`k562 0h rep`) %>% 
  mutate(expression =  (`k562 0h` + `k562 0h rep`) / 2) %>% 
  rename(gene_id = GeneID) %>% 
  select(gene_id, expression) %>% 
  mutate(
    expression = log2(expression),
    cell_type = "k562"
  )

c293t <- c293t %>% 
  select(GeneID:`293t 0h rep`) %>% 
  mutate(expression =  (`293t 0h` + `293t 0h rep`) / 2) %>% 
  rename(gene_id = GeneID) %>% 
  select(gene_id, expression) %>% 
  mutate(
    expression = log2(expression),
    cell_type = "293t"
  )

rna_expression_human <- bind_rows(rpe, hela, k562, c293t)

# load the predicted stability --------------------------------------------
predicted_stability <- 
  read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-05-01-TrainModels/results_data/val_gbm.csv") %>% 
  filter(specie == "human", datatype == "endogenous")


dset <- inner_join(predicted_stability, rna_expression_human, by = c("gene_id", "cell_type")) %>% 
  select(-specie, -datatype, -mdlname) %>% 
  filter(expression > -4) # genes probably not expressed
  

dset %>% 
  gather(key = estimate, value = stability, -gene_id, -expression, -cell_type) %>% 
  ggplot(aes(x=stability, y=expression)) +
  geom_point(size=3/5, shape=16, alpha=.7) +
  facet_grid(~estimate) +
  ggpubr::stat_cor(size=2) +
  geom_rangeframe() +
  labs(
    y = "homeostasis mRNA level\nRPKM",
    x = "mRNA stability",
    subtitle = "the model improves inferences\nof mRNA levelbassed on mRNA stability"
  )

ggsave("corelation_decay_mrna_level.pdf", height  = 2, width = 3.5)

# fit poly linear model ---------------------------------------------------
fit <- lm(expression ~ cell_type + poly(predicted, 3), data = dset)

fit %>% 
  broom::augment() %>% 
  ggplot(aes(x=.fitted, y=expression)) +
  geom_point(shape=16, alpha=.75, size=3/5) +
  geom_abline(linetype=2.5) +
  ggpubr::stat_cor(size=1.5) +
  facet_grid(~cell_type, scales = "free_x") +
  scale_x_continuous(breaks = c(3, 5, 7)) +
  labs(
    title = "Codon Optimality Affects Homeostasis mRNA level in human cells",
    x = "expected mRNA level based on predicted decay",
    subtitle = "test data",
    y= "mRNA expression\nRPKM"
  )
ggsave("homeostasis_mrna_level_human.pdf", width = 7, height = 2)
