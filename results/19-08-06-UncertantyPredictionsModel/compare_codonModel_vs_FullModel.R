# here i compare the prediction for the codon model
# and the ful model (full model: codon + MicroRNAs+ m6A)
library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

testdata <- read_csv('../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/val_gbm.csv')

testdata %>% 
  ggplot(aes(x=predicted, y=observed)) + 
  geom_hex(bins=50) +
  scale_fill_continuous(low="grey90", high="black") +
  labs(
    x = "predicted mRNA stability",
    y = "observed mRNA stability",
    title = "codon + 3' UTR regulation predictive model"
  ) +
  ggpubr::stat_cor()

ggsave("figures/testdataCodonModelANd3utrMode-1.pdf", height = 4, width = 4.5)


mse_codon <- read_csv("results_data/MeanSquareError_boostrap_by_specie.csv") %>% 
  mutate(model = "codon content")

mse_codon_utr <- read_csv("../19-07-18-PredictiveModelWithM6AandMicroRNAs/results_data/MeanSquareError_boostrap_by_specie.csv") %>% 
  mutate(model = "codon content + 3' UTR regulation")

bind_rows(mse_codon, mse_codon_utr) %>% 
  mutate(specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))) %>% 
  ggplot(aes(x=estimate, y=id, color = model)) +
  geom_point(shape=16, alpha = .9) +
  scale_color_manual(values = c("black", "grey70")) +
  geom_errorbarh(aes(xmin=ci_l, xmax=ci_u), height=0) +
  facet_grid(specie~., space="free_y", scales="free_y")  +
  labs(
    y = NULL,
    x = "mean squared erorr"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/model_error_by_specie_ComparingCondonModelVSFullModel.pdf", width = 5, height = 3)

