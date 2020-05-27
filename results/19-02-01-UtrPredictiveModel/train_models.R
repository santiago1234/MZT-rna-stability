library(tidyverse)
library(broom)
library(rsample)


# get train data ----------------------------------------------------------

utrs <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_3utr6mer_composition.csv") %>% 
  rename(Gene_ID = ensembl_gene_id)

data <- codonr::load_decay_aa_codon_composition_data(
  "../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_codon_composition.csv",
  "../19-01-11-GetDecayRateFromTimeCourse/results_data/estimated_decay_rates.csv"
) %>% 
  select(Gene_ID, decay_rate) %>%
  inner_join(utrs) %>% 
  select(-Gene_ID) %>% 
  select(-`3utr`)

dsplit <- initial_split(data)
train_data <- training(dsplit)
test_data <- assessment(dsplit)
  
mdl <- lm(decay_rate ~ ., data = train_data)

test_data$y_pred <- predict(mdl, test_data)

test_data$y_pred %>% hist()

test_data %>% 
  ggplot(aes(y_pred, decay_rate)) +
  geom_point()

tidy_fit <- broom::tidy(mdl)
tidy_fit %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(aes(p.value, estimate, color = estimate)) +
  geom_vline(xintercept = 0.01) +
  scale_x_log10() +
  geom_text(aes(label = term), , size = 3) +
  scale_color_viridis_c()
