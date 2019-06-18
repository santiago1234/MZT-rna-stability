# feature importance and partial dependence analysis
# interpretable machine learning
library(tidyverse)
library(iml)
library(randomForest)

# for this analysis I will take a random sample of 20K genes

dta <- read_csv("../19-04-30-EDA/results_data/training_set.csv")

set.seed(2019)
dta <- dta %>% 
  sample_n(30000) %>% 
  mutate(id_col = paste(gene_id, specie, cell_type, datatype))

# get the codon composition

cc <- codonr::codon_composition(dta, orf_col = "coding", id_col = "id_col") # my R package

# add the data together

dta <- 
  dta %>% 
  inner_join(cc) %>% 
  select(
    -coding,
    -id_col,
    -gene_id
  ) %>% 
  replace_na(list(utrlenlog = median(.$utrlenlog, na.rm = T))) %>% 
  mutate(
    specie = as.factor(specie),
    cell_type = as.factor(cell_type),
    datatype = as.factor(datatype)
  )


rf <- randomForest(decay_rate ~ ., data = dta, ntree = 50)

# feature importance ------------------------------------------------------

X <- select(dta, -decay_rate)
predictor <- Predictor$new(rf, data = X, y = dta$decay_rate)

imp <- FeatureImp$new(predictor, loss = "mae")

fimp <- plot(imp)$data %>% 
  as_tibble()

write_csv(fimp, "results_data/feature_importance.csv")

# get the effect for the top features -------------------------------------

top_features <- 
  fimp %>% 
  arrange(-importance) %>% 
  filter(str_length(feature) == 3) %>% 
  pull(feature) %>% 
  .[1:8] %>% 
  as.character()

top_features <- c(top_features, "utrlenlog", "cdslenlog")

get_ale_feature_effect <- function(feature) {
  ale <- FeatureEffect$new(predictor, feature = feature)
  efect <- ale$plot()$data
  colnames(efect) <- c(".ale", ".type", "x")
  efect$feature <- feature
  efect
}

map_df(top_features, get_ale_feature_effect) %>% 
  write_csv("results_data/feature_effects.csv")
