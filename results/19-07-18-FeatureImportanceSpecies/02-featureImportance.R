library(tidyverse)
library(iml)
library(randomForest)
library(parallel)

# in this script for each specie I run the feature importance analysis
# then I do a global feature analysis where I model all the data together
# i save the importance score to a file

dset <- read_csv("results_data/sample_data_for_feature_importance_analysis.csv")

get_mdldata <- function(specie_) {
  ## preprocess the data to fit the random forest
  # drop features with 0 variance
  
  grp <- filter(dset, specie == specie_)
  
  # drop all the categorical features
  select(grp, -specie, -cell_type, -datatype, -gene_id, -coding)
  
}

feature_importance <- function(mdldata) {
  # fit random forest with 50 trees
  rf <- randomForest(decay_rate ~ ., data = mdldata, ntree = 50)
  
  # get the importance scores
  X <- select(mdldata, -decay_rate)
  
  predictor <- Predictor$new(rf, data = X, y = mdldata$decay_rate)
  
  imp <- FeatureImp$new(predictor, loss = "mae")
  
  plot(imp)$data %>% 
    as_tibble()
  
}


# run the analysis ---------------------------------------------------------

results <- 
  select(dset, specie) %>% 
  distinct() %>% 
  mutate(
    imp = map(specie, function(x) feature_importance(get_mdldata(x)))
  )

results <- 
  results %>% 
  unnest(imp)

# fit a population model that takes into account all the vars -------------

mda_all <- select(dset, -gene_id, -coding) %>% 
  mutate(
  specie = as.factor(specie),
  cell_type = as.factor(cell_type),
  datatype = as.factor(datatype)
)

imp_all <- feature_importance(mda_all)

## save results to file

write_csv(imp_all, "results_data/feature_imp_all.csv")
write_csv(results, "results_data/feature_imp_by_specie.csv")
