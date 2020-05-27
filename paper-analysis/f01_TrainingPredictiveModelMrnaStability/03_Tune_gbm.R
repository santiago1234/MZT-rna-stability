library(tidyverse)
library(caret)
library(doParallel)
library(codonr) # project package
# set parallel processing
cl <- makePSOCKcluster(30)
registerDoParallel(cl)


datos_preprocessed <- read_rds("results_data/preprocessed_data.rds")
# cross-validation
# We used groupedKfold so a gene is not in two groups at the same time
cros_validation = groupKFold(train_set$gene_id, k = 10)
cros_validation <- trainControl(method = "cv", index = cros_validation)


# define the parameter grid
gbmGrid <- expand.grid(
  interaction.depth = seq(1, 7, by = 1),
  n.trees = seq(100, 1000, by = 50),
  shrinkage = c(0.01, 0.1, 0.5),
  n.minobsinnode = 10
)

set.seed(100)
gbmTune <- train(
  x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
  method = "gbm",
  tuneGrid = gbmGrid,
  metric = "Rsquared",
  trControl = cros_validation
  
)

write_rds(gbmTune, "results_data/tunning_gbm.rds")


# same test data with best model predictions 
test_set %>% 
  mutate(
    predicted = predict(gbmTune, datos_preprocessed$X_test)
  ) %>% 
  select(-coding) %>% 
  write_csv("results_data/test_data_predictions.csv")
stopCluster(cl)
