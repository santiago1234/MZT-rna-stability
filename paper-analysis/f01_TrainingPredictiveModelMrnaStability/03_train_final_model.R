library(tidyverse)

# this will use the whole
# training set for training
final_model <- train_final_model()
readr::write_rds(final_model, "results/gbm_final.rds")
