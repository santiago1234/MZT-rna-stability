library(tidyverse)
library(recipes)
library(rsample)
library(caret)
library(doParallel)

# set parallel processing -------------------------------------------------

cl <- makePSOCKcluster(27)
registerDoParallel(cl)

# load training data ------------------------------------------------------

decay <- codonr::load_decay_aa_codon_composition_data(
  cc_dp = "../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_codon_composition.csv",
  decay_dp = "../19-01-11-GetDecayRateFromTimeCourse/results_data/estimated_decay_rates.csv"
)


# initial train/test split ------------------------------------------------

set.seed(10)
dsplit <- rsample::initial_split(decay)
train_data <- rsample::analysis(dsplit)
test_data <- rsample::assessment(dsplit)

write_csv(train_data, "results_data/train_data.csv")
write_csv(test_data, "results_data/test_data.csv")


# data pre-processing -----------------------------------------------------

codon_recipe <- 
  recipe(decay_rate ~ ., data = train_data)  %>%
  update_role(Gene_ID, new_role = "id variable") %>%
  step_rm(starts_with("TAG")) %>% # remove colinear stop codon
  step_spatialsign(all_predictors()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors())

saveRDS(codon_recipe, "results_data/preprocessing_recipe.rds")

# I will use repeated 10 fold cross validation to estimate error r --------

controlObject <- trainControl(
  method = "repeatedcv",
  repeats = 5,
  number = 10
)

## here i specify a function to save the trained models to R objects

mdl_saver <- function(mdl_object, mdl_name) {
  outname <- str_c("results_data/trained_models/", mdl_name, ".rds")
  saveRDS(mdl_object, outname)

}

# linear models -----------------------------------------------------------

# set.seed(669)
# linearregModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "lm",
#   trControl = controlObject
# )
# mdl_saver(linearregModel, "linearregModel")
# 
# 
# set.seed(669)
# plsModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "pls",
#   tuneLength = 15,
#   trControl = controlObject
# )
# mdl_saver(plsModel, "plsModel")
# 
# 
# enetGrid <- expand.grid(
#   .lambda = c(0, .001, .01, .1),
#   .fraction = seq(0.05, 1, length.out = 20)
# )
# set.seed(669)
# enetModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "enet",
#   tuneGrid = enetGrid,
#   trControl = controlObject
# )
# mdl_saver(enetModel, "enetModel")

# mars, nns, svms ---------------------------------------------------------

# set.seed(669)
# marsModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "earth",
#   tuneGrid = expand.grid(.degree = c(1, 2), .nprune = 2:25),
#   trControl = controlObject
# )
# mdl_saver(marsModel, "marsModel")

# set.seed(669)
# svmRModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "svmRadial",
#   tuneLength = 15,
#   trControl = controlObject
# )
# mdl_saver(svmRModel, "svmRModel")


# nnetGrid <- expand.grid(
#   .decay = c(0.001, 0.01, 0.1),
#   .size = seq(1, 27, by = 2),
#   .bag = FALSE
# )
# set.seed(669)
# nnetModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "avNNet",
#   tuneGrid = nnetGrid,
#   linout = TRUE,
#   trace = FALSE,
#   maxit = 1000,
#   trControl = controlObject
# )
# mdl_saver(nnetModel, "nnetModel")

# regression and model trees ----------------------------------------------

# set.seed(669)
# rpartModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "rpart",
#   tuneLength = 10,
#   trControl = controlObject
# )
# mdl_saver(rpartModel, "rpartModel")
# 
# 
# set.seed(669)
# ctreeModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "ctree",
#   tuneLength = 10,
#   trControl =controlObject
# )
# mdl_saver(ctreeModel, "ctreeModel")
# 
# 
# set.seed(669)
# mtModel <- train(
#   codon_recipe,
#   data = train_data,
#   method = "M5",
#   trControl = controlObject
# )
# mdl_saver(mtModel, "mtModel")

# bagged, random forest, ensembls -----------------------------------------

set.seed(669)
treebagModel <- train(
  codon_recipe,
  data = train_data,
  method = "treebag",
  trControl = controlObject
)
mdl_saver(treebagModel, "treebagModel")


set.seed(669)
rfModel <- train(
  codon_recipe,
  data = train_data,
  method = "rf",
  tuneLength = 10,
  ntrees = 1000,
  importance = TRUE,
  trControl = controlObject
)
mdl_saver(rfModel, "rfModel")


gbmGrid <- expand.grid(
  .interaction.depth = seq(1, 7, by = 2),
  .n.trees = seq(100, 1000, by = 50),
  .shrinkage = c(0.01, 0.1)
)
set.seed(669)
gbmModel <- train(
  codon_recipe,
  data = train_data,
  method = "gbm",
  tuneGrid = gbmGrid,
  verbose = FALSE,
  trControl = controlObject
)
mdl_saver(gbmModel, "gbmModel")


cubistGrid <- expand.grid(
  .comittees = c(1, 5, 10, 50, 75, 100),
  .neighbors = c(0, 1, 3, 5, 7, 9)
)
cbModel <- train(
  codon_recipe,
  data = train_data,
  method = "cubist",
  tuneGrid = cubistGrid,
  trControl = controlObject
)
mdl_saver(cbModel, "cbModel")
