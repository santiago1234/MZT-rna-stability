library(caret)
library(recipes)
library(gbm)
library(tidyverse)

# get the data with the codon compositon
datos <- prepare_train_and_test_sets()

  # apply the data pre-processing
rcipe <- preprocessing(datos$X_train)

X_train <- bake(rcipe, datos$X_train)
X_test <- bake(rcipe, datos$X_test)


list(
  X_train = X_train,
  X_test = X_test,
  y_train = datos$y_train,
  y_test = datos$y_test
) %>% 
  write_rds("results_data/preprocessed_data.rds")
