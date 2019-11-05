library(tidyverse)
library(glmnet)
library(recipes)
library(doMC)

datos <- read_csv("results-data/data.csv")

# preprocessing -----------------------------------------------------------


preprocessing <- function(datos) {
  stopifnot(!any(c("specie", "cell_type", "datatype") %in% colnames(datos)))
  
  preprocessing_pipeline <- recipe(stability ~ ., data = datos) %>% 
    step_medianimpute(len_log10_coding, utrlenlog) %>% 
    step_rm(gene_id) %>% 
    step_spatialsign(starts_with("c_")) %>% 
    # I normalize so all the predictors are in the same scale
    step_normalize(all_predictors())
  
  
  training_rec <- prep(preprocessing_pipeline, training = datos)
  bake(training_rec, new_data = datos)
}



# cv lasso ----------------------------------------------------------------

fit_glmnet_Coefs <- function(datos_preproc) {
  registerDoMC(cores=3)
  y <- datos_preproc$stability
  X <- as.matrix(select(datos_preproc, -stability))
  
  cvfit = cv.glmnet(X, y, type.measure = "mse", nfolds = 10, alpha=0)
  
  # get the coefficients maximixing the cross validation MSE
  coef.cv.glmnet(cvfit, s="lambda.min") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "predictor") %>% 
    rename(coef = `1`)
}

# pipeline composition ----------------------------------------------------

analysis <- function(datos) {
  preprocessing(datos) %>% 
    fit_glmnet_Coefs()
    
}


# run the analysis pipeline -----------------------------------------------

results <- 
  datos %>% 
  group_by(specie, cell_type, datatype) %>% 
  nest() %>% 
  mutate(
    lasso_coefs = map(data, analysis)
  ) %>% 
  unnest(lasso_coefs) %>% 
  ungroup() %>% 
  mutate(
    predictor = str_replace(predictor, "c_", "")
  )

write_csv(results, "results-data/lasso-coefficients.csv")  

