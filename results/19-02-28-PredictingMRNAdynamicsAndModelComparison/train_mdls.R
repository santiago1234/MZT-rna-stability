library(tidyverse)
library(caret)
library(doParallel)

# load data ---------------------------------------------------------------

predictorsD <- read_csv("results_data/predictors.csv")

epsilon <- 0.0001
response <- read_csv("../../data/19-01-09-RNAseqProfilesFish/rna-seq-profiles/RNAseq_tidytimecourse.csv") %>% 
  select(-is_maternal) %>% 
  mutate(
    logTPM = log(TPM + epsilon)
  )


get_mdl_data <- function(condition = "wt_ribo", params = c("class_grp")) {
  ## prepares data to fit predictive model & returns the model formula
  x <- predictorsD[, c("Gene_ID", params)]
  
  mdldata <- filter(response, sample_condition == condition) %>% 
    inner_join(x)
  
  fml <- 
    str_c(
    "logTPM ~ ",
    c("time", params) %>% reduce(str_c, sep = " + ")
  ) %>% 
    as.formula()
  
  list(data = mdldata, fml = fml)

  
}

controlObject <- trainControl(method = "cv", number = 5)

train_mdl <- function(fml, data) {
  cl <- makePSOCKcluster(25)
  registerDoParallel(cl)
  set.seed(669)
  rfModel <- train(fml,
                   data = data,
                   method = "rf",
                   tuneLength = 10,
                   ntrees = 1000,
                   importance = FALSE,
                   trControl = controlObject)
  stopCluster(cl)
  return(rfModel)
}


fit_to_data <- function(params, condition, mdl_tag = "tag") {
  # fit the model to the specified data
  input <- get_mdl_data(condition, params)
  mdl <- train_mdl(input$fml, input$data)
  outname <- str_c("rf_", condition, "_", mdl_tag, ".rds")
  write_rds(mdl, str_c("results_data/trained_mdls/", outname))
}

fit_to_data(c("class_grp", "PLS1", "PLS2"), condition = "wt_ribo", "optimality")


mdls <- tribble(
  ~mdl_tag, ~params,
  "null", "class_grp",
  "optimality", c("class_grp", "PLS1", "PLS2"),
  "MiR430", c("class_grp", "MiR430"),
  "m6A", c("class_grp", "m6A"),
  "all_pathways", c("class_grp", "m6A", "MiR430", "PLS1", "PLS2")
)

conditions <- tibble(conditions = response$sample_condition %>% unique())
mdls <- crossing(conditions, mdls)

pmap(
  mdls,
  .f = function(x, y, z) paste0(x, y, z)
)
