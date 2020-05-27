# Here, I produce a learning curve and evaluate the models
# in a train and test set
library(codonr)
library(tidyverse)
library(parallel)

genes_for_training <- unique(train_set$gene_id)

# define the number of genes that will be used for
# traing
set.seed(123)
lcurve <- tibble(
  n = seq(from=1000, to=length(genes_for_training), length.out = 20) %>% 
    as.integer()
)
# get the gene ids
lcurve <- 
  lcurve %>% 
  mutate(
    gene_ids = map(n, ~sample(genes_for_training, size = ., replace = F)),
    # get the index from the training data that contains the genes
    index = map(gene_ids, ~which(train_set$gene_id %in% .)),
    # logical indexes to retrieve the genes
    X_train = map(index, function(x) datos_preprocessed$X_train[x, ]),
    y_train = map(index, function(x) datos_preprocessed$y_train[x]),
  ) %>% 
 select(-index, -gene_ids)

# fit model to each of the data sets --------------------------------------

lcurve <- 
  lcurve %>% 
  mutate(
    model = map2(.x = X_train, .y = y_train, .f = function(x, y) train_final_model(x, y))
  )

