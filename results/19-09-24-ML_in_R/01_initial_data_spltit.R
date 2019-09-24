library(tidyverse)
library(caret)

# initial data split ------------------------------------------------------
# Now, I split the data in a way such that a gene is not in the trainin and
# testing data at the same time

inital_data <- 
  bind_rows(
  read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/training_set.csv"),
  read_csv("../19-04-30-PredictiveModelDecayAllSpecies/19-04-30-EDA/results_data/validation_set.csv")
) %>% 
  select(gene_id, specie, cell_type, datatype, decay_rate, utrlenlog, coding)

# the partition for the group is defined bassed on the gene id
set.seed(123)
folds <- groupKFold(group = inital_data$gene_id, k = 10)

training_data <- inital_data[folds$Fold01, ]
testing_data <- inital_data[-folds$Fold01, ]

# test that a gene id is not in the two groups
stopifnot(
  sum(testing_data$gene_id %in% training_data$gene_id) == 0
)

# save training and testing data
write_csv(training_data, "results_data/train.csv")
write_csv(testing_data, "results_data/testing.csv")
