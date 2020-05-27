library(codonr)
library(tidyverse)

# sample the data (this is extremly important) ----------------------------
set.seed(123)

sample_frac(train_set, size = 1, replace = FALSE) %>% 
  mutate(cdslenlog = log(str_length(coding) + 1)) %>% 
  write_csv("training_set.csv")


sample_frac(test_set, size = 1, replace = FALSE) %>% 
  mutate(cdslenlog = log(str_length(coding) + 1)) %>% 
  write_csv("validation_set.csv")
