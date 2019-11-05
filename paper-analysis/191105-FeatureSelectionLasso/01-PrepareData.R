library(codonr)
library(tidyverse)

dataCo <- prepare_train_and_test_sets()

all_data <- bind_rows(
  mutate(dataCo$X_train, stability = dataCo$y_train),
  mutate(dataCo$X_test, stability = dataCo$y_test),
) %>% 
  select(stability, everything())

write_csv(all_data, "results-data/data.csv")