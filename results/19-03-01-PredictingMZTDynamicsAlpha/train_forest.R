library(tidyverse)
library(caret)
library(randomForest)
library(rsample)

# get time course ---------------------------------------------------------

tc <- read_csv("../19-01-11-GetDecayRateFromTimeCourse/results_data/aamanitin_time_course_tidy_only_active_genes.csv")

predictors <- read_csv("../../data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_codon_composition.csv") %>% 
  rename(Gene_ID = id_col)

data <- inner_join(tc, predictors) %>% 
  mutate(TPM = log(TPM + 0.001))

# models ------------------------------------------------------------------

dsplit <- initial_split(data)

train_data <- training(dsplit)
test_data <- assessment(dsplit)


mdl1 <- randomForest(TPM ~ time, data = train_data)

test_data$pred1 <- predict(mdl1, test_data)

test_data %>% 
  ggplot(aes(pred1, TPM, color = time)) +
  geom_point(shape = ".")

# codon composition model -------------------------------------------------

mdl2 <- randomForest(TPM ~ ., select(train_data, -Gene_ID), ntree = 100, do.trace=T)
write_rds(mdl2i, "random_forest.rds")

