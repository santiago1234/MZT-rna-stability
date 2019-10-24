# Bootstrap 95% CI for R-Squared
library(boot)
library(tidyverse)
library(parallel)

# function to obtain R-Squared from the data
fc <- function(data, indices) {
  d <- data[indices,]
  cor(d$predicted, d$observed)**2
  
}

get_boostrap_results <- function(datos) {
  bstats <- boot(data=datos, statistic=fc, R = nrow(datos))
  bstats_ci <- boot.ci(bstats, type = "norm") 
  tibble(
    r2_estimate = bstats$t0,
    r2_ci_low = bstats_ci$normal[2],
    r2_ci_high = bstats_ci$normal[3],
    `confidence %` = bstats_ci$normal[1] * 100
  )
}


# TODO: replace with the best model results

predictions <- read_csv("results-data/test_data_predictions.csv")


# population boostrap -----------------------------------------------------

get_boostrap_results(predictions) %>% 
  write_csv("results-data//boostrap_r2_all_boostrap.csv")

predictions <- 
  predictions %>% 
  group_by(specie, cell_type, datatype) %>% 
  nest() %>% 
  mutate(
    boostrap = mclapply(data, get_boostrap_results, mc.cores = 4)
  )

predictions %>% 
  unnest(boostrap) %>% 
  mutate(
    n = map_dbl(data, nrow)
  ) %>% 
  ungroup() %>% 
  select(-data) %>% 
  write_csv("results-data//boostrap_r2_byProfile.csv")
