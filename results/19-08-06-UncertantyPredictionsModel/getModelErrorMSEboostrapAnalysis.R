library(tidyverse)
library(boot)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))


set.seed(123)
# estimate MSE using boostrap ---------------------------------------------
# idea from: https://www.statmethods.net/advstats/bootstrapping.html

dta <- read_csv("results_data/prediction_interavals_test_data.csv")

dta <- dta %>% 
  select(-lower_ci, -upper_ci) %>% 
  mutate(mse = (y_observed - median_prediction)**2)


get_boostrap_interval_mse <- function(dgrp, R = 2000) {
  # run boostrap to estimate MSE
  # sample mean 
  samplemean <- function(x, d) {
    return(mean(x[d]))
  }
  
  b_res <- boot(data = dgrp$mse, statistic = samplemean, R = R)
  b_res_ci <- boot.ci(b_res, type="bca")
  tribble(
    ~estimate,~ci_l,~ci_u,
    b_res$t0, b_res_ci$bca[4], b_res_ci$bca[5]
  )
}

# roon boostraping

bres <- 
  dta %>% 
  group_by(specie, cell_type, datatype) %>% 
  nest() %>% 
  mutate(boostrapres = map(data, get_boostrap_interval_mse))
  
## make plot

bres <- 
  bres %>%
  unnest(boostrapres) %>% 
  mutate(
    id = paste(cell_type, datatype, sep = " | "),
    specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))
  )

bres %>% 
  ggplot(aes(
    x = estimate,
    y = reorder(id, -estimate)
  )) +
  geom_point(shape=16, alpha=.99) +
  geom_rangeframe(sides="l", alpha=3/4) +
  geom_errorbarh(aes(xmin=ci_l, xmax=ci_u, height=0), alpha=3/4) +
  facet_grid(specie~., scales = "free_y", space = "free_y") +
  labs(
    x = "mean squared error",
    y = "specie, data type",
    title = "Gradient Boostin Regressor error"
  )

ggsave("figures/model_error_by_specie.pdf", width = 4.5, height = 3.5)

bres %>% 
  select(-data) %>% 
  write_csv("results_data/MeanSquareError_boostrap_by_specie.csv")


