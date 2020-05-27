library(tidyverse)
library(ggthemes)
library(brms)

theme_set(theme_tufte(base_family = "Helvetica"))

# get tc data -------------------------------------------------------------

amanitin <- 
  read_csv("../../data/19-01-09-RNAseqProfilesFish/rna-seq-profiles/RNAseq_tidytimecourse.csv") %>% 
  filter(sample_condition == "aamanitin_polya")

## plot maternal genes higly expressed
set.seed(123)

genes_to_plot <- 
  amanitin %>% 
  filter(is_maternal, time == 2.5) %>% 
  arrange(-TPM) %>% 
  slice(2:10) %>% 
  sample_n(5) %>% 
  pull(Gene_ID)

dta <- amanitin %>% 
  filter(Gene_ID %in% genes_to_plot)


# fit a model to get confidence interval and show in the plot -------------

fit_decay_model <- function(grp) {
  # fit the bayesian linear model to estimate the decay
  brm(formula = bf(log10(TPM) ~ time), data = grp, family = gaussian(),chains = 2, cores = 2)
}

get_predictions_and_intervals <- function(grp, fit) {
  # get the predictins an
  predict(fit, newdata = grp) %>% 
    as_tibble() %>% 
    bind_cols(grp)
  
}

dta <- 
  dta %>% 
  group_by(Gene_ID) %>% 
  nest() %>% 
  mutate(fit = map(data, fit_decay_model))

preds_intervals <- 
  dta %>% 
  mutate(
    predictions = map2(data, fit, get_predictions_and_intervals)
  ) %>% 
  unnest(predictions)


# make visualziation ------------------------------------------------------

preds_intervals %>% 
  ggplot(aes(x=time, y=TPM)) +
  geom_ribbon(aes(ymin=500, ymax=TPM), alpha=1/7) +
  geom_rug(sides="l", color = "grey") +
  geom_line(aes(y = 10**Estimate)) +
  geom_line(aes(y = 10**Q2.5), linetype=2, alpha=2/5) +
  geom_line(aes(y = 10**Q97.5), linetype=2, alpha=2/5) +
  facet_grid(~Gene_ID) +
  scale_x_continuous(breaks = c(3, 5, 7), labels = c("3 hrs", "5 hrs", "7 hrs")) +
  labs(
    y = "expression TPM",
    x = "time post-fertilization"
  ) +
  theme(
    axis.ticks.y  = element_blank(),
    strip.text.x = element_text(size = 5)
  )

ggsave("figures/examples_decay.pdf", height = 2, width = 6)  

# plot histogram of stability ---------------------------------------------

stability <- read_csv("results_data/estimated_decay_rates.csv")

stability %>% 
  filter(term == "beta") %>% 
  ggplot(aes(x=estimate)) +
  geom_histogram(bins = 50, fill="grey50", color = "white") +
  geom_hline(yintercept = c(300, 600, 900, 1200),color="white", size=1/9) +
  scale_y_continuous(breaks = c(300, 600, 900, 1200)) +
  theme(axis.ticks = element_blank()) +
  labs(
    y = "frequency",
    x = "decay rate"
  )
ggsave("figures/decay_rates_hist.pdf", height = 1.5, width = 3.5)

