library(tidyverse)
library(ggridges)

theme_set(theme_light(base_family = "Helvetica"))
results <- read_csv("results_data/results_lm_10_foldCV.csv")
results <- select(results, sample_condition, pathway, contains("test")) %>% 
  mutate(sample_condition = factor(sample_condition, levels = c("ribo0", "poly-A")))

# split into test results and mean results --------------------------------

folds <- 
  results %>% 
  select(sample_condition, pathway, contains("split")) %>% 
  gather(key = fold, value = r2, -sample_condition, -pathway)

results %>% 
  ggplot(aes(x=mean_test_score, y=pathway)) +
  geom_point(
    data = folds, aes(x=r2, y=pathway),
    position = "jitter",
    color = "grey",
    size = 1,
    alpha = 1/2
  ) +
  geom_point() +
  geom_errorbarh(aes(
    xmin = mean_test_score - 2*std_test_score,
    xmax = mean_test_score + 2*std_test_score,
    ), height=0) +
  facet_grid(~sample_condition) +
  labs(
    title = "Which pathway explains more variation in maternal mRNA stability?",
    subtitle = "Linear Model (Repeated 10-fold cross validation)",
    x = "explained variation",
    y = NULL
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))

ggsave("figures/fig2c_wich_pathway_is_the_strongest.pdf", height = 2.5, width = 6)
ggsave("figures/fig2c_wich_pathway_is_the_strongest.png", height = 2.5, width = 6)