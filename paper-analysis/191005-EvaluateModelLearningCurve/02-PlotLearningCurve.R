library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

learning_curve <- 
  list.files('results-data', full.names = T, pattern = "lc-*") %>% 
  map_df(read_csv)


model_by_score <- learning_curve %>% 
  filter(score == "test_scores_mean", train_sizes == max(train_sizes)) %>% 
  arrange(-r2) %>% 
  pull(model)

best_r2 <- learning_curve %>% 
  filter(test_scores_mean == max(test_scores_mean)) %>% 
  pull(test_scores_mean)


learning_curve %>% 
  mutate(model = factor(model, levels = model_by_score)) %>% 
  ggplot(aes(x=train_sizes, y=test_scores_mean)) +
  geom_hline(yintercept = best_r2, color="grey10", size=1/10) +
  geom_point(shape=16, alpha=.99, color="steelblue", size=.3) +
  geom_point(aes(y=train_scores_mean), shape=16, alpha=.99, size=.3) +
  geom_line(color="steelblue", size=1/5) +
  geom_ribbon(aes(ymin=test_scores_mean, ymax=train_scores_mean), alpha=.11) +
  geom_line(aes(y=train_scores_mean), linetype=2, size=1/5) +
  facet_grid(.~model) +
  geom_rangeframe(sides = "b", size=1/10) +
  coord_cartesian(ylim = c(0, .9)) +
  scale_x_continuous(breaks = c(5000, 27000, 50000), labels = c("5k", "27k", "50k")) +
  labs(
    title = "learning curves",
    x = "number of training data points (genes)",
    y = "R2"
  )

ggsave("figures/learning_curve.pdf", height = 2, width = 7)  
