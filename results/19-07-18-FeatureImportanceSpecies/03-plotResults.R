library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

dset <- read_csv("results_data/sample_data_for_feature_importance_analysis.csv")

# fit a linear model to get the direction of effect

directions <- lm(decay_rate ~ ., data = select(dset, -gene_id, -coding)) %>% 
  broom::tidy() %>% 
  filter(!str_detect(term, "specie"), !str_detect(term, "cell"), !str_detect(term, "data"), !str_detect(term, "Intercept")) %>% 
  mutate(
    effect = sign(estimate)
  ) %>% 
  select(term, effect) %>% 
  rename(feature = term)


importance_all <- read_csv("results_data/feature_imp_all.csv") %>% 
  inner_join(directions) %>% 
  select(feature, effect, importance) %>% 
  rename(global_imp = importance)


top_30_featutes <- arrange(importance_all, -global_imp) %>% 
  pull(feature) %>% 
  .[1:30]

importance_specie <- read_csv("results_data/feature_imp_by_specie.csv") %>% 
  inner_join(importance_all)


importance_specie %>% 
  mutate(specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))) %>% 
  filter(feature %in% top_30_featutes) %>% 
  ggplot(aes(x=importance, y=reorder(feature, global_imp), color=specie)) +
  geom_point(aes(shape = as.character(effect)), alpha=.9) +
  geom_errorbarh(aes(xmin = importance.05, xmax = importance.95), alpha=1/4, height = 0) +
  scale_color_colorblind() +
  scale_shape_manual("effect", values = c(16, 17)) +
  geom_rangeframe(alpha=1/2, sides = "b") +
  labs(
    y = NULL,
    x = "Feature Importance (loss:mae)",
    title = "Top 30 predictive features"
  ) +
  theme(legend.position = "bottom")

ggsave("feature_imp.pdf", height = 7, width = 5)
