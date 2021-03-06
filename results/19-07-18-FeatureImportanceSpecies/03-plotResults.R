library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

dset <- read_csv("results_data/sample_data_for_feature_importance_analysis.csv")

# fit a linear model to get the direction of effect
# the model is fit to each specie
directions <- 
  select(dset, -gene_id, -coding, -datatype, -cell_type) %>% 
  group_by(specie) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~lm(decay_rate ~ ., data = .)),
    tidy = map(fit, broom::tidy)
  ) %>% 
  unnest(tidy) %>% 
  mutate(effect = sign(estimate)) %>% 
  filter(term != "(Intercept)") %>%
  select(specie, term, effect) %>% 
  rename(feature = term)



importance_all <- read_csv("results_data/feature_imp_all.csv")

## top 30 most important

top_30_featutes <- arrange(importance_all, -importance) %>% 
  filter(!feature %in% c("cell_type", "datatype", "specie")) %>% 
  pull(feature) %>% 
  .[1:22]
  

importance_all <- importance_all %>% 
  select(feature, importance) %>% 
  rename(global_imp = importance)



importance_specie <- read_csv("results_data/feature_imp_by_specie.csv") %>% 
  inner_join(directions, by = c("feature", "specie")) %>% 
  inner_join(importance_all)


importance_specie %>% 
  mutate(specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))) %>% 
  filter(
    feature %in% top_30_featutes,
    !feature %in% c("cdslenlog", "utrlenlog")
    ) %>% 
  ggplot(aes(x=importance, y=reorder(feature, global_imp), color=as.character(effect))) +
  geom_point(aes(shape = specie), alpha=.9) +
  scale_color_colorblind() +
  labs(
    y = NULL,
    x = "Feature Importance (loss:mae)",
    title = "Most predictive codons"
  ) +
  theme(legend.position = "bottom")

ggsave("feature_imp.pdf", height = 5.5, width = 4.5)


# small multiples ---------------------------------------------------------


importance_specie %>% 
  mutate(specie = factor(specie, levels = c("human", "mouse", "fish", "xenopus"))) %>% 
  filter(
    feature %in% top_30_featutes,
    !feature %in% c("cdslenlog", "utrlenlog")
) %>% 
  group_by(
    specie
  ) %>% 
  ggplot(aes(x=reorder(feature, -global_imp), y=importance, color=as.character(effect))) +
  geom_point(shape=16, alpha=.99) +
  geom_errorbar(aes(ymin=importance.05, ymax=importance.95), width=0) +
  geom_rangeframe(color="black", sides="b", size=1/5) +
  scale_color_manual(values = c("grey", "black")) +
  facet_grid(specie~., scale="free_y") +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1),
    legend.position = "bottom"
) +
  labs(
    x = NULL,
    y = "Feature Importance (loss:mae)"
  )
ggsave("feature_imp_small_multiple.pdf", height = 5.5, width = 4.5)
  