library(tidyverse)

d <- read_csv("results-data/predicted-log2fold-timepoints.csv")

# bring optimality to add plot colors -------------------------------------

opt <- read_csv("../191010-PredictStabilityInMZT/results-data/mzt_predictionsResidualsLog2Fc.csv") %>% 
  select(gene_id, optimality_ratio)

colors <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')

d <- 
  d %>% 
  filter(time > 4, time < 11) %>% 
  mutate(
    time = paste(time, " hrs"),
    time = factor(time, levels = paste(5:10, " hrs"))
  ) %>% 
  inner_join(opt)


d %>% 
  mutate(
    optimality = ntile(optimality_ratio, 10),
    optimality = factor(optimality, levels = rev(1:10))
  ) %>% 
  ggplot(aes(x = pred_log2FC, y = log2FC, color = optimality)) +
  geom_point(size = .1, shape = 19, alpha = .9) +
  scale_color_manual(values = colors) +
  facet_grid(time ~ specie, scales = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Predicted log2 fold-change",
    y = "Observed log2 fold-change",
    subtitle = "MZT log2 X hrs / 2hrs"
  ) +
  theme(legend.position = "none")

ggsave("figures/prediction-by-timepoint.pdf", height = 6, width = 3)
