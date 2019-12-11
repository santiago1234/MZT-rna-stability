library(tidyverse)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))


results <- read_csv("results-data/lasso-coefficients.csv") %>% 
  filter(predictor != "(Intercept)") %>% 
  filter(!predictor %in% c("utrlenlog", "len_log10_coding"))


# get the average to sort the heatmap labels ------------------------------

predictors <- 
  results %>% 
  group_by(predictor, specie) %>% 
  mutate(
    avg = mean(coef)
  ) %>% 
  ungroup() %>% 
  group_by(predictor) %>% 
  summarise(avg_total = mean(avg)) %>% 
  arrange(avg_total) %>% 
  pull(predictor)


results <- results %>% 
  mutate(
    predictor = factor(predictor, levels = predictors),
    id = paste(cell_type, datatype, sep = " ")
  )


# get the median effect by codon ------------------------------------------


classify <- function(x, cutoff = .01) {
  if (abs(x) < cutoff) return("neutral")
  if (x > cutoff) return("optimal")
  if (x < cutoff) return("non-optimal")
}

effects <- results %>% 
  group_by(predictor, specie) %>% 
  summarise(mf = median(coef)) %>% 
  ungroup() %>% 
  group_by(predictor) %>% 
  group_by(predictor) %>% 
  summarise(mf = median(mf)) %>% 
  mutate(
    opt = map_chr(mf, classify),
    opt = factor(opt, levels = c("non-optimal", "neutral", "optimal"))
  )




# plot heatmap ------------------------------------------------------------

results %>% 
  mutate(
    specie = factor(specie, levels = c("xenopus", "fish", "mouse", "human"))
  ) %>% 
  filter(!predictor %in% c("TAG", "TAA", "TGA")) %>% 
  inner_join(effects) %>% 
  ggplot(aes(
    x = predictor,
    y = id,
    fill = coef
  )) +
  geom_tile() +
  scale_fill_viridis_c(option = "B",limits = c(-.08, .08), oob = scales::squish) +
  facet_grid(specie~opt, scales = "free", space = "free") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 80, hjust = 1, size = 6),
    axis.ticks.x = element_line(size=1/4),
    axis.ticks.y = element_blank()
  ) +
  labs(x= NULL, y = NULL)
ggsave("lasso-coefs.pdf", height = 2, width = 9.5)





