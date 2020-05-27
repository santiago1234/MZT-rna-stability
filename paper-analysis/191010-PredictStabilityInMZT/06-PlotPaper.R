library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = "Helvetica"))

mztdatos <- read_csv("results-data/mzt_data_predictionsAndObservedStability.csv") %>% 
  filter(!is.infinite(log2FC), !is.na(log2FC))


# scaled the fold change --------------------------------------------------

mztdatos <- mztdatos %>% 
  group_by(specie) %>% 
  mutate(log2FC = as.numeric(scale(log2FC)))

fit <- lm(log2FC ~ specie + predicted_stability + predicted_stability:specie, data = mztdatos)

mztdatos$predicted <- predict(fit)

mztdatos$resid <- resid(fit)


mztdatos %>% 
  mutate(mir430 = str_detect(`3utr`, "GCACTT")) %>% 
  ggplot(aes(x=predicted, y=log2FC, color=mir430)) +
  geom_point(shape=16, size=1/3, alpha=.9) +
  facet_grid(~specie) +
  geom_rangeframe(color="black", size=1/5) +
  scale_color_manual(values = c("grey", "forestgreen")) +
  ggpubr::stat_cor(color="black") +
  coord_cartesian(xlim = c(-1, 1.3)) +
  theme(axis.ticks = element_line(size=1/5))

ggsave("figures/mzt-predictions.pdf", height = 2, width = 5)

mztdatos %>% 
  filter(specie == "fish") %>% 
  ggplot(aes(predicted, optimality_ratio)) +
  geom_point(shape=16, size=1/3, alpha=.9) +
  geom_rangeframe() +
  ggpubr::stat_cor(color="black") +
  theme(axis.ticks = element_line(size=1/5))

ggsave("figures/Noptimal-predictions.pdf", height = 2, width = 3)


# plot residuals against predictions --------------------------------------

mztdatos %>% 
  mutate(mir430 = str_detect(`3utr`, "GCACTT")) %>% 
  count(specie, mir430)

mztdatos %>% 
  mutate(mir430 = str_detect(`3utr`, "GCACTT")) %>% 
  ggplot(aes(x=predicted, y=resid, color=mir430)) +
  geom_point(shape=16, size=2/3, alpha=.99) +
  facet_wrap(~specie, scales = "free_x") +
  geom_hline(yintercept = 0, size=1/4) +
  geom_rangeframe(color="black", size=1/5) +
  scale_color_manual(values = c("grey70", "#009E73"))

ggsave("figures/residuals_plot.pdf", height = 3, width = 6)

