library(tidyverse)
library(ggthemes)
library(scales)
library(ggforce)
theme_set(theme_tufte(base_family = "Helvetica"))

datum <- readxl::read_excel("/Volumes/bazzinilab/Members/previous members/Majo/FINAL_RESULTS/miR-430/Human/Final_Data/Human_final_table.xlsx")

datum <- 
  datum %>% 
  mutate(
    Optimality = map_chr(Optimality, ~if_else(. == 1, "Optimal", "Non-optimal")),
    Optimality = factor(Optimality, levels = c("Optimal", "Non-optimal")),
    miR = map_chr(miR, ~if_else(. == 1, "Seed", "no Seed")),
    UTR = str_detect(ID, "205") %>% map_chr(~if_else(., "weak seed", "strong seed")),
    Replicate = map_chr(Replicate, ~if_else(. == 1, "replicate A", "replicate B"))
  )


## scale values by the median in each grp

datum <- 
  datum %>% 
  group_by(UTR, Replicate) %>% 
  mutate(value_scaled = as.numeric(scale(value)) )


datum %>% 
  filter(UTR != "weak seed") %>% 
  mutate(
    Optimality = factor(Optimality, levels = c("Non-optimal", "Optimal")),
    miR = factor(miR) %>% fct_rev()
  ) %>% 
  ggplot(aes(x=Optimality, y=value_scaled, fill=Optimality)) +
  geom_boxplot(size=1/5, outlier.shape = NA) +
  facet_grid(miR~.) +
  coord_flip() +
  scale_fill_manual(values = c('blue', 'red')) +
  geom_rangeframe(size = 1/7, color="black", sides = "b") +
  labs(
    y = "mCherry/GFP (scaled)",
    x = NULL
  ) +
  theme(legend.position = 'none', axis.ticks.y = element_blank())
ggsave("figures/majo_res_scaled-STRONG.pdf", height = 1.5, width = 2.5)

datum %>% 
  filter(UTR == "weak seed") %>% 
  mutate(
    Optimality = factor(Optimality, levels = c("Non-optimal", "Optimal")),
    miR = factor(miR) %>% fct_rev()
  ) %>% 
  ggplot(aes(x=Optimality, y=value_scaled, fill=Optimality)) +
  geom_boxplot(size=1/5, outlier.shape = NA) +
  facet_grid(miR~.) +
  coord_flip() +
  scale_fill_manual(values = c('blue', 'red')) +
  geom_rangeframe(size = 1/7, color="black", sides = "b") +
  labs(
    y = "mCherry/GFP (scaled)",
    x = NULL
  ) +
  theme(legend.position = 'none', axis.ticks.y = element_blank())
ggsave("figures/majo_res_scaled-WEAK.pdf", height = 1.5, width = 2.5)

# get pvalues and stats ---------------------------------------------------



# Get-PVALUES
# micro-RNA effect
datum %>% 
  group_by(miR, UTR) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~t.test(value_scaled ~ Optimality, data = .)),
    tidy = map(fit, broom::tidy)
  ) %>% 
  unnest(tidy) %>% 
  select(p.value)


datum %>% 
  group_by(UTR) %>% 
  nest() %>% 
  mutate(
    fit = map(data, ~t.test(value_scaled ~ miR, data = .)),
    tidy = map(fit, broom::tidy)
  ) %>% 
  unnest(tidy) %>% 
  select(p.value)


