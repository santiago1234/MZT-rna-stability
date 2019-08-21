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


datum %>% 
  ggplot(aes(x=miR, y=value, color=Optimality)) +
  scale_y_log10() +
  geom_sina(shape=16, size=1/2, alpha=.9) +
  scale_color_manual(values = c("red", "blue")) +
  geom_rug(sides="l") +
  geom_rangeframe(size=1/5, color="black") +
  scale_x_discrete(labels = c("no microRNA site", "microRNA site")) +
  facet_grid(Replicate~UTR, scales="free_y") +
  theme(legend.position = "none") +
  labs(
    y = "log10 intensity\nmCherry/GFP",
    x = NULL
  )
ggsave("figures/majo_res.pdf", width = 3.5, height = 3)
