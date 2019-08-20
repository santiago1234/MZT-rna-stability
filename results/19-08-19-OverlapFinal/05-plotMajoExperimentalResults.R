library(tidyverse)
library(ggthemes)
library(scales)

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


# DATA PROCESSING
# 1. scale the data
# 2. get the median values

plot_data <- datum %>% 
  group_by(UTR, Replicate) %>% 
  mutate(value = (value - min(value)) / (max(value) - min(value)) ) %>% 
  ungroup() %>% 
  group_by(UTR, Replicate, miR, Optimality) %>% 
  summarise(
    med_int = median(value),
    std = sd(value)
  )

plot_data %>% 
  filter(UTR == "strong seed") %>% 
  filter(Replicate == "replicate B") %>% 
  ggplot(aes(x = miR, y = med_int, fill=Optimality)) +
  facet_grid(~Replicate) +
  geom_bar(stat="identity", position = "dodge", width=.4, alpha=2/3) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_hline(yintercept=1:4, col="white", lwd=1) +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none"
    )

  

plot_data %>% 
  ggplot(aes(x=miR, y=med_int, fill = Optimality)) +
  geom_bar(stat="identity", position = "dodge", alpha=2/3) +
  geom_hline(yintercept=c(0, .25, .5, .75 ,1), col="white", lwd=1/3) +
  facet_grid(Replicate~UTR) +
  scale_x_discrete(labels = c("no microRNA site", "microRNA site")) +
  labs(
    y = "median intensity\nmCherry/GFP (scaled)",
    x = NULL
  ) +
  scale_fill_manual(values = c("red", "blue")) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 30, hjust=1),
    legend.position = "none"
  )

ggsave("figures/majo_res.pdf", width = 3, height = 4)
