library(tidyverse)
library(ggthemes)
theme_set(ggthemes::theme_tufte(base_family = "Helvetica"))

quants <- read_csv("datos-intensity-8hrs.csv") %>% 
  mutate(log2ratio = -log2(GFP / RFP))

quants %>% 
  mutate(reporter = factor(reporter, levels = c("non-optimal", "optimal"))) %>% 
  ggplot(aes(x=reporter, y=log2ratio, fill=reporter)) +
  geom_boxplot(outlier.shape = NA, size=1/5) +
  labs(
    y = "intensity log2(Cherry/GFP)"
  ) +
  geom_rangeframe() +
  scale_fill_manual(values = c("#0571b0", "#ca2027")) +
  coord_flip(ylim = c(-2.3, -0.7)) +
  theme(legend.position = "none")

t.test(log2ratio ~ reporter, data = quants)
ggsave("reporter-boxplot.pdf", width = 2.5, height = 1.3)
