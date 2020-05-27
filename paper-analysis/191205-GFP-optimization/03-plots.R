library(tidyverse)
library(ggthemes)
library(ggforce)
library(scales)
library(optimalcodonR)

intensities <- read_csv("results-data/gfps-experiments-intensities.csv")
predictions <- read_csv("results-data/gfp_sequences.csv")

gfps <- 
  intensities %>% 
  inner_join(predictions)

gfps %>% 
  ggplot(aes(x=half_life, y=intensity, group =sequence, color=log10(intensity))) +
  geom_jitter(shape=16, alpha=.7, size=2) +
  geom_rangeframe(size=1/5, color = 'black') +
  geom_mark_rect(expand = unit("3", "mm"), size=1/4) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10**-2, 10**2.5),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = c(3, 4, 5), labels = c("3 hrs", "4 hrs", "5 hrs")) +
  scale_color_viridis_c() +
  theme_tufte() +
  theme(
    legend.position = 'none'
  ) +
  labs(
    x = "Predicted mRNA half-life",
    y = "log10 GFP/mCherry"
  )
ggsave("gfp_plot.pdf", height = 3, width = 4)


# analysis to compare pairwaise gfp intensities
# I am calling these results in the manuscript
intensity_eGFP <- filter(gfps, sequence == "eGFP") %>% pull(intensity)
intensity_neutralGFP <- filter(gfps, sequence == "neutral") %>% pull(intensity)
intensity_OptimizedGFP <- filter(gfps, sequence == "suprema") %>% pull(intensity)
intensity_DeoptimizedGFP <- filter(gfps, sequence == "infima") %>% pull(intensity)


t.test(log2(intensity_eGFP),  log2(intensity_neutralGFP))
t.test(log2(intensity_neutralGFP),  log2(intensity_DeoptimizedGFP))
t.test(log2(intensity_OptimizedGFP),  log2(intensity_eGFP))
t.test(log2(intensity_OptimizedGFP),  log2(intensity_DeoptimizedGFP))


t.test(intensity_OptimizedGFP,  intensity_DeoptimizedGFP)
