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


# I got this pallete form here http://colorbrewer2.org/#type=diverging&scheme=RdBu&n=101
colors <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
mztdatos %>% 
  mutate(
    optimality = ntile(optimality_ratio, 10),
    optimality = factor(optimality, levels = rev(1:10))
  ) %>% 
  filter(abs(log2FC) < 5) %>% 
  ggplot(aes(x=predicted, y = log2FC, color=optimality)) +
  geom_point(shape = 16, size=.5, alpha=.9) +
  facet_grid(~specie, scales = "free_x") +
  ggpubr::stat_cor(color = "black", size=3) +
  scale_color_manual(values = colors) +
  geom_rangeframe(color = "black") +
  coord_cartesian(ylim = c(-5.5, 5.5)) +
  labs(
    color = "optimality quantile\nBazzini et al., 2016",
    y = "MZT mRNA stability\nlog2 fold change scaled",
    x = "predicted stability"
  )
ggsave("figures/predictions-mzt.pdf", height = 3, width = 7)

# save predictions and residuals ------------------------------------------

write_csv(mztdatos, "results-data/mzt_predictionsResidualsLog2Fc.csv")

