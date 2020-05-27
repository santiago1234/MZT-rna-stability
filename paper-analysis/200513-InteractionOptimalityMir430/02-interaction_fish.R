# do the analysis for fish

library(tidyverse)
library(gam)
library(ggthemes)
library(tidybayes)
library(cowplot)
library(ggforce)

theme_set(theme_tufte(base_family = 'Helvetica'))

fish <- read_csv("results-data/mir_optimality_fc.csv") %>% 
  filter(specie == "fish")

# create interacting variables fish and frog

fish <- 
  fish %>% 
  mutate(
    interaction1 = PLS1 * MiR430,
    interaction2 = PLS2 * MiR430
  )


# fit model and get p value comparing the significance of the inte --------

fit <- gam(log2FC ~  PLS1 + PLS2  + s(interaction1, 5) + s(interaction2, 5) + MiR430, data = fish)
fit0 <-  gam(log2FC ~ PLS1 + PLS2 + MiR430, data = fish)
anova(fit0, fit, test="F")

par(mfrow=c(1,5))
plot(fit, se=TRUE,col="blue")


# plot data grid (counterfactual predictions) -----------------------------

x_new <- tibble( 
  PLS1 = seq(from = -5.5 ,to = 6, length.out = 30),
  PLS2 = seq(from = -5.5, to = 6, length.out = 30),
  id = 1:30
) %>% 
  crossing(tibble(MiR430 = c(0, 1.5))) %>% 
  mutate(
    interaction1 = MiR430 * PLS1,
    interaction2 = MiR430 * PLS2
  )

x_new$prediction <- predict(fit, newdata = x_new)

# plot implications -------------------------------------------------------

plot1 <- x_new %>% 
  ggplot(aes(x = PLS1, y = prediction, color = MiR430 > 0)) +
  geom_point(data = fish, aes(y = log2FC), alpha = .7, shape = 16, color = "grey") +
  geom_line(size = 1) +
  scale_color_manual(values = c("black", "forestgreen")) +
  labs(
    x = "Codon optimality level",
    y = "log2 fold change (6h / 2h)"
  ) +
  tidybayes::theme_tidybayes() +
  scale_x_continuous(expand = c(0, 0))

x_new %>% 
  select(-contains('interaction')) %>% 
  pivot_wider(names_from = MiR430, values_from = prediction) %>% 
  mutate(diff = `0` - `1.5`) %>% 
  ggplot(aes(x = PLS1, y = diff)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2)

# bootstrap analysis ------------------------------------------------------
# the pourpose of the bostrap is to provide confidence interval for the line
# estimated above

boostrap <- function() {
  d_b <- sample_frac(fish, replace = T)
  fit <- gam(log2FC ~  PLS1 + PLS2  + s(interaction1, 5) + s(interaction2, 5) + MiR430, data = d_b)
  
  x_new$prediction <- predict(fit, newdata = x_new)
  
  x_new %>% 
    select(-contains('interaction')) %>% 
    pivot_wider(names_from = MiR430, values_from = prediction) %>% 
    mutate(diff = `1.5` - `0`)
  
}


results <- tibble(
  bid = as.character(1:50)
) %>% 
  mutate(b = map(bid, function(i) boostrap())) %>% 
  unnest(b)

set.seed(42)
plot2 <- results %>% 
  ggplot(aes(x = PLS1, y = diff)) +
  stat_lineribbon() +
  scale_fill_brewer("Confidence\n", palette = 3) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  coord_cartesian(ylim = c(-1.3, .5)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    subtitle = "microRNA repression depends\n on codon optimality level",
    y = "miR-430 repression\n(counterfactual effect)",
    x = "Codon optimality level"
  ) +
  theme_tidybayes()

pdf('fish_interaction.pdf', height = 5, width = 4)
plot_grid(plot1, plot2, ncol = 1, align = 'v')
dev.off()

# boxplot -----------------------------------------------------------------

category <- function(x) {
  if (between(x, -Inf, -3)) return('non-optimal')
  if (between(x, -.5, -.3)) return('neutral')
  if (between(x, 3, Inf)) return('optimal')
  return('out')
  
}

mir_datos_boxplot <- 
  fish %>% 
  mutate(
    opt = map_chr(PLS1, category)
  ) %>% 
  filter(opt != 'out') %>% 
  mutate(opt = factor(opt, levels = c('non-optimal', 'neutral', 'optimal'))) %>% 
  select(gene_id, log2FC, opt, MiR430)


mir_datos_boxplot %>% 
  mutate(MiR430 = if_else(MiR430 > 0, 'seed', 'no seed')) %>% 
  ggplot(aes(y = log2FC, x = MiR430)) +
  coord_cartesian(ylim = c(-5, 2)) +
  geom_boxplot(aes(fill = opt), outlier.shape = NA, alpha = .9, size = .3) +
  scale_fill_manual(values = c('blue', 'grey', 'red')) +
  theme_tidybayes() +
  facet_grid(~opt) +
  labs(
    x = 'miR-430 seed',
    y = "log2 fold change\n(6h / 2h)"
  ) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 20, hjust = 1))
ggsave('fish-boxplot.pdf', height = 2, width = 3)

plot2 +
  coord_cartesian(ylim = c(-1, .5)) +
  labs(subtitle = 'p = 0.001')
ggsave('fish-boostrap.pdf', height = 3, width = 5)
