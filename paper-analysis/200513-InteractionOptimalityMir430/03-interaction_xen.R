# do the analysis for xenopus

library(tidyverse)
library(gam)
library(ggthemes)
library(tidybayes)
library(cowplot)

theme_set(theme_tufte(base_family = 'Helvetica'))

xen <- read_csv("results-data/mir_optimality_fc.csv") %>% 
  filter(specie == "xenopus")

# create interacting variables fish and frog

xen <- 
  xen %>% 
  mutate(
    interaction1 = PLS1 * MiR430,
    interaction2 = PLS2 * MiR430
  )


# fit model and get p value comparing the significance of the inte --------

fit <- gam(log2FC ~  PLS2  +  s(interaction2, 3) + MiR430, data = xen)
fit0 <-  gam(log2FC ~ PLS2 + MiR430, data = xen)
anova(fit0, fit, test="F")

par(mfrow=c(1,3))
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
  geom_point(data = xen, aes(y = log2FC), alpha = .7, shape = 16, color = "grey") +
  geom_line(size = 1) +
  scale_color_manual(values = c("black", "forestgreen")) +
  labs(
    x = "Codon optimality level",
    y = "log2 fold change (9h / 1h)"
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
  d_b <- sample_frac(xen, replace = T)
  fit <- gam(log2FC ~  PLS2 + s(interaction2, 3) + MiR430, data = d_b)
  
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

plot2 <- results %>% 
  ggplot(aes(x = PLS1, y = diff)) +
  stat_lineribbon() +
  scale_fill_brewer("Confidence\n", palette = 3) +
  geom_hline(yintercept = 0, color = "grey20", linetype = 2) +
  coord_cartesian(ylim = c(-1.3, .5)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    subtitle = "microRNA repression depends\n on codon optimality level",
    y = "miR-427 repression\n(counterfactual effect)",
    x = "Codon optimality level"
  ) +
  theme_tidybayes()

pdf('xen_interaction.pdf', height = 5, width = 4)
plot_grid(plot1, plot2, ncol = 1, align = 'v')
dev.off()


plot2 +
  coord_cartesian(ylim = c(-1, .5)) +
  labs(subtitle = 'p = 0.0896')
ggsave('frog-boostrap.pdf', height = 3, width = 5)
