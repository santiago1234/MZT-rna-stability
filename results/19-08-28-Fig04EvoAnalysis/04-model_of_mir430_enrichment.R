# binomial model
# see the readme file

## Modeling
#see the previous analysis: __*../19-03-20-EvolutionaryPreassure/EvoPressure2.html*__
#
#The following model will be used to evaluate our hypothesis:
#  
#  $$
#  \begin{aligned}
#\text{miR-430 sites} &\sim \text{Binomial}(n=\text{3' UTR length}, p = p_i) \\
#  logit(p_i) &= \alpha + \beta\text{mRNA stability} + \beta\text{optimality}
#  \end{aligned}
#  $$

library(tidyverse)
library(ggthemes)
library(ggforce)

theme_set(theme_tufte(base_family = "Helvetica"))

opt <- read_csv("results_data/optimality_information.csv")
foldchange <- read_csv("results_data/fold_change_late_vs_early.csv")

datum <- inner_join(opt, foldchange) %>% 
  mutate(utrlen = str_length(`3utr`)) %>% 
  filter(!is.na(`3utr`), !is.na(coding), !is.infinite(log2FC))


# fit the model to fish data ----------------------------------------------


dt2 <- datum %>% 
  mutate(
    utrlen = str_length(`3utr`)
  ) %>% 
  filter(specie == "fish")

model <- brm(data = dt2, family = binomial,
             miR430 | trials(utrlen) ~ 1 + PLS1 + PLS2 + log2FC,
             prior(normal(0, 10), class = Intercept),
             seed = 10, cores = 2, chains = 2)



# predict some fact data to show the model surface ------------------------

topredict <- tibble(
  PLS1 = seq(-6, 6, length.out = 50),
  PLS2 = seq(-6, 6, length.out = 50)
  
) %>% 
  crossing(
    tibble(log2FC = seq(-6, 6, length.out = 50))
  ) %>% 
  mutate(utrlen = quantile(dt2$utrlen)[3])


topredict <- 
  fitted(model, newdata = topredict) %>% 
  as_tibble() %>% 
  bind_cols(topredict)


# higligh smarca2 and show the predictive surface -------------------------


smarca2 <- dt2 %>% 
  filter(gene_id == "ENSDARG00000008904") %>% 
  mutate(Estimate=10) %>% 
  mutate(name = "smarca2")


topredict %>% 
  ggplot(aes(x=log2FC, y=PLS1, color=Estimate)) +
  geom_point(shape=15, size=4) +
  geom_text(data = smarca2, aes(x = log2FC, y = PLS1, label=name), color="white", size=2) +
  geom_point(data = smarca2, aes(x = log2FC, y = PLS1), color="white", shape=16, alpha=.9) +
  scale_color_viridis_c(option = "A") +
  scale_x_continuous(breaks = c(-5, 0, 5)) +
  coord_cartesian(xlim = c(-6, 6), ylim=c(-6, 6)) +
  labs(
    x = "mRNA stability\n(log2 fold change)",
    y = "optimality level"
  )
ggsave("figures/04-miR430enrichmentAsfunctionOfOptimality.pdf", height = 2, width = 3)
